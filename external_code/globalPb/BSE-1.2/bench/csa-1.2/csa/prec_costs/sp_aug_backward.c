#include	<math.h>
#include	"csa_types.h"
#include	"csa_defs.h"

extern	rhs_ptr	*bucket, deq_list();
extern	char	*st_pop();
extern	long	num_buckets;
extern	void	st_reset();
extern	rhs_ptr	head_rhs_node, tail_rhs_node;
extern	lhs_ptr	head_lhs_node, tail_lhs_node;
extern	unsigned	total_e;
extern	double	epsilon;
extern	ACTIVE_TYPE	active;

extern	unsigned	sp_augs, a_scans;
extern	unsigned	sp_aug_time, myclock();

lhs_ptr	closest_node;
unsigned	long	closest_dist;
rhs_ptr	scanned;
unsigned	long	level;	/* level currently being scanned */

/*
augment() moves a unit of excess from v along an augmenting path to
cancel a deficit.
*/

void	augment(v)

lhs_ptr	v;

{
lhs_ptr	x;
rhs_ptr	w;

#ifdef	DEBUG
(void) printf("augment(%ld):\n", v - head_lhs_node + 1);
#endif

do
  {
  w = v->aug_path->head;
#if	defined(DEBUG) || defined(LOG_PATHS)
  (void) printf("%ld %ld, ", v - head_lhs_node + 1,
		w - head_rhs_node + 1);
#endif
  v->matched = v->aug_path;
  x = w->matched;
  w->matched = v;
  v = x;
  }
while (v);
#if	defined(DEBUG) || defined(LOG_PATHS)
putchar('\n');
#endif
}

/*
Doing an a_scan on w updates the current estimate of required price
changes on nodes adjacent (in the rhs sense) to w to establish deficit
reachability in the admissible graph for all excesses.
*/

void	a_scan(w)

rhs_ptr	w;

{
register	rl_aptr	b, b_stop;
register	lr_aptr	a;
register	rhs_ptr	u;
lhs_ptr		v;
register	long	wk, uk;
register	double	p;
double	u_to_w_cost;

#ifdef	DEBUG
(void) printf("doing a_scan(%ld) key=%ld\n", w - head_rhs_node + 1, w->key);
#endif

a_scans++;
b_stop = (w+1)->priced_out;
p = w->p;
wk = w->key;
for (b = w->back_arcs; b != b_stop; b++)
  if (a = b->tail->matched)
    {
    if (((u = a->head) != w) && u->node_info.priced_in && (u->key > level))
      {
      u_to_w_cost = u->p - a->rev->c + b->c - p;
      if (u_to_w_cost < 0.0)
	uk = wk;
      else
#if	defined(STRONG_PO) || !defined(USE_PRICE_OUT)
	{
	/*
	It could happen, in the case of strong price-outs, that
	priced-in arcs cause violation of the condition that node
	price changes are bounded in each iteration. In this case, or
	simply when arc costs are very widely distributed and no
	price-outs are used, some priced-in arcs can have truly huge
	costs, causing the following computation of uk to overflow.
	When this happens, ignore the key. There will be a smaller
	one for the same node.
	*/
	if (epsilon * (u->key - wk) > u_to_w_cost)
#endif
	/*
	No ceiling in the following line; such an operation could
	make the gap between a's reduced cost and that of the matching
	arc greater than epsilon, thus violating epsilon optimality.
	*/
	uk = wk + 1 + (long) (u_to_w_cost / epsilon);
#if	defined(STRONG_PO) || !defined(USE_PRICE_OUT)
	else uk = u->key;
#endif
	}
      if (u->key > uk)
	{
	if (u->key != num_buckets)
	  delete_list(u, &bucket[u->key]);
	u->key = uk;
	insert_list(u, &bucket[uk]);
	/* Keep track of to-be-admissible path through this node */
	b->tail->aug_path = b->rev;
	}
      }
    }
  else
    /*
    Encountered an excess -- b's tail isn't matched. Determine what
    price decrease on b's tail would be required to make the edge to w
    admissible. Recall that back arcs costs are offset so preferred
    arcs have zero partial reduced cost at this point, so we need only
    examine the stored cost of the present edge, rather than compute
    the minimum.
    */
    {
#if	defined(STRONG_PO) || !defined(USE_PRICE_OUT)
    if (epsilon * (closest_dist - wk) > b->c - p)
#endif
    uk = wk + (long) ceil((b->c - p) / epsilon);
#if	defined(STRONG_PO) || !defined(USE_PRICE_OUT)
    else
      uk = closest_dist;
#endif
    if (uk < closest_dist)
      {
      (v = b->tail)->aug_path = b->rev;
      closest_dist = uk;
      closest_node = v;
      if (uk == 0)
	{
#ifdef	DEBUG
	(void) printf("claiming excess at node %ld\n",
		      v - head_lhs_node + 1);
#endif
	break;
	}
      }
    }

insert_list(w, &scanned);
}

void	sp_aug()

{
rhs_ptr	w;
lhs_ptr	v;
double	delta_c, this_cost;
lr_aptr	a, a_stop;
lhs_ptr	save_active[EXCESS_THRESH];	/* Fix this. It should be */
					/* malloc'ed once, possibly at */
					/* full size, and set up by */
					/* the initialization */
					/* routines. */
lhs_ptr	*save_top = save_active,
	*active_node;
void	exit();

sp_aug_time -= myclock();
sp_augs++;

#ifdef	DEBUG
(void) printf("Doing sp_aug(): epsilon = %lg, total_e = %lu\n",
	      epsilon, total_e);
for (level = 0; level < num_buckets; level++)
  if (bucket[level] != tail_rhs_node)
     {
     (void) printf("Bucket init failure!\n");
     exit(1);
     }
#endif

for (level = 0; level < total_e; level++)
  {
  get_active_node(v);
  *(save_top++) = v;
  }

scanned = tail_rhs_node;
while (total_e > 0)
  {
  for (active_node = save_active; active_node != save_top; active_node++)
    if (!(v = *active_node)->matched)
      {
      v = *active_node;
      a_stop = (v+1)->priced_out;
#ifdef	DEBUG
      (void) printf("excess at node %ld\n", v - head_lhs_node + 1);
#endif
      a = v->first;
#ifdef	NO_FEAS_PROMISE
      if (a == a_stop)
	{
	(void) printf("Infeasible problem\n");
	exit(9);
	}
#endif
      delta_c = a->rev->c - a->head->p;
      for (a++; a != a_stop; a++)
	if ((this_cost = a->rev->c - a->head->p) < delta_c)
	  delta_c = this_cost;
#ifdef	STRONG_PO
      a_stop = v->priced_out - 1;
#else
      a_stop = v->first - 1;
#endif
      for (a--; a != a_stop; a--)
	a->rev->c -= delta_c;
      }

  /*
  Fix the following so that keys are initialized to the right thing,
  and left there when we're done. Never make the code go through all
  the nodes here. This may mean keeping track of a global list of
  nodes with deficits.
  */
  for (w = head_rhs_node; w != tail_rhs_node; w++)
    if (w->matched)
      w->key = num_buckets;
    else
      {
      w->key = 0;
      insert_list(w, &bucket[0]);
      }

  level = 0;
  closest_dist = num_buckets;
  closest_node = tail_lhs_node;

  while (level < closest_dist)
    if (bucket[level] == tail_rhs_node)
      level++;
    else
      {
      w = deq_list(&bucket[level]);
      a_scan(w);
      }

#ifdef	DEBUG
  (void) printf("level=%ld, num_buckets=%ld,\n",
		level, num_buckets);
  (void) printf("closest_dist=%ld, closest_node=%ld\n", closest_dist,
		closest_node - head_lhs_node + 1);
#endif
  /* Augment from the node with smallest required price change. */
  if (closest_node == tail_lhs_node)
    {
    (void) printf("Error: scanning failure.\n");
    exit(1);
    }
  augment(closest_node);
  while (scanned != tail_rhs_node)
    {
    w = deq_list(&scanned);
    w->p += epsilon * (closest_dist - w->key);
    w->key = num_buckets;
    }
  for (level = closest_dist; level != num_buckets; level++)
    bucket[level] = tail_rhs_node;

  total_e--;
  }

#ifdef	DEBUG
for (level = 0; level < num_buckets; level++)
  if (bucket[level] != tail_rhs_node)
     {
     (void) printf("Bucket check failure!\n");
     exit(1);
     }
#endif

sp_aug_time += myclock();
}
