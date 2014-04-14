#include	<math.h>
#include	"csa_defs.h"
#include	"csa_types.h"

extern	rhs_ptr	head_rhs_node, tail_rhs_node;
#ifdef	LOG_PATHS
extern	lhs_ptr	head_lhs_node;
extern	rhs_ptr	head_rhs_node;
#endif

#ifdef	CHECK_EPS_OPT
extern	void		check_e_o();
extern	double	scale_factor;
#endif

extern	double	epsilon;
extern	unsigned	total_e;
extern	ACTIVE_TYPE	active;

extern	rhs_ptr	*bucket;
extern	long	num_buckets;

extern	char	*st_pop();
extern	rhs_ptr deq_list();

extern	unsigned	myclock();
extern	unsigned	sp_augs, a_scans, sp_aug_time;

void	augment(w)

rhs_ptr	w;

{
rhs_ptr	x;
lhs_ptr	v;
lr_aptr	a;

do
  {
  v = w->aug_path->tail;
#ifdef	LOG_PATHS
  (void) printf("%ld %ld", v - head_lhs_node + 1,
		w - head_rhs_node + 1);
#endif
  w->matched = v;
  a = v->matched;
  v->matched = w->aug_path;
  if (a)
    {
#ifdef	LOG_PATHS
    (void) printf(", ");
#endif
    x = a->head;
    w = x;
    }
  }
while (a);
#ifdef	LOG_PATHS
putchar('\n');
#endif
}

void	a_scan(w)

rhs_ptr	w;

{
lhs_ptr	v = w->matched;
rhs_ptr	u;
lr_aptr	a, a_stop;
double	delta_c = w->p - v->matched->c,
	w_to_u_cost;
long	wk = w->key,
	uk;

a_scans++;

#ifdef	DEBUG
(void) printf("a_scan(%ld) key=%lu\n", w - head_rhs_node + 1, w->key);
#endif

a_stop = (v+1)->priced_out;
for (a = v->first; a != a_stop; a++)
  if (a != v->matched)
    {
    u = a->head;
    w_to_u_cost = delta_c + a->c - u->p;
    if (w_to_u_cost < 0.0)
      uk = wk;
    else
      if (epsilon * (u->key - wk) > w_to_u_cost)
	uk = wk + 1 + (long) (w_to_u_cost / epsilon);
      else
	uk = u->key;
    if (u->key > uk)
      {
      if (u->key != num_buckets)
	delete_list(u, &bucket[u->key]);
      u->key = uk;
      insert_list(u, &bucket[uk]);
      u->aug_path = a;
      }
    }
}

/*
We assume that on entry to sp_aug(), all rhs nodes have their key
fields set to num_buckets.
*/

void	sp_aug()

{
lhs_ptr	v;
rhs_ptr	w;
lr_aptr	a, a_stop;
unsigned	long	level;
long	k;
double	delta_c, this_cost;
void	exit();

sp_aug_time -= myclock();
sp_augs++;

#ifdef	DEBUG
(void) printf("Doing sp_aug(): epsilon = %lg, total_e = %lu\n",
	      epsilon, total_e);
#endif

#ifdef	CHECK_EPS_OPT
  check_e_o(epsilon);
#endif

while (total_e > 0)
  {
  /*
  Get neighbors of this active node into the proper buckets.
  */
  get_active_node(v);
  a_stop = (v+1)->priced_out;
  a = v->first;
  delta_c = a->c - a->head->p;
  for (a++; a != a_stop; a++)
    if ((this_cost = a->c - a->head->p) < delta_c)
      delta_c = this_cost;
  a_stop = v->first - 1;
  for (a--; a != a_stop; a--)
    {
    /*
    Insert a's head into the proper bucket with the right key.
    */
    w = a->head;
    w->aug_path = a;
    this_cost = a->c - w->p - delta_c;
    if ((this_cost /= epsilon) < (double) num_buckets)
      {
      k = (long) this_cost;
      if (!w->matched && (k == 0))
	{
	augment(w);
	total_e--;
	break;
	}
      else if (k < num_buckets)
	{
	/*
	Here we make the (very reasonable) assumption that there are
	no multiple arcs.
	*/
	w->key = k;
	insert_list(w, &bucket[k]);
	}
      }
    }

  level = 0;
  if (a == a_stop)	/* If we didn't find a deficit and augment already */
    {
    while (level < num_buckets)
      if (bucket[level] == tail_rhs_node)
	level++;
      else
	{
	w = deq_list(&bucket[level]);
	if (w->matched)
	  a_scan(w);
	else
	  {
	  augment(w);
	  w->key = num_buckets;
	  total_e--;
	  break;
	  }
	}

    if (level == num_buckets)
      {
      (void) printf("Error: scanning failure\n");
      exit(-1);
      }
    }

  /*
  Adjust prices and clean out remaining buckets.
  */
#ifdef	DEBUG
  (void) printf("level = %lu\n", level);
#endif
  /*
  Nodes at this level and higher need no price adjustment; empty this
  level's bucket pronto.
  */
  bucket[level] = tail_rhs_node;
  /*
  Nodes at levels lower than this need price adjustment; others simply
  need to be deleted from their buckets. This is quicker than going
  through all the buckets in order and cleaning them out.
  */
  for (w = head_rhs_node; w != tail_rhs_node; w++)
    if ((k = w->key) != num_buckets)
      {
      if (k < level)
	{
#ifdef	DEBUG
	(void) printf("%ld->p -= %ld * epsilon\n", w - head_rhs_node + 1,
		      level - k);
#endif
	w->p -= (level - k) * epsilon;
	}
      else if (k > level)
	delete_list(w, &bucket[k]);
      w->key = num_buckets;
      }
#ifdef	CHECK_EPS_OPT
  check_e_o(epsilon);
#endif
  }

sp_aug_time += myclock();
}
