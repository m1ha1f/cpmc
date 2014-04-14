#include	<math.h>
#include	"csa_types.h"
#include	"csa_defs.h"

extern	rhs_ptr	*bucket, deq_list();
extern	long	num_buckets;
extern	void	st_reset();
extern	rhs_ptr	head_rhs_node, tail_rhs_node;
extern	lhs_ptr	head_lhs_node, tail_lhs_node;
extern	unsigned	total_e;
extern	double	epsilon;

extern	unsigned	p_updates, u_scans;
extern	unsigned	p_update_time, myclock();

/*
Doing a u_scan on w updates the current estimate of required price
changes on nodes adjacent (in the rhs sense) to w to establish deficit
reachability in the admissible graph for all excesses.
*/

unsigned	u_scan(w)

rhs_ptr	w;

{
register	rl_aptr	b, b_stop;
register	lr_aptr	a;
register	rhs_ptr	u;
lhs_ptr		v;
register	long	wk, uk;
register	double	p;
double	u_to_w_cost;
unsigned	excess_found = 0;

u_scans++;
b_stop = (w+1)->priced_out;
p = w->p;
wk = w->key;
for (b = w->back_arcs; b != b_stop; b++)
  if (a = b->tail->matched)
    {
    if (((u = a->head) != w) && u->node_info.priced_in)
      {
#ifdef	P_U_ZERO_BACK_MCH_ARCS
      u_to_w_cost = u->p + b->c - p;
#else
      u_to_w_cost = u->p - a->rev->c + b->c - p;
#endif
      if (u->key >= 0)
	{
	if (u_to_w_cost < 0.0)
	  uk = wk;
	else
	  /*
	  Preliminary check to make sure we're in the ballpark to
	  avoid overflow and to avoid costly double-to-long casts if
	  we don't need them.
	  */
	  if (epsilon * (u->key - wk) > u_to_w_cost)
	    uk = wk + 1 + (long) (u_to_w_cost / epsilon);
	  else
	    uk = u->key;
	if (u->key > uk)
	  {
	  if (u->key != num_buckets)
	    delete_list(u, &bucket[u->key]);
	  u->key = uk;
	  insert_list(u, &bucket[uk]);
	  }
	}
      }
    }
  else
    /*
    Encountered an excess -- b's tail isn't matched. Determine what
    price decrease on b's tail would be required to make the edge to w
    admissible. Recall that back arc costs are offset so preferred
    arcs have zero partial reduced cost at this point, so we need only
    examine the stored cost of the present edge, rather than compute
    the minimum. Avoid costly ceiling and cast when possible; also
    avoid overflows.
    */
    if ((u_to_w_cost = b->c - p) < epsilon * ((v = b->tail)->delta_reqd - wk))
      {
      uk = wk + (long) ceil(u_to_w_cost / epsilon);
      if (uk < v->delta_reqd)
	{
	if (uk == 0)
	  {
	  excess_found++;
#ifdef	DEBUG
	  (void) printf("claiming excess at node %ld\n",
			v - head_lhs_node + 1);
#endif
	  }
	v->delta_reqd = uk;
	}
      }

w->p -= epsilon * w->key;
w->key = -1;
return(excess_found);
}

void	p_update()

{
rhs_ptr	w;
lhs_ptr	v;
double	delta_c, this_cost;
long	balance, level;
lr_aptr	a, a_stop;

p_update_time -= myclock();
p_updates++;

#ifdef	DEBUG
(void) printf("Doing p_update(): epsilon = %lg, total_e = %lu\n",
	      epsilon, total_e);
#endif

for (v = head_lhs_node; v != tail_lhs_node; v++)
  {
  a_stop = (v+1)->priced_out;
#ifdef	P_U_ZERO_BACK_MCH_ARCS
  if (v->matched)
    {
    if (v->matched->head->node_info.priced_in)
      {
      /*
      Offset back arc costs so back matching arc has zero stored cost
      */
      delta_c = v->matched->rev->c;
#ifdef	STRONG_PO
      /*
      In the case of strong price-outs, we could price in a back arc
      later that's priced out now. So offset the costs of all the
      incident back arcs.
      */
      a = v->priced_out;
#else
      a = v->first;
#endif
      for (; a != a_stop; a++)
	a->rev->c -= delta_c;
      }
    }
  else
#else
  if (!v->matched)
#endif
    {
#ifdef	DEBUG
    (void) printf("excess at node %ld\n", v - head_lhs_node + 1);
#endif
    v->delta_reqd = num_buckets;
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
  }

for (w = head_rhs_node; w != tail_rhs_node; w++)
  if (w->matched)
    w->key = num_buckets;
  else
    {
    w->key = 0;
    insert_list(w, &bucket[0]);
    }

balance = -total_e;
level = 0;

while ((balance < 0) && (level < num_buckets))
  if (bucket[level] == tail_rhs_node)
    level++;
  else
    {
    w = deq_list(&bucket[level]);
    balance += u_scan(w);
    }

/*
Now figure out by how much we need to decrease prices of nodes that
didn't get scanned.
*/
for (v = head_lhs_node; v != tail_lhs_node; v++)
  {
  if (!v->matched)
    {
    if (v->delta_reqd == num_buckets)
      (void) printf("%u : excess at node %ld unclaimed after scans!\n",
		    p_updates, v - head_lhs_node + 1);
    if (v->delta_reqd > level)
      level = v->delta_reqd;
    }
  }

delta_c = level * epsilon;
for (w = head_rhs_node; w != tail_rhs_node; w++)
  {
  if ((w->key != num_buckets) && (w->key >= 0))
    delete_list(w, &bucket[w->key]);
  if (w->key >= 0)
    w->p -= delta_c;
  }

p_update_time += myclock();
}
