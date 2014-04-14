#include	<stdio.h>
#include	"csa_types.h"
#include	"csa_defs.h"

#ifdef	USE_P_UPDATE
extern	WORK_TYPE	upd_work_thresh;
extern	void		p_update();
#endif
#ifdef	STRONG_PO
extern	WORK_TYPE	po_work_thresh;
extern	int		check_po_arcs();
#endif

extern	unsigned	double_pushes, pushes, relabelings, refines;
extern	unsigned	total_e;
extern	double		po_cost_thresh;
extern	lhs_ptr		head_lhs_node, tail_lhs_node;
extern	ACTIVE_TYPE	active;

extern	char		*st_pop(), *deq();

extern	unsigned	refine_time;

/*
All costs and prices in units of epsilon.
*/

/* Assume v has excess (is unassigned) and do a double push from v. */

void	double_push(v)

lhs_ptr	v;

{
long	mini, adm_gap, red_cost;
lr_aptr	a, a_stop, adm;
rhs_ptr	w;
lhs_ptr	u;

extern	rhs_ptr	head_rhs_node;
#ifdef	DEBUG
(void) printf("double_push(%ld) ", v - head_lhs_node + 1);
#endif
/* Find an admissible (minimum residual) arc out of v. */
a_stop = (v+1)->priced_out;
a = v->first;
#ifdef	NO_FEAS_PROMISE
if (a == a_stop)
  {
  (void) printf("Infeasible problem\n");
  exit(9);
  }
#endif
mini = a->c - a->head->p;
adm = a;
adm_gap = (long) po_cost_thresh + 1;
/*
After this loop, mini is the minimum reduced cost of an edge out of v,
and adm_gap is the difference between the second-to-minimum and the
minimum.
*/
for (a++; a != a_stop; a++)
  if (mini > (red_cost = a->c - a->head->p))
    {
    adm_gap = mini - red_cost;
    mini = red_cost;
    adm = a;
    }
  else if (adm_gap > red_cost - mini)
    adm_gap = red_cost - mini;
/* Match v through adm */
w = adm->head;
if (u = w->matched)
  /*
  If w's matched arc is priced in, go ahead and unmatch (u, w) and
  match (v, w). If w's matched arc is priced out, abort the double
  push and relabel w so v no longer prefers w.
  */
  if (w->node_info.priced_in)
    {
    pushes += 2;
    double_pushes++;
#ifdef	DEBUG
    (void) printf("matching (%ld, %ld).\n", v - head_lhs_node + 1,
		  w - head_rhs_node + 1);
    (void) printf("unmatching (%ld, %ld), ", u - head_lhs_node + 1,
		  w - head_rhs_node + 1);
#endif
    u->matched = NULL;
    make_active(u);
    v->matched = adm;
    w->matched = v;
    }
  else
    {
#ifdef	DEBUG
    (void) printf("aborting double_push\n");
#endif
    adm_gap = (long) po_cost_thresh;
    make_active(v);
    }
else
  {
  total_e--;
  pushes++;
#ifdef	DEBUG
  (void) printf("matching (%ld, %ld).\n", v - head_lhs_node + 1,
		w - head_rhs_node + 1);
#endif
  v->matched = adm;
  w->matched = v;
  }
/*
Relabel w: v's implicit price is chosen to make the implicit reduced
cost of v's new preferred arc (mini + adm_gap) equal to zero. Then w's
price is chosen so that the arc just matched has implicit reduced cost
-1.
*/
relabelings++;
w->p -= adm_gap + 1;
}

/*
void	check_matching()

{
lhs_ptr	v;
lr_aptr	a;
extern	rhs_ptr	head_rhs_node;

for (v = head_lhs_node; v != tail_lhs_node; v++)
  if (v->matched)
    {
    if (v->matched->head->matched != v)
      {
      (void) printf("Inconsistent matching. (%ld, %ld, ",
		    v - head_lhs_node + 1,
		    v->matched->head - head_rhs_node + 1);
      if (v->matched->head->matched)
        (void) printf("%ld)\n", v->matched->head->matched - head_lhs_node + 1);
      else
	(void) printf("NULL)\n");
      }
    if ((v->matched->head->node_info.priced_in &&
	 (((long) v->matched - (long) v->first < 0) ||
	  ((long) (v+1)->priced_out - (long) v->matched <= 0))) ||
	(!v->matched->head->node_info.priced_in &&
	 (((long) v->matched - (long) v->priced_out < 0) ||
	  ((long) v->first - (long) v->matched <= 0))))
      (void) printf("Inconsistent matched arc address.\n");
    }
}
*/

void	refine()

{
lhs_ptr	v;
unsigned	myclock();
#ifdef	USE_P_UPDATE
WORK_TYPE	old_refine_work_upd;
#endif
#ifdef	STRONG_PO
WORK_TYPE	old_refine_work_po;
#endif

refine_time -= myclock();
refines++;
/*
Saturate all matching arcs. The negative arcs are a subset of these,
and it seems faster to saturate them all than to infer lhs node prices
and saturate only the negative ones.
*/
total_e = 0;
for (v = head_lhs_node; v != tail_lhs_node; v++)
  {
  if (v->matched && v->matched->head->node_info.priced_in)
    {
    v->matched->head->matched = NULL;
    v->matched = NULL;
    }
  if (!v->matched)
    {
    total_e++;
    make_active(v);
    }
  }

#ifdef	USE_P_UPDATE
old_refine_work_upd = REFINE_WORK;
#endif
#ifdef	STRONG_PO
old_refine_work_po = REFINE_WORK;
#endif

#ifdef	STRONG_PO
while ((total_e > 0) || (old_refine_work_po = REFINE_WORK,
			 !check_po_arcs()))
#else
while (total_e > 0)
#endif
  {
#ifdef	USE_P_UPDATE
  if (REFINE_WORK - old_refine_work_upd > upd_work_thresh)
    {
    old_refine_work_upd = REFINE_WORK;
    p_update();
    }
#endif
#ifdef	STRONG_PO
  if (REFINE_WORK - old_refine_work_po > po_work_thresh)
    {
    old_refine_work_po = REFINE_WORK;
    (void) check_po_arcs();
    }
#endif
  get_active_node(v);
  double_push(v);
  }

refine_time += myclock();
}
