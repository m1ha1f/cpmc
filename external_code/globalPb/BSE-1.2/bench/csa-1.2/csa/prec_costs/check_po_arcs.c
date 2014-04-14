#include	<stdio.h>
#include	"csa_types.h"
#include	"csa_defs.h"

extern	unsigned	total_e;
extern	ACTIVE_TYPE	active;
extern	unsigned	fix_ins;
extern	lhs_ptr		head_lhs_node, tail_lhs_node;
#ifdef	DEBUG
extern	rhs_ptr		head_rhs_node;
#endif
extern	double		po_cost_thresh, epsilon;
#ifdef	QUICK_MIN
extern	void	best_build();
#endif

int	check_po_arcs()

{
lhs_ptr	v;
lr_aptr	a, a_start, a_stop;
int	one_priced_in, fix_this_node, fix_in = FALSE;
double	match_rc, this_cost, v_price, this_price, po_cutoff, thresh,
	fix_in_thresh;
#ifdef	QUICK_MIN
int	need_best_rebuild;
#endif

#ifdef	DEBUG
(void) printf("Checking priced-out arcs. total_e=%lu\n", total_e);
#endif

po_cutoff = po_cost_thresh * epsilon;
for (v = head_lhs_node; v != tail_lhs_node; v++)
  {
#ifdef	QUICK_MIN
  need_best_rebuild = FALSE;
#endif
  /*
  All routines that incorporate prices into stored costs must update
  stored costs of priced-out arcs so the following code correctly
  computes reduced costs of priced-out arcs. At the present time,
  there are no such routines.
  */
  a_stop = (v+1)->priced_out;
  if (a = v->matched)
    {
    /*
    Node v is matched. Price in any arcs not far costlier than the
    matching arc, and if there are any such arcs, make sure the
    matching arc is priced in, too.
    */
    a_start = v->first;
    one_priced_in = (a_start != a_stop);
    fix_this_node = FALSE;
    match_rc = a->c - a->head->p;
    thresh = match_rc + po_cutoff;
    fix_in_thresh = match_rc - epsilon;
    for (a = v->priced_out; a != v->first; a++)
      if ((a != v->matched) && ((this_cost = a->c - a->head->p) < thresh))
	{
	price_in_unm_arc(v, a);
	one_priced_in = TRUE;
	if (this_cost < fix_in_thresh)
	  {
	  /*
	  Epsilon-optimality violated by priced-in arc.
	  */
	  fix_in = TRUE;
	  fix_this_node = TRUE;
#ifdef	DEBUG
	  (void) printf("Fixing in arc (%ld, %ld)\n", v - head_lhs_node + 1,
			a->head - head_rhs_node + 1);
#endif
	  }
#ifdef	QUICK_MIN
	need_best_rebuild = TRUE;
#endif
	/*
	If we priced in the last priced-out arc in the list,
	a == v->first, and we need to keep a from advancing too far.
	*/
	if (a == v->first) break;
	}
    /*
    Now if matching arc is priced out and there is some arc now priced
    in that has a reduced cost not far enough above that of the
    matching arc, price in the matching arc. We already know this
    condition on arcs we priced in, of course. Don't check them.
    */
    if (!v->matched->head->node_info.priced_in)
      if (one_priced_in)
	{
	a = v->matched;
	price_in_mch_arc(v, a);
#ifdef	QUICK_MIN
	need_best_rebuild = TRUE;
#endif
	}
    /*
    If a fix-in occurred on this node, unmatch it to preserve
    epsilon-optimality.
    */
    if (fix_this_node)
      {
#ifdef	DEBUG
      (void) printf("Fix-in -- unmatching (%ld, %ld)\n",
		    v - head_lhs_node + 1,
		    v->matched->head - head_rhs_node + 1);
#endif
      v->matched->head->matched = NULL;
      v->matched = NULL;
      total_e++;
      make_active(v);
      }
    }
  else
    {
    /*
    Node v is unmatched. Price any arc in whose reduced cost is less
    than po_cutoff above the minimum priced-in arc.
    */
    a = v->first;
    if (a != a_stop)
      {
#ifdef	EXPLICIT_LHS_PRICES
      v_price = v->p;
#else
      v_price = a->head->p - a->c;
      for (a++; a != a_stop; a++)
	if (v_price < (this_price = a->head->p - a->c))
	  v_price = this_price;
#endif
      for (a = v->priced_out; a != v->first; a++)
	if (v_price - (this_price = a->head->p - a->c) < po_cutoff)
	  {
	  price_in_unm_arc(v, a);
	  /*
	  If (this_price > v_price), we might have priced in some arcs
	  unnecessarily because the reduced costs of arcs incident to
	  v turn out to be higher than we thought. This is OK, but do
	  the right thing for the rest of the priced-out arcs.
	  */
	  if (this_price > v_price)
	    v_price = this_price;
#ifdef	QUICK_MIN
	  need_best_rebuild = TRUE;
#endif
	  /*
	  If we priced in the last priced-out arc in the list,
	  a == v->first, and we need to keep a from advancing too far.
	  */
	  if (a == v->first) break;
	  }
      }
    }
#ifdef	QUICK_MIN
  /*
  Make sure v->node_info.few_arcs reflects the priced-in degree of v.
  */
  if (a_stop - v->first < NUM_BEST + 1)
    v->node_info.few_arcs = TRUE;
  else
    {
    v->node_info.few_arcs = FALSE;
    if (need_best_rebuild)
      best_build(v);
    }
#endif
  }

#ifdef	DEBUG
(void) printf("Checked priced-out arcs. total_e=%lu\n", total_e);
#endif

if (fix_in) fix_ins++;
return(!fix_in);
}
