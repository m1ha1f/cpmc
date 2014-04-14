#include	<stdio.h>
#include	"csa_types.h"
#include	"csa_defs.h"

extern	unsigned	total_e;
extern	ACTIVE_TYPE	active;
extern	unsigned	fix_ins;
extern	lhs_ptr		head_lhs_node, tail_lhs_node;
extern	double		po_cost_thresh, epsilon;
#ifdef	DEBUG
extern	rhs_ptr		head_rhs_node;
#endif

int	check_po_arcs()

{
lhs_ptr	v;
rhs_ptr	w;
lr_aptr	a, a_start, a_stop;
int	one_priced_in, fix_this_node, fix_in = FALSE;
double	match_rc, this_cost, v_price, this_price, po_cutoff, thresh,
	fix_in_thresh;

po_cutoff = po_cost_thresh * epsilon;
for (v = head_lhs_node; v != tail_lhs_node; v++)
  {
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
    w = a->head;
#ifdef	DEBUG
    (void) printf("c_p_a: matching arc (%ld, %ld), c = %ld\n",
		  v - head_lhs_node + 1, w - head_rhs_node + 1,
		  a->c - w->p);
#endif
    match_rc = a->c_init - w->base_p - epsilon * w->p;
    thresh = match_rc + po_cutoff;
    fix_in_thresh = match_rc - epsilon;
    for (a = v->priced_out; a != v->first; a++)
      {
      w = a->head;
      if ((a != v->matched) &&
	  ((this_cost = a->c_init - w->base_p - epsilon * w->p) < thresh))
	{
#ifdef	DEBUG
	(void) printf("c_p_a: matched, pricing in (%ld, %ld), c = %ld\n",
		      v - head_lhs_node + 1,
		      w - head_rhs_node + 1, a->c - w->p);
#endif
	price_in_unm_arc(v, a);
	one_priced_in = TRUE;
	if (this_cost < fix_in_thresh)
	  {
	  /*
	  Epsilon-optimality violated by priced-in arc.
	  */
	  fix_in = TRUE;
	  fix_this_node = TRUE;
	  }
	/*
	If we priced in the last priced-out arc in the list,
	a == v->first, and we need to keep a from advancing too far.
	*/
	if (a == v->first) break;
	}
      }
    /*
    Now if matching arc is priced out and there is some arc now priced
    in, price in the matching arc.
    */
    if (!(a = v->matched)->head->node_info.priced_in)
      {
      if (one_priced_in)
	{
#ifdef	DEBUG
	(void) printf("c_p_a: pricing in matching arc: (%ld, %ld), c = %ld\n",
		      v - head_lhs_node + 1, a->head - head_rhs_node + 1,
		      a->c - a->head->p);
#endif
	price_in_mch_arc(v, a);
	a = v->matched;
	}
      }
    /*
    If a fix-in occurred on this node, unmatch it to preserve
    epsilon-optimality.
    */
    if (fix_this_node)
      {
      a->head->matched = NULL;
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
      w = a->head;
      v_price = w->base_p + epsilon * w->p - a->c_init;
      for (a++; a != a_stop; a++)
	{
	w = a->head;
	if (v_price < (this_price = w->base_p + epsilon * w->p - a->c_init))
	  v_price = this_price;
	}
      for (a = v->priced_out; a != v->first; a++)
	{
	w = a->head;
	if (v_price -
	    (this_price = epsilon * w->p + w->base_p - a->c_init) < po_cutoff)
	  {
#ifdef	DEBUG
	  (void) printf("c_p_a: unmatched, pricing in (%ld, %ld), c = %ld\n",
			v - head_lhs_node + 1, a->head - head_rhs_node + 1,
			a->c - w->p);
#endif
	  price_in_unm_arc(v, a);
	  /*
	  If (this_price > v_price), we might be pricing in some arcs
	  unnecessarily because the reduced costs of arcs incident to
	  v turn out to be higher than we thought. This is OK, but do
	  the right thing for the rest of the priced-out arcs.
	  */
	  if (this_price > v_price)
	    v_price = this_price;
	  /*
	  If we priced in the last priced-out arc in the list,
	  a == v->first, and we need to keep a from advancing too far.
	  */
	  if (a == v->first) break;
	  }
	}
      }
    }
  }
if (fix_in) fix_ins++;
return(!fix_in);
}
