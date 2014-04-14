#include	<stdio.h>
#include	<math.h>
#include	"csa_types.h"
#include	"csa_defs.h"

extern	double	epsilon, min_epsilon, scale_factor;
extern	double	po_cost_thresh;
#ifdef	STRONG_PO
extern	double	banish_thresh;
#endif
extern	lhs_ptr	head_lhs_node, tail_lhs_node;
extern	rhs_ptr	head_rhs_node, tail_rhs_node;

void	incorporate_price_changes()

{
rhs_ptr	v;

for (v = head_rhs_node; v != tail_rhs_node; v++)
  {
  v->base_p += epsilon * v->p;
  v->p = 0;
  }
}

/*
Compute rounded reduced costs, and price out high-cost arcs. Operate
under the assumption that price changes have already been incorporated
into base prices.
*/

int	compute_rnd_red_costs()

{
lhs_ptr	v;
lr_aptr	a, a_start, a_stop;
double	v_price, this_price, red_cost, absolute_po_thresh, thresh,
	b_thresh;
#ifdef	STRONG_PO
double	absolute_banish_thresh, fix_in_thresh;
#endif
extern	void	exit();
int	fix_in = FALSE, one_priced_in;
#ifdef	BACK_PRICE_OUT
lhs_ptr	u;
rhs_ptr	w;
rl_aptr	b, b_stop;
#endif

absolute_po_thresh = po_cost_thresh * epsilon;
#ifdef	STRONG_PO
absolute_banish_thresh = banish_thresh * epsilon;
#endif
for (v = head_lhs_node; v != tail_lhs_node; v++)
  {
  a_stop = (v+1)->priced_out;
  if (a = v->matched)
    v_price = a->head->base_p - a->c_init;
  else
    {
    /*
    Figure out v's implicit price via the observation that its minimum
    reduced cost outgoing arc has reduced cost zero. Consider
    priced-in and priced-out arcs on an equal footing here. If there
    is no arc incident to this node, we don't need to do anything.
    */
    for (a = v->priced_out; a == v->matched; a++);
    if (a != a_stop)
      {
      v_price = a->head->base_p - a->c_init;
      for (a++; a != a_stop; a++)
	if ((a != v->matched) &&
	    (v_price < (this_price = a->head->base_p - a->c_init)))
	  v_price = this_price;
      }
    }
  /*
  For each arc incident to v, either
  - price it in and store its rounded reduced cost, or
  - price it out, or
  - store its rounded reduced cost.
  */
  a_start = v->first;
  thresh = absolute_po_thresh - v_price;
#ifdef	STRONG_PO
  b_thresh = absolute_banish_thresh - v_price;
#else
  b_thresh = thresh;
#endif
  one_priced_in = FALSE;
#ifdef	STRONG_PO
  /*
  Check for arcs to price in.
  */
  fix_in_thresh = -epsilon - v_price;
  for (a = v->priced_out; a != v->first; a++)
    {
    red_cost = a->c_init - a->head->base_p;
    if ((a != v->matched) && (red_cost < b_thresh))
      {
      a->c = (long) ((v_price + red_cost) / epsilon);
#ifdef	USE_P_UPDATE
      a->rev->c = a->c;
#endif
      if (red_cost < thresh)
	{
#ifdef	DEBUG
	(void) printf("upd_e: Pricing in (%ld, %ld), c = %ld\n",
		      v - head_lhs_node + 1,
		      a->head - head_rhs_node + 1, a->c);
#endif
	price_in_unm_arc(v, a);
	one_priced_in = TRUE;
	if (red_cost < fix_in_thresh)
	  fix_in = TRUE;
	if (a == v->first) break;
	}
      }
    }
#endif
  /*
  For each priced-in arc incident to v, either price it out or store
  its rounded reduced cost.
  */
  for (a = a_start; a != a_stop; a++)
    {
    if (a != v->matched)
      {
      red_cost = a->c_init - a->head->base_p;
      if (red_cost < b_thresh)
	{
	a->c = (long) ((v_price + red_cost) / epsilon);
#ifdef	USE_P_UPDATE
	a->rev->c = a->c;
#endif
#ifdef	CHECK_EPS_OPT
	if (a->c <= -3 * scale_factor)
	  {
	  (void) printf("Epsilon optimality violation : (%ld, %ld), c = %ld\n",
			v - head_lhs_node + 1, a->head - head_rhs_node + 1,
			a->c);
	  show_lhs_node(v - head_lhs_node + 1);
	  show_rhs_node(a->head - head_rhs_node + 1);
	  }
#endif
	}
#ifdef	DEBUG
      else
	{
	a->c = MAGIC_MARKER;
#ifdef	USE_P_UPDATE
	a->rev->c = -MAGIC_MARKER;
#endif
	}
#endif
      if (red_cost >= thresh)
	{
	/*
	Price this arc out.
	*/
#ifdef	DEBUG
	(void) printf("upd_e: pricing out (%ld, %ld), c = %ld\n",
		      v - head_lhs_node + 1,
		      a->head - head_rhs_node + 1, a->c);
#endif
	price_out_unm_arc(v, a);
	}
      else
	one_priced_in = TRUE;
      }
    }
  /*
  Now deal with the matching arc, if any.
  */
  if (a = v->matched)
    {
    a->c = 0;
#ifdef	USE_P_UPDATE
    a->rev->c = 0;
#endif
#ifdef	STRONG_PO
    if (one_priced_in)
      {
      if (!a->head->node_info.priced_in)
	{
	/*
	Matching arc is priced out and shouldn't be. Price it in.
	*/
#ifdef	DEBUG
	(void) printf("upd_e: Pricing in matched (%ld, %ld), c = %ld\n",
		      v - head_lhs_node + 1,
		      a->head - head_rhs_node + 1, a->c);
#endif
	price_in_mch_arc(v, a);
	}
      }
    else
#else
    if (!one_priced_in)
#endif
      if (a->head->node_info.priced_in)
	{
	/*
	No arcs are priced in except the matching arc. Price it out,
	too, and if we use back-arc price-outs, price out all the arcs
	incident to its head.
	*/
#ifdef	DEBUG
	(void) printf("upd_e: pricing out matched (%ld, %ld), c = %ld\n",
		      v - head_lhs_node + 1,
		      a->head - head_rhs_node + 1, a->c);
#endif
	price_out_mch_arc(v, a);
#ifdef	BACK_PRICE_OUT
	w = a->head;
	b_stop = (w+1)->priced_out;
	for (b = w->back_arcs; b != b_stop; b++)
	  {
	  u = b->tail;
	  a = b->rev;
#ifdef	DEBUG
	  (void) printf("upd_e: pricing out rev (%ld, %ld), c = %ld\n",
			u - head_lhs_node + 1,
			a->head - head_rhs_node + 1, a->c);
#endif
	  price_out_unm_arc(u, a);
	  }
#endif
	}
    }
  }
return(!fix_in);
}

int	update_epsilon()

{
incorporate_price_changes();

epsilon /= scale_factor;
if (epsilon < min_epsilon) epsilon = min_epsilon;

return(compute_rnd_red_costs());
}
