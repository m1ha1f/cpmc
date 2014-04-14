#include	<stdio.h>
#include	"csa_types.h"

extern	lhs_ptr	head_lhs_node, tail_lhs_node;
extern	rhs_ptr	head_rhs_node, tail_rhs_node;
extern	lr_aptr	head_arc;

void	show_lhs_node(lhs_id)

int	lhs_id;

{
lhs_ptr	v = &head_lhs_node[lhs_id - 1];
int	rhs_id;
rhs_ptr	w;
lr_aptr	a;

(void) printf("Lhs node %d ", lhs_id);
if (v->matched)
  {
  w = v->matched->head;
  rhs_id = w - head_rhs_node + 1;
  (void) printf("matched thru cost %lg (stored cost %ld) to rhs node %d",
		v->matched->c_init, v->matched->c, rhs_id);
  if (w->matched == v)
    (void) putchar('\n');
  else
    {
    lhs_id = w->matched - head_lhs_node + 1;
    (void) printf(", matched back to lhs node %d\n", lhs_id);
    }
  (void) printf("\tMatching arc is priced ");
  if ((long) v->matched - (long) v->first >= 0)
    (void) printf("in\n");
  else
    (void) printf("out\n");
  }
else
  (void) printf("unmatched\n");
(void) printf("\t%d arcs priced out, %d arcs priced in\n",
	      v->first - v->priced_out, (v+1)->priced_out - v->first);
if ((v+1)->priced_out - v->first > 0)
  {
  (void) printf("\tPriced in arcs:\n");
  for (a = v->first; a != (v+1)->priced_out; a++)
    {
    rhs_id = a->head - head_rhs_node + 1;
    (void) printf("\t\t(%d, %d) cost %lg, stored cost %ld, cmp cost %ld\n",
		  lhs_id, rhs_id, a->c_init, a->c, a->c - a->head->p);
    }
  }
}

void	show_rhs_node(rhs_id)

int	rhs_id;

{
rhs_ptr	v = &head_rhs_node[rhs_id - 1];
int	lhs_id;

(void) printf("Rhs node %d, base_p %lg, delta p %ld ",
	      rhs_id, v->base_p, v->p);
if (v->matched)
  {
  lhs_id = v->matched - head_lhs_node + 1;
  if (v->matched->matched->head == v)
    (void) printf("matched thru cost %lg to lhs node %d\n",
		  v->matched->matched->c_init, lhs_id);
  else
    (void) printf("matched inconsistently to lhs node %d\n", lhs_id);
  }
else
  (void) printf("unmatched\n");
}

void	show_lhs()

{
int	id;

for (id = 1; id <= tail_lhs_node - head_lhs_node; id++)
  show_lhs_node(id);
}

void	show_rhs()

{
int	id;

for (id = 1; id <= tail_rhs_node - head_rhs_node; id++)
  show_rhs_node(id);
}
