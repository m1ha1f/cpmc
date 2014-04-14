#include	<stdio.h>
#include	<malloc.h>
#include	"csa_types.h"
#include	"csa_defs.h"

#define ERRBASE		1000	/* Base number for user-defined errors*/
#define BADINPUT1	1001	/* Bad input file format */
#define BADINPUT2	1002	/* Bad input file format */
#define BADINPUT3	1003	/* Bad input file format */
#define BADINPUT4	1004	/* Bad input file format */
#define BADINPUT5	1005	/* Bad input file format */
#define BADCOUNT	1006	/* Arc count discrepancy */
#define	NONCONTIG	1007	/* Node id numbers not contiguous */
#define NOMEM		1008	/* Not enough memory */
  
char *err_messages[] =
{
"Can't read from the input file.",
"Not a correct assignment problem line.",
"Error reading a node descriptor from the input.",
"Error reading an arc descriptor from the input.",
"Unknown line type in the input",
"Inconsistent number of arcs in the input.",
"Parsing noncontiguous node ID numbers not implemented.",
"Can't obtain enough memory to solve this problem.",
};

void parse_error(err_index)

int err_index;

{
void	exit();

(void) fprintf(stderr,"CSA: Error while parsing the input: %s \n",
	      err_messages[(err_index % ERRBASE) - 1]);
exit(1);
}

typedef	struct	temp_arc	{
				lhs_ptr	tail;
				rhs_ptr	head;
				long	cost;
				}	*ta_ptr;

extern	char	*banner;
extern	unsigned	m, n;
extern	lhs_ptr	head_lhs_node, tail_lhs_node;
extern	rhs_ptr	head_rhs_node, tail_rhs_node;
extern	lr_aptr	head_lr_arc, tail_lr_arc;
#ifdef	STORE_REV_ARCS
extern	rl_aptr	head_rl_arc, tail_rl_arc;
#endif

unsigned long	parse()

{
char	line_type, prob_type[4], in_line[MAXLINE];
unsigned	lhs_known = FALSE, arc_count, tail, head, lhs_n, node_id,
		swap, id_offset, temp;
long	cost, *lhs_degree;
#ifdef	STORE_REV_ARCS
long	*rhs_degree;
rl_aptr	b;
#endif
unsigned long	max_cost = 0;
lr_aptr	a;
ta_ptr	temp_a, temp_arcs;
lhs_ptr	l_v;
rhs_ptr	r_v;
extern	int	abs();

/* skip initial comments */
do
  if (fgets(in_line,MAXLINE,stdin) == NULL)
    parse_error(BADINPUT1);
while (in_line[0] == 'c');

if ((sscanf(in_line, "%c%3s%d%d", &line_type, prob_type, &n, &m) != 4) ||
    (line_type != 'p') ||
    (prob_type[0] != 'a') ||
    (prob_type[1] != 's') ||
    (prob_type[2] != 'n'))
  parse_error(BADINPUT2);

arc_count = 0;
lhs_n = 0;

while (fgets(in_line,MAXLINE,stdin) != NULL) {
  switch (in_line[0])
    {
    case 'c': break;

    case 'n':
      if (sscanf(in_line, "%*c%d", &node_id) == 0)
	parse_error(BADINPUT3);
/*fprintf(stderr,"n %d\n", node_id);*/
      if (node_id != ++lhs_n)
	parse_error(NONCONTIG);
      break;

    case 'a':
      if (!lhs_known)
	{
	lhs_known = TRUE;
	head_lr_arc = (lr_aptr) malloc((m + 1) * sizeof(struct lr_arc));
	tail_lr_arc = head_lr_arc + m;
#ifdef	STORE_REV_ARCS
	head_rl_arc = (rl_aptr) malloc((m + 1) * sizeof(struct rl_arc));
	tail_rl_arc = head_rl_arc + m;
#endif
	id_offset = lhs_n;
	if (lhs_n > n - lhs_n)
	  {
	  lhs_n = n - lhs_n;
	  swap = TRUE;
	  }
	else
	  swap = FALSE;
	head_lhs_node = (lhs_ptr) malloc((lhs_n + 1) *
					 sizeof(struct lhs_node));
	tail_lhs_node = head_lhs_node + lhs_n;
	head_rhs_node = (rhs_ptr) malloc((n - lhs_n + 1) *
					 sizeof(struct rhs_node));
	tail_rhs_node = head_rhs_node + n - lhs_n;
	lhs_degree = (long *) malloc(lhs_n * sizeof(long));
#ifdef	STORE_REV_ARCS
	rhs_degree = (long *) malloc((n - lhs_n) * sizeof(long));
	if ((rhs_degree == NULL) || (head_rl_arc == NULL))
	  parse_error(NOMEM);
	for (tail = 0; tail < n - lhs_n; tail++)
	  rhs_degree[tail] = 0;
#endif
	temp_arcs = (ta_ptr) malloc(m * sizeof(struct temp_arc));
	if ((head_lhs_node == NULL) || (head_lr_arc == NULL) ||
	    (lhs_degree == NULL) || (temp_arcs == NULL))
	  parse_error(NOMEM);
	temp_a = temp_arcs;
	for (tail = 0; tail < lhs_n; tail++)
	  lhs_degree[tail] = 0;
	(void) puts(banner);
	}
      if (sscanf(in_line, "%*c%d%d%ld", &tail, &head, &cost) == 0)
	parse_error(BADINPUT4);
/*fprintf(stderr,"a %d %d %ld\n", tail,head,cost);*/

      head -= id_offset;
      if (swap)
	{
	temp = head;
	head = tail;
	tail = temp;
	}

      if ((tail < 1) || (tail > lhs_n) ||
	  (head < 1) || (head > n - lhs_n))
	parse_error(BADINPUT4);

      arc_count++;

      if (arc_count > m)
	parse_error(BADCOUNT);
      head--; tail--;
      temp_a->head = head_rhs_node + head;
      temp_a->tail = head_lhs_node + tail;
      temp_a->cost = cost;
      if ((cost = abs((int) cost)) > max_cost) max_cost = cost;
      temp_a++;
      lhs_degree[tail]++;
#ifdef	STORE_REV_ARCS
      rhs_degree[head]++;
#endif
      break;

    case '\n': break;
    case 0: break;

    default:
      parse_error(BADINPUT5);
      break;
    }
}
fprintf(stderr,"CSA: Done parsing input.\n");

if (arc_count != m)
  parse_error(BADCOUNT);

a = head_lr_arc;
for (tail = 0, l_v = head_lhs_node; l_v != tail_lhs_node; l_v++, tail++)
  {
  l_v->priced_out = l_v->first = a;
  l_v->matched = NULL;
  a += lhs_degree[tail];
#ifdef	QUICK_MIN
  if (lhs_degree[tail] < NUM_BEST + 1)
    l_v->node_info.few_arcs = TRUE;
  else
    l_v->node_info.few_arcs = FALSE;
#endif
  }
tail_lhs_node->priced_out = a;

#ifdef	STORE_REV_ARCS
tail = 0;
b = head_rl_arc;
#endif
for (r_v = head_rhs_node; r_v != tail_rhs_node; r_v++)
  {
  r_v->node_info.priced_in = TRUE;
  r_v->matched = NULL;
#ifdef	STORE_REV_ARCS
  r_v->priced_out = r_v->back_arcs = b;
  b += rhs_degree[tail];
  tail++;
#endif
#ifdef	ROUND_COSTS
  r_v->base_p = 0.0;
  r_v->p = 0;
#else
  r_v->p = 0.0;
#endif
  }
#ifdef	STORE_REV_ARCS
tail_rhs_node->priced_out = b;
#endif

for (temp_a--; temp_a != temp_arcs - 1; temp_a--)
  {
  a = temp_a->tail->first + (--lhs_degree[temp_a->tail - head_lhs_node]);
  a->head = temp_a->head;
#ifdef	ROUND_COSTS
#ifdef	MIN_COST
  a->c_init = (double) temp_a->cost;
#else
  a->c_init = (double) -temp_a->cost;
#endif
#else	/* PREC_COSTS */
#ifdef	MIN_COST
  a->c = (double) temp_a->cost;
#else
  a->c = (double) -temp_a->cost;
#endif
#endif	/* ROUND_COSTS */
#ifdef	USE_SP_AUG_FORWARD
  a->tail = temp_a->tail;
#endif
#ifdef	STORE_REV_ARCS
  b = temp_a->head->back_arcs + (--rhs_degree[temp_a->head - head_rhs_node]);
  a->rev = b;
#if	defined(ROUND_COSTS) || defined(USE_PRICE_OUT) || \
	defined(USE_SP_AUG_BACKWARD)
  b->rev = a;
#endif
  b->tail = temp_a->tail;
#if	defined(USE_P_UPDATE) || defined(USE_SP_AUG_BACKWARD)
/*
In the ROUND_COSTS case, update_epsilon() takes care of b->c.
*/
#ifdef	PREC_COSTS
  b->c = a->c;
#endif
#endif
#endif
  }

(void) free((char *) temp_arcs);
(void) free((char *) lhs_degree);
#ifdef	STORE_REV_ARCS
(void) free((char *) rhs_degree);
#endif

return(max_cost);
}
