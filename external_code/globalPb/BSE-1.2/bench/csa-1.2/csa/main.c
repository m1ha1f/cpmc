/*
This program solves assignment problems posed in DIMACS format.
Written by Robert Kennedy
	   Department of Computer Science
	   Building 460
	   Stanford University
	   Stanford, CA 94305-2140
	   USA

in cooperation with Andrew V. Goldberg.

Revision history:
July, 1994	Option for (approximate) shortest-augmenting path
		finish-up stage (either "forward" or "backward" search
		added, along with a year's minor changes and speedups.
July, 1993	First distribution.

Copyright 1993, 1994 by Robert Kennedy
Permission to copy, exerpt, or prepare derivative works from this code
is granted provided this notice is included.

To obtain a copy of this program, send electronic mail to
ftp-request@theory.stanford.edu with the subject line "send csas.tar";
you will be mailed a uuencoded .tar file.

Please report any problems to robert@cs.stanford.edu.
*/

#include	<stdio.h>
#include	<math.h>
#include	<malloc.h>
#include	"csa_types.h"
#include	"csa_defs.h"

/* ------------------------- Problem size variables -------------------- */
unsigned	n, m;

/* --------------- Data structures describing the problem -------------- */
lhs_ptr	head_lhs_node, tail_lhs_node;
rhs_ptr	head_rhs_node, tail_rhs_node;
lr_aptr	head_lr_arc, tail_lr_arc;
#ifdef	STORE_REV_ARCS
rl_aptr	head_rl_arc, tail_rl_arc;
#endif

/* ------------------- Bookkeeping/profiling variables ----------------- */
unsigned	double_pushes = 0,
		pushes = 0,
		relabelings = 0,
		refines = 0,
		refine_time = 0;
#ifdef	USE_P_REFINE
unsigned	p_refines = 0,
		r_scans = 0,
		p_refine_time = 0;
#endif
#ifdef	USE_P_UPDATE
unsigned	p_updates = 0,
		u_scans = 0,
		p_update_time = 0;
#endif
#ifdef	USE_SP_AUG
unsigned	sp_augs = 0,
		a_scans = 0,
		sp_aug_time = 0;
#endif
#ifdef	STRONG_PO
unsigned	fix_ins = 0;
#endif
#ifdef	QUICK_MIN
unsigned	rebuilds = 0,
		scans = 0,
		non_scans = 0;
#endif

/* ------------------------- Tunable variables ------------------------- */
/*
Cost threshhold for pricing out: used even when price-outs are
switched off to make bounding reduced-cost differences easier in
double_push(). In principle this is not necessary, but it makes the
code better.
*/
double		po_cost_thresh;
double		scale_factor;	/* scaling factor */
#ifdef	USE_P_UPDATE
WORK_TYPE	upd_work_thresh;/* work threshhold for global update */
#endif
#ifdef	STRONG_PO
#ifdef	ROUND_COSTS
/*
Cost threshhold for certainty that an arc will never be priced in,
when strong price-outs are used. We need this because the possibility
of pricing an arc in requires us to maintain the rounded cost for that
arc in units of epsilon, but doing so for arcs with high reduced cost
generates integer overflows. Hence only those priced-out arcs that
aren't banished entirely have their rounded reduced costs calculated.
*/
double	banish_thresh;
#endif
WORK_TYPE	po_work_thresh;	/* work threshhold for price-in checks */
#endif

/*
Processing variables.
*/
double		epsilon;	/* scaling parameter */
double		min_epsilon;	/* snap to this value when epsilon small */
unsigned	total_e;	/* total excess */
ACTIVE_TYPE	active;		/* list of active nodes */
#ifdef	USE_P_REFINE
stack		reached_nodes;	/* nodes reached in topological ordering */
#endif
#if	defined(USE_P_REFINE) || defined(USE_P_UPDATE) || defined(USE_SP_AUG)
rhs_ptr		*bucket;	/* buckets for use in price refinements */
long		num_buckets;	/* number of buckets */
#endif

/*
Miscellaneous variables.
*/
char	*banner =
"==========================================================================";

void	show_usage(name)

char	*name;

{
void	exit();

(void) printf("Usage: %s [ scale [ update thresh [ price out thresh ] ] ]\n",
	      name);
exit(1);
}

void	parse_cmdline(argc, argv)

unsigned	argc;
char		*argv[];

{
char	*cmd = argv[0];

if (argc > 1)
  {
  if (sscanf(argv[1], "%lg", &scale_factor) == 0) show_usage(cmd);
  argc--; argv++;
  }
else
  scale_factor = DEFAULT_SCALE_FACTOR;

#ifdef	USE_P_UPDATE
if (argc > 1)
  {
  double	upd_fac;

  if (sscanf(argv[1], "%lg", &upd_fac) == 0) show_usage(cmd);
  argc--; argv++;
  upd_work_thresh = (unsigned) (upd_fac * (double) n);
  }
else
  upd_work_thresh = DEFAULT_UPD_FAC * n;
#endif

#ifdef	STRONG_PO
if (argc > 1)
  {
  if (sscanf(argv[1], "%lg", &po_cost_thresh) == 0) show_usage(cmd);
  argc--; argv++;
  }
else
  po_cost_thresh = DEFAULT_PO_COST_THRESH;
#ifdef	ROUND_COSTS
banish_thresh = 10.0 * (double) n * (scale_factor + 1);
#endif

if (argc > 1)
  {
  double	po_fac;

  if (sscanf(argv[1], "%lg", &po_fac) == 0) show_usage(cmd);
  argv--; argv++;
  po_work_thresh = (unsigned) (po_fac * (double) n);
  }
else
  po_work_thresh = DEFAULT_PO_WORK_THRESH * n;
#else
po_cost_thresh = 2.0 * (double) n * (scale_factor + 1);
#endif
}

void	describe_self()

{
static	char	*desc[20];
int	i = 0;
#ifdef	QUICK_MIN
char	minstr[40];
#endif

#ifdef	ROUND_COSTS
desc[i++] = "Rounded costs";
#endif
#ifdef	PREC_COSTS
desc[i++] = "Precise costs";
#ifdef	USE_PRICE_OUT
desc[i++] = "Price-outs";
#endif
#endif
#ifdef	STRONG_PO
desc[i++] = "Strong price-outs";
#endif
#ifdef	BACK_PRICE_OUT
desc[i++] = "Back price-outs";
#endif
#ifdef	USE_P_UPDATE
desc[i++] = "Global updates";
#endif
#ifdef	USE_SP_AUG_FORWARD
desc[i++] = "Forward SAP cleanup";
#endif
#ifdef	USE_SP_AUG_BACKWARD
desc[i++] = "Backward SAP cleanup";
#endif
#ifdef	STORE_REV_ARCS
desc[i++] = "Reverse arcs";
#endif
#ifdef	USE_P_REFINE
desc[i++] = "Price refinement";
#endif
#ifdef	QUEUE_ORDER
desc[i++] = "Queue ordering";
#else
desc[i++] = "Stack ordering";
#endif
#ifdef	QUICK_MIN
(void) sprintf(minstr, "Quick minima; NUM_BEST = %d", NUM_BEST);
desc[i++] = minstr;
#endif

desc[i] = NULL;

  (void) fprintf(stderr,"CSA: ");
for (i = 5; i > 0; i--)
  (void) fprintf(stderr,"=");

(void) fprintf(stderr," %s", desc[0]);
for (i = 1; desc[i]; i++)
  (void) fprintf(stderr,"; %s", desc[i]);
(void) fprintf(stderr," ");
for (i = 5; i > 0; i--)
  (void) fprintf(stderr,"=");
(void) fprintf(stderr,"\n");
}

void	init(argc, argv)

unsigned	argc;
char		*argv[];

{
void	exit();
extern	unsigned long	parse();
#ifdef	QUEUE_ORDER
extern	queue		q_create();
#endif
#if	defined(USE_P_REFINE) || !defined(QUEUE_ORDER)
/*
was #if	defined(USE_P_REFINE) || defined(USE_P_UPDATE) || !defined(QUEUE_ORDER)
*/
extern	stack		st_create();
#endif
#if	defined(USE_P_REFINE) || defined(USE_P_UPDATE) || defined(USE_SP_AUG)
rhs_ptr	r_v;
long	i;
#endif
#ifdef	QUICK_MIN
lhs_ptr	l_v;
void	best_build();
#endif

describe_self();
 fprintf(stderr,"CSA: parsing input...\n");
epsilon = parse();
 fprintf(stderr,"CSA: parsing cmd line...\n");
parse_cmdline(argc, argv);
 fprintf(stderr,"CSA: buiding...\n");

create_active(n);
#ifdef	USE_P_REFINE
reached_nodes = st_create(n);
#endif
#if	defined(USE_P_REFINE) || defined(USE_P_UPDATE) || defined(USE_SP_AUG)
#ifdef	PREC_COSTS
num_buckets = scale_factor * n + 1;
#else
num_buckets = 2 * scale_factor * n + 1;
#endif
bucket = (rhs_ptr *) malloc((unsigned) num_buckets * sizeof(rhs_ptr));
if (bucket == NULL)
  {
  (void) printf("Insufficient memory.\n");
  exit(9);
  }
for (i = 0; i < num_buckets; i++)
  bucket[i] = tail_rhs_node;
for (r_v = head_rhs_node; r_v != tail_rhs_node; r_v++)
  r_v->key = num_buckets;
#endif
#ifdef	QUICK_MIN
for (l_v = head_lhs_node; l_v != tail_lhs_node; l_v++)
  if (!l_v->node_info.few_arcs)
    best_build(l_v);
/*
Count only those builds that take place after initialization; first
setup is free.
*/
rebuilds = 0;
#endif
}

double	compute_cost()

{
double	cost = 0.0;
lhs_ptr	v;

for (v = head_lhs_node; v != tail_lhs_node; v++)
  if (v->matched)
#ifdef	ROUND_COSTS
    cost += v->matched->c_init;
#else
    cost += v->matched->c;
#endif

return(cost);
}

void	display_results(time)

unsigned	time;

{
#ifdef	SAVE_RESULT
lhs_ptr	v;
FILE	*f;
double	edge_cost;
#endif

(void) fprintf(stderr,"CSA: |>   cost %17.0f,    time %10.3f seconds\n",
	      compute_cost(), (double) time / 60.0);
/*
Avoid division by zero.
*/
if (time == 0) time = 1;
(void) fprintf(stderr,"CSA: |>   %u refines:     %lg%%     %u relabelings\n",
	      refines, 100.0 * (double) refine_time / (double) time,
	      relabelings);
(void) fprintf(stderr,"CSA: |>                   %u double pushes, %u pushes\n",
	      double_pushes, pushes);
#ifdef	USE_P_REFINE
(void) fprintf(stderr,"CSA: |>   %u p_refines: %lg%%      %u r_scans\n",
	      p_refines, 100.0 * (double) p_refine_time / (double) time,
	      r_scans);
#endif
#ifdef	USE_P_UPDATE
(void) fprintf(stderr,"CSA: |>   %u p_updates: %lg%%      %u u_scans\n",
	      p_updates, 100.0 * (double) p_update_time / (double) time,
	      u_scans);
#endif
#ifdef	USE_SP_AUG
(void) fprintf(stderr,"CSA: |>   %u sp_augs:   %lg%%      %u a_scans\n",
	      sp_augs, 100.0 * (double) sp_aug_time / (double) time,
	      a_scans);
#endif
#ifdef	STRONG_PO
(void) fprintf(stderr,"CSA: |>   %u fix-ins\n", fix_ins);
#endif
#ifdef	QUICK_MIN
(void) fprintf(stderr,"CSA: |>   %u list rebuilds, %u full scans, %u avoided scans\n",
	      rebuilds, scans, non_scans);
#endif
(void) puts(banner);
#ifdef	SAVE_RESULT
/*f = fopen("output.flow", "w");*/
f = fdopen(1,"w");
for (v = head_lhs_node; v != tail_lhs_node; v++)
  {
#ifdef	ROUND_COSTS
  edge_cost = -v->matched->c_init;
#else
  edge_cost = -v->matched->c;
#endif
  (void) fprintf(f, "f %lu %lu %.0lf\n",
		 v - head_lhs_node + 1,
		 v->matched->head - head_rhs_node + 1 +
		 tail_lhs_node - head_lhs_node,
		 edge_cost);
  }
(void) fclose(f);
#endif
}

int	main(argc, argv)

unsigned	argc;
char		*argv[];

{
unsigned	time, myclock();
extern	int	update_epsilon();
extern	void	refine();
extern	int	p_refine();

 fprintf(stderr,"CSA: initializing...\n");
init(argc, argv);
 fprintf(stderr,"CSA: solving...\n");

(void) fprintf(stderr,"CSA: |>  n = %u,  m = %u,  sc_f = %lg", n, m, scale_factor);
#if	defined(USE_PRICE_OUT) || defined(ROUND_COSTS)
(void) fprintf(stderr,",  po_thr = %lg", po_cost_thresh);
#endif
(void) fprintf(stderr,"\n");

#ifdef	PREC_COSTS
min_epsilon = 2.0 / (double) (n + 1);
#else
min_epsilon = 1.0 / (double) (n + 1);
#endif

time = myclock();

#ifdef	USE_P_REFINE
(void) update_epsilon();
refine();
#endif

while (epsilon > min_epsilon)
  {
#ifdef	VERBOSE_TIME
  (void) fprintf(stderr,"CSA: |>   Epsilon = %lg; time = %lg\n",
		epsilon, ((double) (myclock() - time)) / 60.0);
#endif
#ifdef	USE_P_REFINE
  if (!update_epsilon() || !p_refine())
    refine();
#else
  (void) update_epsilon();
  refine();
#endif
  }

time = myclock() - time;

display_results(time);
return(0);
}
