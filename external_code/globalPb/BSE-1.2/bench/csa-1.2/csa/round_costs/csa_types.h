#define	ROUND_COSTS

#if	defined(USE_P_UPDATE) || defined(BACK_PRICE_OUT)
#define	STORE_REV_ARCS
#endif

typedef	struct	lhs_node	{
				/*
				first arc in the arc array associated
				with this node.
				*/
				struct	lr_arc	*priced_out;
				/*
				first priced-in arc in the arc array
				associated with this node.
				*/
				struct	lr_arc	*first;
				/*
				matching arc (if any) associated with
				this node; NULL if this node is
				unmatched.
				*/
				struct	lr_arc	*matched;
#ifdef	USE_P_UPDATE
				/*
				price change required on this node (in
				units of epsilon) to ensure that its
				excess can reach a deficit in the
				admissible graph. computed and used in
				p_update().
				*/
				long	delta_reqd;
#endif
				}	*lhs_ptr;

typedef	struct	rhs_node	{
				struct	{
#ifdef	USE_P_REFINE
					/*
					depth-first search flags.
					dfs is to determine whether
					admissible graph contains a
					cycle in p_refine().
					*/
					unsigned	srchng : 1;
					unsigned	srched : 1;
#endif
					/*
					flag to indicate this node's
					matching arc (if any) is
					priced in.
					*/
					unsigned	priced_in : 1;
					}	node_info;
				/*
				lhs node this rhs node is matched to.
				*/
				lhs_ptr	matched;
				/*
				price of this node at end of previous
				iteration.
				*/
				double	base_p;
				/*
				price change (in units of epsilon) to
				this node during the current
				iteration.
				*/
				long	p;
#if	defined(USE_P_REFINE) || defined(USE_P_UPDATE)
				/*
				number of epsilons of price change
				required at this node to accomplish
				p_refine()'s or p_update()'s goal.
				*/
				long	key;
				/*
				fields to maintain buckets of nodes as
				lists in p_refine() and p_update().
				*/
				struct	rhs_node	*prev, *next;
#endif
#ifdef	STORE_REV_ARCS
				/*
				first back arc in the arc array
				associated with this node.
				*/
				struct	rl_arc	*priced_out;
				/*
				first priced-in back arc in the arc
				array associated with this node.
				*/
				struct	rl_arc	*back_arcs;
#endif
				}	*rhs_ptr;

#ifdef	STORE_REV_ARCS
typedef	struct	rl_arc		{
				/*
				lhs node associated with this back
				arc. some would have liked the name
				head better.
				*/
				lhs_ptr	tail;
#ifdef	USE_P_UPDATE
				/*
				cost of this back arc (in units of
				epsilon). this cost gets modified to
				incorporate other arc costs in
				p_update(), while forward arc costs
				remain constant throughout an
				iteration.
				*/
				long	c;
#endif
				/*
				this arc's reverse in the forward arc
				list. we have to store this even if we
				don't do back-arc price-outs because
				when we price a reverse arc out, we
				need to be able to find the forward
				arc corresponding to the back arc we
				swap it with.
				*/
				struct	lr_arc	*rev;
				}	*rl_aptr;
#endif

typedef	struct	lr_arc		{
				/*
				rhs node associated with this arc.
				*/
				rhs_ptr	head;
				/*
				initial arc cost.
				*/
				double	c_init;
				/*
				rounded reduced arc cost at beginning
				of current iteration.
				*/
				long	c;
#ifdef	STORE_REV_ARCS
				/*
				this arc's reverse in the back arc
				list.
				*/
				struct	rl_arc	*rev;
#endif
				}	*lr_aptr;

typedef	struct	stack_st	{
				/*
				Sometimes stacks have lhs nodes, and
				other times they have rhs nodes. So
				there's a little type clash;
				everything gets cast to (char *) so we
				can use the same structure for both.
				*/
				char	**bottom;
				char	**top;
				}	*stack;

typedef	struct	queue_st	{
				/*
				Sometimes queues have lhs nodes, and
				other times they have rhs nodes. So
				there's a little type clash;
				everything gets cast to (char *) so we
				can use the same structure for both.
				*/
				char		**head;
				char		**tail;
				char		**storage;
				char		**end;
				unsigned	max_size;
				}	*queue;
