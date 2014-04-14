#
# Makefile for rounded-costs (bit-scaling) codes implementing the cost
# scaling push-relabel assignment algorithm.
# Filenames of executable images are derived from compile-time macro
# definitions as follows (for more information on the compile-time
# macros used here, see the README file).
# Filename fragment "_s" denotes stack ordering.
# Filename fragment "_q" denotes queue ordering.
# Filename fragment "_pr" denotes price refinement.
# Filename fragment "_pu" denotes global price updates.
# Filename fragment "_tpo" denotes theoretically-justified price-outs.
# Filename fragment "_rtpo" denotes theoretically-justified price-outs
#   with additional arcs priced out based on incidence with a node
#   whose matching arc is priced out.
# Filename fragment "_spo" denotes strong price-outs (speculative
#   arc-fixing).
# Filename fragment "_rspo" denotes strong price-outs (speculative
#   arc-fixing) with additional arcs priced out based on incidence
#   with a node whose matching arc is priced out.
#
BASEFILES=main.c refine.c update_epsilon.c parse.c stack.c timer.c debug.c
HEADERS=csa_types.h csa_defs.h
CFLAGS=-O3

links:		#main.c parse.c stack.c timer.c list.c queue.c csa_defs.h
	ln -s ../main.c main.c
	ln -s ../parse.c parse.c
	ln -s ../stack.c stack.c
	ln -s ../timer.c timer.c
	ln -s ../list.c list.c
	ln -s ../queue.c queue.c
	ln -s ../csa_defs.h csa_defs.h
	touch links

csa_s_tpo:	$(BASEFILES) $(HEADERS)
	cc $(CFLAGS) -o $@ $(BASEFILES)

csa_q_tpo:	$(BASEFILES) $(HEADERS) queue.c
	cc $(CFLAGS) -DQUEUE_ORDER -o $@ $(BASEFILES) queue.c

csa_s_rtpo:	$(BASEFILES) $(HEADERS)
	cc $(CFLAGS) -DBACK_PRICE_OUT -o $@ $(BASEFILES)

csa_q_rtpo:	$(BASEFILES) $(HEADERS) queue.c
	cc $(CFLAGS) -DBACK_PRICE_OUT -DQUEUE_ORDER -o $@ $(BASEFILES) queue.c

csa_s_spo:	$(BASEFILES) $(HEADERS) check_po_arcs.c
	cc $(CFLAGS) -DSTRONG_PO -o $@ $(BASEFILES) check_po_arcs.c -lm

csa_q_spo:	$(BASEFILES) $(HEADERS) queue.c check_po_arcs.c
	cc $(CFLAGS) -DSTRONG_PO -DQUEUE_ORDER -o $@ $(BASEFILES) queue.c check_po_arcs.c -lm

csa_s_rspo:	$(BASEFILES) $(HEADERS) check_po_arcs.c
	cc $(CFLAGS) -DSTRONG_PO -DBACK_PRICE_OUT -o $@ $(BASEFILES) check_po_arcs.c -lm

csa_q_rspo:	$(BASEFILES) $(HEADERS) queue.c check_po_arcs.c
	cc $(CFLAGS) -DSTRONG_PO -DBACK_PRICE_OUT -DQUEUE_ORDER -o $@ $(BASEFILES) queue.c check_po_arcs.c -lm

csa_s_tpo_pr:	$(BASEFILES) $(HEADERS) p_refine.c list.c
	cc $(CFLAGS) -DUSE_P_REFINE -o $@ $(BASEFILES) p_refine.c list.c -lm

csa_q_tpo_pr:	$(BASEFILES) $(HEADERS) p_refine.c queue.c list.c
	cc $(CFLAGS) -DUSE_P_REFINE -DQUEUE_ORDER -o $@ $(BASEFILES) p_refine.c queue.c list.c -lm

csa_s_rtpo_pr:	$(BASEFILES) $(HEADERS) p_refine.c list.c
	cc $(CFLAGS) -DBACK_PRICE_OUT -DUSE_P_REFINE -o $@ $(BASEFILES) p_refine.c list.c -lm

csa_q_rtpo_pr:	$(BASEFILES) $(HEADERS) p_refine.c queue.c list.c
	cc $(CFLAGS) -DBACK_PRICE_OUT -DUSE_P_REFINE -DQUEUE_ORDER -o $@ $(BASEFILES) p_refine.c queue.c list.c -lm

csa_s_spo_pr:	$(BASEFILES) $(HEADERS) p_refine.c check_po_arcs.c list.c
	cc $(CFLAGS) -DSTRONG_PO -DUSE_P_REFINE -o $@ $(BASEFILES) p_refine.c check_po_arcs.c list.c -lm

csa_q_spo_pr:	$(BASEFILES) $(HEADERS) p_refine.c check_po_arcs.c queue.c list.c
	cc $(CFLAGS) -DSTRONG_PO -DUSE_P_REFINE -DQUEUE_ORDER -o $@ $(BASEFILES) p_refine.c check_po_arcs.c queue.c list.c -lm

csa_s_rspo_pr:	$(BASEFILES) $(HEADERS) p_refine.c check_po_arcs.c list.c
	cc $(CFLAGS) -DSTRONG_PO -DUSE_P_REFINE -DBACK_PRICE_OUT -o $@ $(BASEFILES) p_refine.c check_po_arcs.c list.c -lm

csa_q_rspo_pr:	$(BASEFILES) $(HEADERS) p_refine.c check_po_arcs.c queue.c list.c
	cc $(CFLAGS) -DSTRONG_PO -DUSE_P_REFINE -DBACK_PRICE_OUT -DQUEUE_ORDER -o $@ $(BASEFILES) p_refine.c check_po_arcs.c queue.c list.c -lm

csa_s_tpo_pu:	$(BASEFILES) $(HEADERS) p_update.c list.c
	cc $(CFLAGS) -DUSE_P_UPDATE -o $@ $(BASEFILES) p_update.c list.c -lm

csa_q_tpo_pu:	$(BASEFILES) $(HEADERS) p_update.c queue.c list.c
	cc $(CFLAGS) -DUSE_P_UPDATE -DQUEUE_ORDER -o $@ $(BASEFILES) p_update.c queue.c list.c -lm

csa_s_rtpo_pu:	$(BASEFILES) $(HEADERS) p_update.c list.c
	cc $(CFLAGS) -DBACK_PRICE_OUT -DUSE_P_UPDATE -o $@ $(BASEFILES) p_update.c list.c -lm

csa_q_rtpo_pu:	$(BASEFILES) $(HEADERS) p_update.c queue.c list.c
	cc $(CFLAGS) -DBACK_PRICE_OUT -DUSE_P_UPDATE -DQUEUE_ORDER -o $@ $(BASEFILES) p_update.c queue.c list.c -lm

csa_s_spo_pu:	$(BASEFILES) $(HEADERS) p_update.c check_po_arcs.c list.c
	cc $(CFLAGS) -DSTRONG_PO -DUSE_P_UPDATE -o $@ $(BASEFILES) p_update.c check_po_arcs.c list.c -lm

csa_q_spo_pu:	$(BASEFILES) $(HEADERS) p_update.c check_po_arcs.c queue.c list.c
	cc $(CFLAGS) -DSTRONG_PO -DUSE_P_UPDATE -DQUEUE_ORDER -o $@ $(BASEFILES) p_update.c check_po_arcs.c queue.c list.c -lm

csa_s_tpo_pr_pu:	$(BASEFILES) $(HEADERS) p_refine.c p_update.c list.c
	cc $(CFLAGS) -DUSE_P_UPDATE -DUSE_P_REFINE -o $@ $(BASEFILES) p_refine.c p_update.c list.c -lm

csa_q_tpo_pr_pu:	$(BASEFILES) $(HEADERS) p_refine.c p_update.c queue.c list.c
	cc $(CFLAGS) -DUSE_P_UPDATE -DUSE_P_REFINE -DQUEUE_ORDER -o $@ $(BASEFILES) p_refine.c p_update.c queue.c list.c -lm

csa_s_rtpo_pr_pu:	$(BASEFILES) $(HEADERS) p_refine.c p_update.c list.c
	cc $(CFLAGS) -DUSE_P_UPDATE -DUSE_P_REFINE -DBACK_PRICE_OUT -o $@ $(BASEFILES) p_refine.c p_update.c list.c -lm

csa_q_rtpo_pr_pu:	$(BASEFILES) $(HEADERS) p_refine.c p_update.c queue.c list.c
	cc $(CFLAGS) -DUSE_P_UPDATE -DUSE_P_REFINE -DBACK_PRICE_OUT -DQUEUE_ORDER -o $@ $(BASEFILES) p_refine.c p_update.c queue.c list.c -lm

csa_s_spo_pr_pu:	$(BASEFILES) $(HEADERS) p_refine.c p_update.c check_po_arcs.c list.c
	cc $(CFLAGS) -DSTRONG_PO -DUSE_P_UPDATE -DUSE_P_REFINE -o $@ $(BASEFILES) p_refine.c p_update.c check_po_arcs.c list.c -lm

csa_q_spo_pr_pu:	$(BASEFILES) $(HEADERS) p_refine.c p_update.c check_po_arcs.c queue.c list.c
	cc $(CFLAGS) -DSTRONG_PO -DUSE_P_UPDATE -DUSE_P_REFINE -DQUEUE_ORDER -o $@ $(BASEFILES) p_refine.c p_update.c check_po_arcs.c queue.c list.c -lm

all:	csa_s_tpo csa_q_tpo csa_s_rtpo csa_q_rtpo csa_s_spo csa_q_spo csa_s_rspo csa_q_rspo csa_s_tpo_pr csa_q_tpo_pr csa_s_rtpo_pr csa_q_rtpo_pr csa_s_spo_pr csa_q_spo_pr csa_s_rspo_pr csa_q_rspo_pr csa_s_tpo_pu csa_q_tpo_pu csa_s_rtpo_pu csa_q_rtpo_pu csa_s_spo_pu csa_q_spo_pu csa_s_tpo_pr_pu csa_q_tpo_pr_pu csa_s_rtpo_pr_pu csa_q_rtpo_pr_pu csa_s_spo_pr_pu csa_q_spo_pr_pu

clean: 
	-rm -f csa_s_tpo csa_q_tpo csa_s_rtpo csa_q_rtpo csa_s_spo csa_q_spo csa_s_rspo csa_q_rspo csa_s_tpo_pr csa_q_tpo_pr csa_s_rtpo_pr csa_q_rtpo_pr csa_s_spo_pr csa_q_spo_pr csa_s_rspo_pr csa_q_rspo_pr csa_s_tpo_pu csa_q_tpo_pu csa_s_rtpo_pu csa_q_rtpo_pu csa_s_spo_pu csa_q_spo_pu csa_s_tpo_pr_pu csa_q_tpo_pr_pu csa_s_rtpo_pr_pu csa_q_rtpo_pr_pu csa_s_spo_pr_pu csa_q_spo_pr_pu
