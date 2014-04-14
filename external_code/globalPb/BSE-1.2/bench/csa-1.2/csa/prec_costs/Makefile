CXX = gcc
BASEFILES=main.c refine.c update_epsilon.c parse.c stack.c timer.c debug.c
HEADERS=csa_types.h csa_defs.h
CFLAGS=-O3

all: links csa

links: main.c parse.c stack.c timer.c list.c queue.c csa_defs.h
	ln -sf ../main.c main.c
	ln -sf ../parse.c parse.c
	ln -sf ../stack.c stack.c
	ln -sf ../timer.c timer.c
	ln -sf ../list.c list.c
	ln -sf ../queue.c queue.c
	ln -sf ../csa_defs.h csa_defs.h
	touch links

csa: $(BASEFILES) $(HEADERS) links
	${CXX} $(CFLAGS) -DQUICK_MIN -DSAVE_RESULT -DMIN_COST -o $@ $(BASEFILES)

install: csa
	cp $< ../../../../

clean:
	rm -f csa
