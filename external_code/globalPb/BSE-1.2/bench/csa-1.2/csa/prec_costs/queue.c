#include	<stdio.h>
#include	<malloc.h>
#include	"csa_types.h"

extern	char	*nomem_msg;

queue	q_create(size)

unsigned	size;

{
queue	q;
void	exit();

q = (queue) malloc(sizeof(struct queue_st));
if (q == NULL)
  {
  (void) fprintf(stderr,nomem_msg);
  exit(9);
  }
q->storage = (char **) malloc(sizeof(lhs_ptr) * (size + 1));
if (q->storage == NULL)
  {
  (void) fprintf(stderr,nomem_msg);
  exit(9);
  }
q->end = q->storage + size;
q->tail = q->head = q->storage;
q->max_size = size;
return(q);
}

char	*deq(q)

queue	q;

{
char	*p;

p = *(q->head);
if (q->head == q->end) q->head = q->storage;
else q->head++;
return(p);
}
