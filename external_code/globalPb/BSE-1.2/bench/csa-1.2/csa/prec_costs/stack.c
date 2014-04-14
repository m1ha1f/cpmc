#include	<stdio.h>
#include	<malloc.h>
#include	"csa_types.h"

char	*nomem_msg = "Insufficient memory.\n";

void	st_reset(s)

stack	s;

{
s->top = s->bottom;
}

char	*st_pop(s)

stack	s;

{
s->top--;
return(*(s->top));
}

stack	st_create(size)

unsigned	size;

{
stack	s;
void	exit();
  
s = (stack) malloc(sizeof(struct stack_st));

if (s == NULL)
  {
  (void) fprintf(stderr,nomem_msg);
  exit(9);
  }
s->bottom = (char **) malloc(size * sizeof(char *));
if (s->bottom == NULL)
  {
  (void) fprintf(stderr,nomem_msg);
  exit(9);
  }
s->top = s->bottom;
return(s);
}
