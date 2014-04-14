#include	"csa_types.h"

extern	rhs_ptr	tail_rhs_node;

rhs_ptr	deq_list(head)

rhs_ptr	*head;

{
rhs_ptr	ans;

ans = *head;
*head = ans->next;
ans->next->prev = tail_rhs_node;
return(ans);
}
