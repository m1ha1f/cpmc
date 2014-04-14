/*---
function [i,j] = cartprod(set1,set2) or
  
Input:    
    set1, set2 = column vectors
Output:
    [i,j] = column vectors of UINT32
        
--*/


# include "mex.h"
# include "math.h"

void mexFunction(
    int nargout,
    mxArray *out[],
    int nargin,
    const mxArray *in[]
) {
  /* declare variables */
  int nr, nc, np, nb, total, to_samp1, to_samp2;    
  double *sz, *set1, *set2, *qi, *qj;
  int nr1, nr2, c;
  unsigned long i, j, count;

  /* check argument */
  if (nargin < 2) {
      mexErrMsgTxt("Two input arguments required");
  }
  if (nargout > 2) {
      mexErrMsgTxt("Too many output arguments.");
  }

  set1 = mxGetData(in[0]);
  nr1 = mxGetM(in[0]);   
  c = mxGetN(in[0]);

  if (c > nr1) {
      mexErrMsgTxt("Expecting column vector!");
  }

  set2 = mxGetData(in[1]);
  nr2 = mxGetM(in[1]);   
  c = mxGetN(in[1]);    

  if (c > nr2) {
      mexErrMsgTxt("Expecting column vector!");
  }
    
  out[0] = mxCreateDoubleMatrix(nr1*nr2, 1, mxREAL);
  out[1] = mxCreateDoubleMatrix(nr1*nr2, 1, mxREAL);
  
  qi = mxGetData(out[0]);
  qj = mxGetData(out[1]);
  if (out[0]==NULL || out[1]==NULL) {
    mexErrMsgTxt("Not enough space for the output matrix.");
  }

  count = 0;
  /* computation */     
  for (i=0; i<nr1; i++) {
    for (j=0; j<nr2; j++) {
      qi[count] = set1[i];
      qj[count] = set2[j];
      count++;
    }
  }
}
