#include <mex.h>
#include "chi2float.h"
/*#include "chi2float.h"*/

/*
  computes the chi??? distance between the input arguments
  d(X,Y) = sum ((X(i)-Y(i))???)/(X(i)+Y(i))
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	float *vecA, *vecB, *pA,*pB, *dist;
	int i,j,dim,ptsA, ptsB, k,kptsA,kptsB, *sym;

	if (nrhs == 0)
	{
		mexPrintf("Usage: d = chi2_mex(X,Y);\n");
		mexPrintf("where X and Y are matrices of dimension [dim,npts]\n");
		mexPrintf("\nExample\n a = rand(2,10);\n b = rand(2,20);\n d = chi2_mex(a,b);\n");
		return;
	}

	if (nrhs != 3){
		mexPrintf("three input arguments expected: A, B, and sym, sym says wether A and B are the same");
		return;
	}

	if (mxGetNumberOfDimensions(prhs[0]) != 2 || mxGetNumberOfDimensions(prhs[1]) != 2)
	{
		mexPrintf("inputs must be two dimensional");
		return;
	}
    
    
    mxClassID the_class = mxGetClassID(prhs[0]);
    
    if(the_class != mxSINGLE_CLASS) {
        mexErrMsgTxt("Histograms should have single precision!\n");
    }
            

	vecA = (float *)mxGetPr(prhs[0]);
	vecB = (float *)mxGetPr(prhs[1]);
    sym = (int *) mxGetPr(prhs[2]);
    
	ptsA = mxGetN(prhs[0]);
	ptsB = mxGetN(prhs[1]);
	dim = mxGetM(prhs[0]);

	if (dim != mxGetM(prhs[1]))
	{
		mexPrintf("Dimension mismatch");
		return;
	}

    const mwSize ndims[2] = {ptsA, ptsB};
    
    mxArray* mxdist = mxCreateNumericArray(2, ndims,the_class,mxREAL);    
    dist = (float *)mxGetData(mxdist);
	/*plhs[0] = mxCreateDoubleMatrix(ptsA,ptsB,mxREAL);*/
	/*dist = (float *)mxGetPr(plhs[0]);*/
    
	/*chi2_distance_float(dim,ptsB,vecB,ptsA,vecA,dist);    */
          
/* printf("get_num_threads: %d\n", omp_get_num_threads()); 
 printf("hello");*/
 
    if(*sym) {
        chi2sym_distance_float(dim,ptsB,vecB, dist); 
    } else {
        chi2_distance_float(dim,ptsB,vecB,ptsA,vecA,dist); 
    }
            
    plhs[0] = mxdist;
    
	return;
}
