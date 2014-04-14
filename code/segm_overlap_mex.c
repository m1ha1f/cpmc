/* Copyright (C) 2010 Joao Carreira

 This code is part of the extended implementation of the paper:
 
 J. Carreira, C. Sminchisescu, Constrained Parametric Min-Cuts for Automatic Object Segmentation, IEEE CVPR 2010
 */

/*---
function [o] = segm_overlap_mex(segms)
Input:
    segms - binary matrix, with segms in columns  

Output:
    o = symmetric overlap matrix
--*/

#include "mex.h"
#include "math.h"
#include "time.h"

extern void overlap(unsigned int *intersections, unsigned int *reunions, mxLogical *segms, int nc, int nr);

void mexFunction(
    int nargout,
    mxArray *out[],
    int nargin,
    const mxArray *in[]) {

    /* declare variables */
    int nr, nc;
    register unsigned int index_ij, index_ik, index_jk;
	register unsigned int i, j, k;
    mxLogical *segms;
    unsigned int *intersections, *reunions;
    float *o;
    
    /* check argument */
    if (nargin<1) {
        mexErrMsgTxt("One input argument required ( the segments in column form )");
    }
    if (nargout>1) {
        mexErrMsgTxt("Too many output arguments");
    }

    nr = mxGetM(in[0]);
    nc = mxGetN(in[0]);

    if (!mxIsLogical(in[0]) || mxGetNumberOfDimensions(in[0]) != 2) {
        mexErrMsgTxt("Usage: segms must be a logical matrix");
    }

    segms = (bool *) mxGetData(in[0]);
    intersections = (unsigned int *)malloc(nc* nc*sizeof(unsigned int));
    reunions = (unsigned int *)malloc(nc* nc*sizeof(unsigned int));        
    
    if (intersections==NULL || reunions == NULL) {
        mexErrMsgTxt("Not enough memory");
    }
    
    out[0] = mxCreateNumericMatrix(nc,nc,mxSINGLE_CLASS, (mxComplexity) 0);
    if (out[0]==NULL) {
	    mexErrMsgTxt("Not enough memory for the output matrix");
    }
    o = (float *) mxGetPr(out[0]);
	


  	/*clock_t begin=clock();*/
    /* intersections and reunions */
/*
    #pragma omp parallel for
    for(i=0; i<nc; i++) {
        for(j=i+1; j<nc; j++) {
            index_ij = i*nc + j;
            intersections[index_ij] = 0;
            reunions[index_ij] = 0;            
            
            for( k=0; k<nr; k++) {
                index_ik = i*nr + k;
                index_jk = j*nr + k;
                intersections[index_ij] = intersections[index_ij] + (segms[index_ik] * segms[index_jk]);
                reunions[index_ij] = reunions[index_ij] + (segms[index_ik] | segms[index_jk]);
            }
        }
    }
*/
    overlap(intersections, reunions, segms, nc, nr);

   /*clock_t end=clock();
    double diffms=((end - begin)*1000)/CLOCKS_PER_SEC;
    mexPrintf("time taken in main loop: %f\n", diffms);*/
    
    /* computation of overlap, for lower triangle */
    for (i=0; i<nc; i++) {
        for(j=i+1; j<nc; j++) {             
            index_ij = i*nc+j;
            o[index_ij] = (float) intersections[index_ij] / reunions[index_ij];            
        }
    }
        
    /* copy lower triangle to upper triangle; matrix is symmetric */
    for (i=0; i<nc; i++) {      
        /*mexPrintf("i: %d\n",i);*/
        for(j=i+1; j<nc; j++) {   
            /*mexPrintf("j: %d\n",j);
            mexPrintf("index: %d\n", j*nc+i);*/
            o[j*nc+i] = o[i*nc+j];
        }
    }
    
    /* fill diagonal with ones */
    for (i=0; i<nc; i++) {
       o[i*nc+i] = 1;
    }
        
/*
    for(i=0;i<nc; i++) {
        for (j=0;j<nr; j++) {
            mexPrintf("o(i,j): %f\n", o[i*nc+j]);
        }
    }
*/

}
