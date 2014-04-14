/*
 * Written by Andreas Mueller
 *
 */

/*
 * intersection = OVERLAP_CARE(GT, MASK, CARE)
 * Calculates the ratio between the intersection of (GT union Mask) to (GT intersect Mask)
 * where both are restricted to CARE. All arguments have to be 2d logical arrays of the same size.
 *
 */
 #include "mex.h"
void mexFunction(int nOut, mxArray *pOut[],
		 int nIn, const mxArray *pIn[])
{
  int rows, cols, len, i, j;
  mxLogical *indataGT, *indataMask, *indataCare;
  int intersect=0;
  int all=0;
  double *outdata;

  if((nIn != 3) || (nOut != 1))
    mexErrMsgTxt("Usage: intersection = overlap_care(GT, mask, care)");
  for (i=0;i<3;i++)
    if (!mxIsLogical(pIn[i]) || mxGetNumberOfDimensions(pIn[i]) != 2) {
            mexErrMsgTxt("Usage: th argument must be a logical matrix");
        }
  
  rows = mxGetM(pIn[0]);
  cols = mxGetN(pIn[0]);
  len=rows*cols;
  indataGT = (bool*) mxGetData(pIn[0]);
  indataMask = (bool*) mxGetData(pIn[1]);
  indataCare = (bool*) mxGetData(pIn[2]);
  
  
  for(i=0;i<len; i++){
        intersect+=indataGT[i] && indataMask[i] && indataCare[i];
        all+=(indataGT[i] || indataMask[i]) && indataCare[i];
      }
  pOut[0]=mxCreateDoubleMatrix(1,1, mxREAL);
  outdata=mxGetPr(pOut[0]);
  outdata[0]=((double)intersect)/all;

}
