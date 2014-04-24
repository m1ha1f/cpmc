/***************************************************************************/
/*      Name:       CMF_Mex.c                               
                                                
        Performing the continuous max-flow algorithm to solve the
        continuous min-cut problem in 3D
 
        Usage: [u, erriter, i, timet] = CMF3D_Mex(penalty, C_s, C_t, para);
 
        Inputs (penalty, C_s, C_t, para): 
 
               - penalty: point to the edge-weight penalty parameters to
                          total-variation function.
 
                 For the case without incorporating image-edge weights, 
                 penalty is given by the constant everywhere. For the case 
                 with image-edge weights, penalty is given by the pixelwise 
                 weight function:
 
                 for example, penalty(x) = b/(1 + a*| grad f(x)|) where b,a > 0 .
   
               - C_s: point to the capacities of source flows ps
 
               - C_t: point to the capacities of sink flows pt
 
               - para: a sequence of parameters for the algorithm
                    para[0,1,2]: rows, cols, heights of the given volume
                    para[3]: the maximum iteration number
                    para[4]: the error bound for convergence
                    para[5]: cc for the step-size of augmented Lagrangian method
                    para[6]: the step-size for the graident-projection step to the
                           total-variation function. Its optimal range is [0.1, 0.17].
 
        Outputs (u, erriter, i, timet):
 
              - u: the final results u(x) in [0,1]. As the following paper,
                   the global binary result can be available by threshholding u
                  by any constant alpha in (0,1):

                  Nikolova, M.; Esedoglu, S.; Chan, T. F. 
                  Algorithms for Finding Global Minimizers of Image Segmentation and Denoising Models 
                  SIAM J. App. Math., 2006, 66, 1632-1648

               - erriter: it returns the error evaluation of each iteration,
                  i.e. it shows the convergence rate. One can check the algorithm
                  performance.

               - i: gives the total number of iterations, when the algorithm converges.

               - timet: gives the total computation time.

           Compile:
 
               >> mex CMF3D_Mex.c
 
           Example:
 
               >> [u, erriter, i, timet] = CMF3D_Mex(single(penalty), single(Cs), single(Ct), single(para));

               >> us = max(u, beta);  % where beta in (0,1)

               >> figure, loglog(erriter,'DisplayName','erriterN');figure(gcf)
 
               >> isosurface(u,0.5), axis([1 rows 1 cols 1 heights]), daspect([1 1 1]);

 
       The original continuous max-flow algorithm was proposed in the following papers:

       [1] Yuan, J.; Bae, E.;  Tai, X.-C. 
           A Study on Continuous Max-Flow and Min-Cut Approaches 
           CVPR, 2010

       [2] Yuan, J.; Bae, E.; Tai, X.-C.; Boycov, Y.
           A study on continuous max-flow and min-cut approaches. Part I: Binary labeling
           UCLA CAM, Technical Report 10-61, 2010

       The mimetic finite-difference discretization method was proposed for 
       the total-variation function in the paper:

       [1] Yuan, J.; Schn{\"o}rr, C.; Steidl, G.
           Simultaneous Optical Flow Estimation and Decomposition
           SIAM J.~Scientific Computing, 2007, vol. 29, page 2283-2304, number 6

       This software can be used only for research purposes, you should cite ALL of
       the aforementioned papers in any resulting publication.

       Please email Jing Yuan (cn.yuanjing@gmail.com) for any questions, suggestions and bug reports

       The Software is provided "as is", without warranty of any kind.


                           Version 1.0
               https://sites.google.com/site/wwwjingyuan/       

               Copyright 2011 Jing Yuan (cn.yuanjing@gmail.com)      

*/
/***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mex.h>
#include <math.h>
#include <time.h>

#define YES 0
#define NO 1

#define PI 3.1415926

#define MAX(a,b) ( a > b ? a : b )
#define MIN(a,b) ( a <= b ? a : b )
#define SIGN(x) ( x >= 0.0 ? 1.0 : -1.0 )
#define ABS(x) ( (x) > 0.0 ? x : -(x) )
#define SQR(x) (x)*(x)

#ifndef HAVE_RINT 
#define rint(A) floor((A)+(((A) < 0)? -0.5 : 0.5)) 
#endif

float SQRT(float number) {
    long i;
    float x, y;
    const float f = 1.5F;
    
    x = number * 0.5F;
    y  = number;
    i  = * ( long * ) &y;
    i  = 0x5f3759df - ( i >> 1 );
    y  = * ( float * ) &i;
    y  = y * ( f - ( x * y * y ) );
    y  = y * ( f - ( x * y * y ) );
    return number * y;
}


/**********************************************/
/************** MAIN FUNCTION *****************/
/**********************************************/

/****************************************/
extern void mexFunction(int iNbOut, mxArray *pmxOut[],
int iNbIn, const mxArray *pmxIn[])
{
    
  /* iNbOut: number of outputs
     pmxOut: array of pointers to output arguments */
    
  /* iNbIn: number of inputs
     pmxIn: array of pointers to input arguments */
    
    
    float   *pfpenalty, *pfu, *pfCs, *pfCt, *pfcvg, *pfVecParameters;
    float   *pfbx, *pfby, *pfbz, *pfps, *pfpt, *pfgk, *tt, *pfdv;
    float   fError, cc, steps;
    float   fpt, fps;
    int     *punum, iNy, iNx, iNz, iNdim, iDim[3], ix, iy, iz, i, iNI;
    int     iNbIters, szImg, index, S2D, idx, idz;
    time_t  start_time, end_time;
    

    /* Inputs */
    pfpenalty = mxGetData(pmxIn[0]); /* Given penalty */
    pfCs = mxGetData(pmxIn[1]); /* bound of source flows */
    pfCt = mxGetData(pmxIn[2]); /* bound of sink flows */
    pfVecParameters = mxGetData(pmxIn[3]); /* Vector of parameters */
    
    /* 
     *pfVecParameters Setting
     * [0] : number of columns 
     * [1] : number of rows
     * [2] : number of heights
     * [3] : the maximum iteration number
     * [4] : error criterion
     * [5] : cc for the step-size of ALM
     * [6] : steps for the step-size of projected-gradient of p
     */
    
    /* Size */
    iNy = (int) pfVecParameters[0];
    iNx = (int) pfVecParameters[1];
    iNz = (int) pfVecParameters[2];
    S2D = iNy*iNx;
    szImg = iNy*iNx*iNz;
    
    
    /* Choice of region segmentation model */
    iNbIters = (int) pfVecParameters[3]; /* total number of iterations */
    fError = (float) pfVecParameters[4]; /* error criterion */
    cc = (float) pfVecParameters[5]; /* cc for ALM */
    steps = (float) pfVecParameters[6]; /* steps for each iteration */
    
   
    /* Outputs */
    /* outputs the computed u(x)  */
    iNdim = 3;
    iDim[0] = iNy;
    iDim[1] = iNx;
    iDim[2] = iNz;
    
    pmxOut[0] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    pfu = mxGetData(pmxOut[0]);
    
   
    /* outputs the convergence rate  */
    iNdim = 2;
    iDim[0] = 1;
    iDim[1] = iNbIters;
    pmxOut[1] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    pfcvg = mxGetData(pmxOut[1]);
    
    /* outputs the iteration number  */
    iNdim = 2;
    iDim[0] = 1;
    iDim[1] = 1;
    pmxOut[2] = mxCreateNumericArray(iNdim,(const int*)iDim,mxUINT16_CLASS,mxREAL);
    punum = mxGetData(pmxOut[2]);
    
    /* outputs the computation time  */
    iNdim = 2;
    iDim[0] = 1;
    iDim[1] = 1;
    pmxOut[3] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    tt = mxGetData(pmxOut[3]);
    
    /* Memory allocation */
    
    /* allocate the memory for px1 */
    pfbx = (float *) calloc( (unsigned)(iNy*(iNx+1)*iNz), sizeof(float) );
    if (!pfbx)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for py1 */
    pfby = (float *) calloc( (unsigned)((iNy+1)*iNx*iNz), sizeof(float) );
    if (!pfby)
        mexPrintf("Memory allocation failure\n");

    /* allocate the memory for pz1 */
    pfbz = (float *) calloc( (unsigned)(iNy*iNx*(iNz+1)), sizeof(float) );
    if (!pfbz)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for ps */
    pfps = (float *) calloc( (unsigned)(iNy*iNx*iNz), sizeof(float) );
    if (!pfps)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for pt */
    pfpt = (float *) calloc( (unsigned)(iNy*iNx*iNz), sizeof(float) );
    if (!pfpt)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for gk */
    pfgk = (float *) calloc( (unsigned)(iNy*iNx*iNz), sizeof(float) );
    if (!pfgk)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for div */
    pfdv = (float *) calloc( (unsigned)(iNy*iNx*iNz), sizeof(float) );
    if (!pfdv)
        mexPrintf("Memory allocation failure\n");
    
    /* Preprocessing initial values */
    for (iz=0; iz < iNz; iz++){
        idz = iz*S2D;
        for (ix=0; ix < iNx; ix++){
            idx = idz + ix*iNy;
            for (iy=0; iy < iNy; iy++){
                
                index = idx + iy;
                
                if (pfCs[index] < pfCt[index]){
                    pfps[index] = pfCs[index];
                    pfpt[index] = pfCs[index];
                    pfdv[index] = pfbx[index+iNy] - pfbx[index] 
                                + pfby[index+1] - pfby[index]
                                + pfbz[index+S2D] - pfbz[index];
                }
                else{
                    pfu[index] = 1;
                    pfps[index] = pfCt[index];
                    pfpt[index] = pfCt[index];
                    pfdv[index] = pfbx[index+iNy] - pfbx[index] 
                                + pfby[index+1] - pfby[index]
                                + pfbz[index+S2D] - pfbz[index];
                }
            }
        }
    }
 
     /*  Main iterations */
    
    iNI = 0;
     
    start_time = clock();
    
    while( iNI<iNbIters ) 
    { 
        /* update p */
        
        for (iz=0; iz < iNz; iz++){
            idz = iz*S2D;
            for (ix=0; ix < iNx; ix++){
                idx = idz + ix*iNy;
                for (iy=0; iy < iNy; iy++){
                    index = idx + iy;
                    pfgk[index] = pfdv[index] - (pfps[index] 
                                - pfpt[index] + pfu[index]/cc);
                }
            }
        }   
        
        /* update px*/ 
        
        for (iz=0; iz < iNz; iz++){
            idz = iz*S2D;
            for (ix=0; ix < iNx-1; ix++){
                idx = idz + (ix+1)*iNy;
                for (iy=0; iy < iNy; iy++){
                    index = idx + iy;
                    pfbx[index] = steps*(pfgk[index] - pfgk[index-iNy]) + pfbx[index];
                }
            }
        }
        
        /* update py */
        
        for (iz=0; iz < iNz; iz++){
            idz = iz*S2D;
            for (ix=0; ix < iNx; ix++){
                idx = idz + ix*iNy;
                for (iy=0; iy < iNy-1; iy++){
                    index = idx + iy + 1;
                    pfby[index] = steps*(pfgk[index] - pfgk[index-1]) + pfby[index];
                }
            }
        }
        
        /* update pz */
        
        for (iz=0; iz < iNz-1; iz++){
            idz = (iz+1)*S2D;
            for (ix=0; ix < iNx; ix++){
                idx = idz + ix*iNy;
                for (iy=0; iy < iNy; iy++){
                    index = idx + iy;
                    pfbz[index] = steps*(pfgk[index] - pfgk[index-S2D]) + pfbz[index];
                }
            }
        }
        
        /* projection step */
        
        for (iz=0; iz < iNz; iz++){
            idz = iz*S2D;
            for (ix=0; ix < iNx; ix++){
                idx = idz + ix*iNy;
                for (iy=0; iy < iNy; iy++){
                    index = idx + iy;
                    fpt = SQRT((SQR(pfbx[index]) + SQR(pfbx[index+iNy]) 
                      + SQR(pfby[index]) + SQR(pfby[index+1]) 
                      + SQR(pfbz[index]) + SQR(pfbz[index+S2D]))*0.5);
                
                    if (fpt > pfpenalty[index])
                        fpt = fpt / pfpenalty[index];
                    else
                        fpt = 1;

                    pfgk[index] = 1/fpt;
                }
            }
        }
        
        for (iz=0; iz < iNz; iz++){
            idz = iz*S2D;
            for (ix=0; ix < iNx-1; ix++){
                idx = idz + (ix+1)*iNy;
                for (iy=0; iy < iNy; iy++){
                    index = idx + iy;
                    pfbx[index] = (pfgk[index] + pfgk[index-iNy])*0.5*pfbx[index];
                }
            }
        }
        
        for (iz=0; iz < iNz; iz++){
            idz = iz*S2D;
            for (ix=0; ix < iNx; ix++){
                idx = idz + ix*iNy;
                for (iy=0; iy < iNy-1; iy++){
                    index = idx + iy + 1;
                    pfby[index] = 0.5*(pfgk[index] + pfgk[index-1])*pfby[index];
                }
            }
        }
        
        for (iz=0; iz < iNz-1; iz++){
            idz = (iz+1)*S2D;
            for (ix=0; ix < iNx; ix++){
                idx = idz + ix*iNy;
                for (iy=0; iy < iNy; iy++){
                    index = idx + iy;
                    pfbz[index] = 0.5*(pfgk[index] + pfgk[index-S2D])*pfbz[index];
                }
            }
        }
        
        fps = 0;
        
        for (iz=0; iz < iNz; iz++){
            idz = iz*S2D;
            for (ix=0; ix < iNx; ix++){
                idx = idz + ix*iNy;
                for (iy=0; iy < iNy; iy++){
                    
                    float fpt; 
                    
                    index = idx + iy;
            
                    /* compute the divergence */

                    pfdv[index] = pfbx[index+iNy] - pfbx[index] 
                         + pfby[index+1] - pfby[index]
                         + pfbz[index+S2D] - pfbz[index];

                    /* update ps */

                    fpt = pfpt[index] - pfu[index]/cc + pfdv[index] + 1/cc;
                    pfps[index] = MIN(fpt, pfCs[index]);

                    /* update pt */

                    fpt = pfps[index] + pfu[index]/cc - pfdv[index];
                    pfpt[index]  = MIN(fpt , pfCt[index]);

                    /* update multipliers */

                    fpt = cc*(pfpt[index] + pfdv[index] - pfps[index]);
                    fps += ABS(fpt);

                    pfu[index] -= fpt;
                }
            }
        }
        
         pfcvg[iNI] = fps / szImg;
        
        if (pfcvg[iNI] <= fError)
            break;
                
        iNI ++;
 
    }   
   
    mexPrintf("number of iterations = %i\n",iNI);
    
    /* Outputs (see above) */
    
    punum[0] = iNI;
    
    /* Free memory */
    
    free( (float *) pfbx );
    free( (float *) pfby );
    free( (float *) pfbz );
    free( (float *) pfps );
    free( (float *) pfpt );
    free( (float *) pfgk );
    free( (float *) pfdv );
    
    end_time = clock();
    tt[0] = difftime(end_time,start_time)/1000;
    
    mexPrintf("\nComputing Time for max-flow = %.4f sec\n \n",tt[0]);
    
}
/****************************************/

