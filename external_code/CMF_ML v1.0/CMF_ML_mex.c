/***************************************************************************/
/*      Name:       CMF_ML_mex.c                               
                                                
        Performing the continuous max-flow algorithm to solve the
        2D continuous-cut problem with multiple labels (Potts Model)
 
        Usage: [u, erriter, i, timet] = CMF_ML_mex(penalty, C_t, para);
 
        Inputs (penalty, C_t, para): 
 
               - penalty(x): point to the edge-weight penalty to
                          the total-variation function.
 
                 For the case without incorporating image-edge weights, 
                 penalty is given by the constant everywhere. For the case 
                 with image-edge weights, penalty is given by the pixelwise 
                 weight function:
 
                 for example, penalty(x) = b/(1 + a*| grad f(x)|) where b,a > 0.
   
               - C_t(x,i=1...nlab): point to the capacities of 
                 sink flows pt(x,i=1...nlab);
 
               - para: a sequence of parameters for the algorithm
                    para[0,1]: rows, cols of the given image
                    para[2]: the number of labels or regions
                    para[3]: the maximum iteration number
                    para[4]: the error bound for convergence
                    para[5]: cc for the step-size of augmented Lagrangian method
                    para[6]: the step-size for the graident-projection step to the
                           total-variation function. Its optimal range is [0.1, 0.17].
 
        Outputs (u, erriter, i, timet):
 
              - u: the computed continuous labeling function u(x,i=1...nlab) in [0,1]. 
                   As the following paper [2], the final label function is 
                   given by the maximum of u(x,i=1...nlab) at each x.
 
               - erriter: it returns the error evaluation of each iteration,
                  i.e. it shows the convergence rate. One can check the algorithm
                  performance.

               - i: gives the total number of iterations, when the algorithm converges.

               - timet: gives the total computation time.

           Compile:
 
               >> mex CMF_ML_mex.c
 
           Example:
 
               >> [u, erriter, i, timet] = CMF_ML_mex(single(penalty), single(Ct), single(para));

               >> To demonstrate the final result, see the related matlab file.
 
               >> figure, loglog(erriter,'DisplayName','erriter');figure(gcf)

 
       The original continuous max-flow algorithm was proposed in the following papers:

       [1] Yuan, J.; Bae, E.;  Tai, X.-C. 
           A Study on Continuous Max-Flow and Min-Cut Approaches 
           CVPR, 2010

       [2] Yuan, J.; Bae, E.; Tai, X.-C.; Boycov, Y.
           A Continuous Max-Flow Approach to Potts Model
           ECCV, 2010

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




/*  compilation command (under matlab): mex CMF_ML_mex.c  */

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
    
    
    float   *pfu, *pfCt, *pfcvg, *pfVecParameters, *pfpenalty;
    float   *pfbx, *pfby, *pfps, *pfpt, *pfgk, *pfft, *pfdv;
    float   fLambda, fError, cc, steps;
    float   fpt, fps, *tt;
    int     *punum, iNy, iNx, iLab, iNdim, iDim[4], ix, iy, iz, id, i, iNI;
    int     iNbIters, szImg, S2D, idx, index, index1, idz;
    time_t  start_time, end_time;
    
    
    /* Inputs */
    pfpenalty = mxGetData(pmxIn[0]);        /* penalty parameters */
    pfCt = mxGetData(pmxIn[1]);             /* bound of sink flows */
    pfVecParameters = mxGetData(pmxIn[2]);  /* Vector of parameters */
    
    
    /* 
     *pfVecParameters Setting
     * [0] : number of columns 
     * [1] : number of rows
     * [2] : number of labels
     * [3] : the maximum iteration number
     * [4] : error criterion
     * [5] : cc for the step-size of ALM
     * [6] : steps for the step-size of projected-gradient of p
     */
    
    /* Size */
    iNy = (int) pfVecParameters[0];
    iNx = (int) pfVecParameters[1];
    iLab = (int) pfVecParameters[2];
    szImg = iNy * iNx * iLab;
    S2D = iNy*iNx;
    
    iNbIters = (int) pfVecParameters[3]; 
    fError = (float) pfVecParameters[4]; 
    cc = (float) pfVecParameters[5]; 
    steps = (float) pfVecParameters[6]; 
    
   
    /* Outputs */
    /* outputs the computed u(x)  */
    iNdim = 3;
    iDim[0] = iNy;
    iDim[1] = iNx;
    iDim[2] = iLab; 
    
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
    
    /* allocate the memory for bx, where p(x,i=1...nlab) = (bx(x,i),by(x,i)) */

    pfbx = (float *) calloc( (unsigned)(iNy*(iNx+1)*iLab), sizeof(float) );
    if (!pfbx)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for by */
    
    pfby = (float *) calloc( (unsigned)((iNy+1)*iNx*iLab), sizeof(float) );
    if (!pfby)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for ps */
    pfps = (float *) calloc( (unsigned)(iNy*iNx), sizeof(float) );
    if (!pfps)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for pt */
    pfpt = (float *) calloc( (unsigned)(iNy*iNx*iLab), sizeof(float) );
    if (!pfpt)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for div */
    pfdv = (float *) calloc( (unsigned)(iNy*iNx*iLab), sizeof(float) );
    if (!pfdv)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for gk */
    pfgk = (float *) calloc( (unsigned)(iNy*iNx), sizeof(float) );
    if (!pfgk)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for ft */
    pfft = (float *) calloc( (unsigned)(iLab), sizeof(float) );
    if (!pfft)
        mexPrintf("Memory allocation failure\n");
    
    /* Preprocessing initial values for flows and u */
    for (ix=0; ix < iNx; ix++){
        idx = ix*iNy;
        for (iy=0; iy < iNy; iy++){
            index = idx + iy;
            
            fpt = pfCt[index];
            
            pfdv[index] = pfby[index+1] - pfby[index] 
                        + pfbx[index+iNy] - pfbx[index];
            id = 0;

            for (iz = 1; iz < iLab; iz++){
                idz = index + iz*S2D;
                
                pfdv[idz] = pfby[idz+1] - pfby[idz] 
                          + pfbx[idz+iNy] - pfbx[idz];
                
                if (fpt >= pfCt[idz]){
                    fpt = pfCt[idz];
                    id = iz;
                }
            }

            pfps[index] = fpt;
            pfu[index + id*S2D] = 1;

            for (iz = 0; iz < iLab; iz++){
                pfpt[index + iz*S2D] = fpt;
            }
        }
    }
    
    /*  Main iterations */
    
    iNI = 0;

    start_time = clock();
    
    
    while( iNI<iNbIters ) 
    { 
        
        /* 
         * update the flow fields p(x,i=1...nlab) and pt(x,i=1...nlab) 
         * at each layer 
         */
        
        for (iz = 0; iz < iLab; iz++){
        
            idz = iz*S2D;
            
            /* update the spatial flow field p(x,iz) = (bx(x,iz),by(x,iz)) */
            
            for (ix=0; ix < iNx; ix++){
                idx = ix*iNy;
                for (iy=0; iy < iNy; iy++){
                    index1 = idx + iy;
                    index = index1 + idz;
                    pfgk[index1] = pfdv[index] - (pfps[index1] 
                                    - pfpt[index] + pfu[index]/cc);
                }
            }
            
            for (ix=0; ix< iNx-1; ix++){
                idx = (ix+1)*iNy; 
                for (iy=0; iy< iNy; iy++){
                    index1 = idx + iy;
                    index = index1 + idz;
                    pfbx[index] = steps*(pfgk[index1] - pfgk[index1-iNy]) + pfbx[index];
                }
            }
            
            for(ix = 0; ix < iNx; ix ++){
                idx = ix*iNy; 
                for(iy = 0; iy < iNy-1; iy ++){
                    index1 = idx + iy + 1;
                    index = index1 + idz;
                    
                    pfby[index] = steps*(pfgk[index1] - pfgk[index1-1]) + pfby[index];
                }
            }
            
            /* projection step */
            
            for (ix=0; ix< iNx; ix++){
                idx = ix*iNy;
                for (iy=0; iy< iNy; iy++){
                    index1 = idx + iy;
                    index = index1 + idz;
                    
                    fpt = SQRT((SQR(pfbx[index+iNy]) + SQR(pfbx[index]) + 
                            SQR(pfby[index+1]) + SQR(pfby[index]))*0.5);

                    if (fpt > pfpenalty[index1])
                        fpt = fpt / pfpenalty[index1];
                    else
                        fpt = 1;

                    pfgk[index1] = 1/fpt;
                }
            }
            
            /* update the component bx(x,iz) */

            for (ix=0; ix< iNx-1; ix++){
                idx = (ix+1)*iNy;
                for (iy=0; iy< iNy; iy++){
                    index1 = idx + iy;
                    index = index1 + idz;
                    
                    pfbx[index] = (pfgk[index1] + pfgk[index1-iNy])
                                  *0.5*pfbx[index];
                }
            }
            
            /* update the component by(x,iz) */
            
            for (ix=0; ix<iNx; ix++){
                idx = ix*iNy;
                for (iy=0; iy< iNy-1; iy++){
                    index1 = idx + iy + 1;
                    index = index1 + idz;
                    
                    pfby[index] = 0.5*(pfgk[index1-1] + pfgk[index1])
                                *pfby[index];
                }
            }      
            
            /* update the sink flow field pt(x,iz)  */

            for (ix=0; ix< iNx; ix++){
                idx = ix*iNy;
                for (iy=0; iy< iNy; iy++){
                    index1 = idx + iy;
                    index = index1 + idz;
                    
                    /* update the divergence field dv(x,iz)  */
                    
                    pfdv[index] = pfby[index+1] - pfby[index] 
                          + pfbx[index+iNy] - pfbx[index];
                    
                    fpt = pfps[index1] + pfu[index]/cc - pfdv[index];
                    
                    pfpt[index] = MIN(fpt , pfCt[index]);
                }
            }
        
        }
        
        /* update the source flow field ps(x)  */
        
        fps = 0;
        
        for (ix=0; ix< iNx; ix++){
            idx = ix*iNy;
            for (iy=0; iy< iNy; iy++){
                index1 = idx + iy;
                
                fpt = 0;
                
                for (iz = 0; iz < iLab; iz++){
                    index = index1 + iz*S2D;
                    
                    pfft[iz] = pfdv[index] + pfpt[index];
                    
                    fpt += (pfft[iz] -pfu[index]/cc);  
                    
                    }
                
                pfps[index1] = fpt/iLab + 1/(cc*iLab);
                
                /* update the multipliers */
                
                for (iz = 0; iz < iLab; iz++){
                    fpt = cc*(pfft[iz] - pfps[index1]);
                    pfu[index1+iz*S2D] -= fpt;
                    fps += ABS(fpt);
                }
            }
        }
        
        /* evaluate the convergence error bound */
        
        pfcvg[iNI] = fps / szImg;
        
        if (pfcvg[iNI] <= fError)
            break;
                
        iNI ++;
     }   
    
    end_time = clock();
    
    /* Outputs (see above) */
    punum[0] = ++iNI;
    mexPrintf("number of iterations = %i\n",iNI);
    
    /* Free memory */
    free( (float *) pfbx );
    free( (float *) pfby );
    free( (float *) pfps );
    free( (float *) pfpt );
    free( (float *) pfgk );
    free( (float *) pfdv );
    free( (float *) pfft );

    tt[0] = difftime(end_time, start_time)/1000000;
   
    mexPrintf("\nComputing Time for max-flow = %.4f sec\n \n", tt[0]);
    
}

/**********************************************/
/************** END MAIN FUNCTION *************/
/**********************************************/
