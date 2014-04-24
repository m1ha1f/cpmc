/***************************************************************************/
/*      Name:       CMF3D_ML_mex.c                               
                                                
        Performing the continuous max-flow algorithm to solve the
        3D continuous-cut problem with multiple labels (Potts Model)
 
        Usage: [u, erriter, i, timet] = CMF3D_ML_mex(penalty, C_t, para);
 
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
                    para[0,1,2]: rows, cols, heights of the given 3D volume data
                    para[3]: the number of labels or regions
                    para[4]: the maximum iteration number
                    para[5]: the error bound for convergence
                    para[6]: cc for the step-size of augmented Lagrangian method
                    para[7]: the step-size for the graident-projection step to the
                           total-variation function. Its optimal range is [0.06, 0.12].
 
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
 
               >> mex CMF3D_ML_mex.c
 
           Example:
 
               >> [u, erriter, i, timet] = CMF3D_ML_mex(single(penalty), single(Ct), single(para));

               >> You can choose a 2D slice to demonstrate the final result.
 
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




/*  compilation command (under matlab): mex CMF3D_ML_mex.c  */

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
    float   *pfbx, *pfby, *pfbz, *pfps, *pfpt, *pfgk, *pfft, *pfdiv;
    float   fLambda, fError, cc, steps;
    float   fpt, fps, *tt;
    int     *punum, iNy, iNx, iNz, iLab, iNdim, iDim[4], ix, iy, iz, id, i,ik, iNI;
    int     iNbIters, szImg, S3D, S2D;
    int index, index1, indz, indd, indy;
    time_t  start_time, end_time;
    
  
    /* Inputs */
    pfpenalty = mxGetData(pmxIn[0]);        /* penalty parameters */
    pfCt = mxGetData(pmxIn[1]); /* bound of sink flows */
    pfVecParameters = mxGetData(pmxIn[2]); /* Vector of parameters */
    
    
    /* 
     *pfVecParameters Setting
     * [0] : number of columns 
     * [1] : number of rows
     * [2] : number of heights
     * [3] : number of labels
     * [4] : the maximum iteration number
     * [5] : error criterion
     * [6] : cc for the step-size of ALM
     * [7] : steps for the step-size of projected-gradient of p
     */
    
    /* Size */
    iNy = (int) pfVecParameters[0];
    iNx = (int) pfVecParameters[1];
    iNz = (int) pfVecParameters[2];
    iLab = (int) pfVecParameters[3];
    S2D = iNy * iNx;
    S3D = S2D * iNz;
    szImg = S3D * iLab;
    
    iNbIters = (int) pfVecParameters[4]; /* total number of iterations */
    fError = (float) pfVecParameters[5]; /* error criterion */
    cc = (float) pfVecParameters[6]; /* cc for ALM */
    steps = (float) pfVecParameters[7]; /* steps for each iteration */
    
   
    /* Outputs */
    /* outputs the computed u(x)  */
    iNdim = 4;
    iDim[0] = iNy;
    iDim[1] = iNx;
    iDim[2] = iNz;
    iDim[3] = iLab; 
    
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
    pfbx = (float *) calloc( (unsigned)(iNy*(iNx+1)*iNz*iLab), sizeof(float) );
    if (!pfbx)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for py1 */
    pfby = (float *) calloc( (unsigned)((iNy+1)*iNx*iNz*iLab), sizeof(float) );
    if (!pfby)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for py1 */
    pfbz = (float *) calloc( (unsigned)(iNy*iNx*(iNz+1)*iLab), sizeof(float) );
    if (!pfbz)
        mexPrintf("Memory allocation failure\n");
   
     /* allocate the memory for ps */
    pfdiv = (float *) calloc( (unsigned)(iNy*iNx*iNz*iLab), sizeof(float) );
    if (!pfdiv)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for ps */
    pfps = (float *) calloc( (unsigned)(iNy*iNx*iNz), sizeof(float) );
    if (!pfps)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for pt */
    pfpt = (float *) calloc( (unsigned)(iNy*iNx*iNz*iLab), sizeof(float) );
    if (!pfpt)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for gk */
    pfgk = (float *) calloc( (unsigned)(iNy*iNx*iNz), sizeof(float) );
    if (!pfgk)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for ft */
    pfft = (float *) calloc( (unsigned)(iLab), sizeof(float) );
    if (!pfft)
        mexPrintf("Memory allocation failure\n");
    
    /* Preprocessing initial values */
    for (iz=0; iz< iNz; iz++){
        indz = iz*S2D;
        for (ix=0; ix < iNx; ix++){
            indy = ix*iNy + indz;
            for (iy=0; iy < iNy; iy++){
                index = indy + iy;

                fpt = pfCt[index];
                ik = 0;
                
                for (id = 1; id < iLab; id++){
                    index1 = index + id*S3D;
                    if (fpt >= pfCt[index1]){
                        fpt = pfCt[index1];
                        ik = id;
                    }
                }
                    
                pfps[index] = fpt;
                pfu[index+ik*S3D] = 1;
                
                for (id = 0; id < iLab; id++){        
                    pfpt[index+id*S3D] = fpt;
                }
            }
        }          
    }
    
    /*  Main iterations */
    
    start_time = clock();
    
    iNI = 0;
    
    while( iNI<iNbIters ) 
    { 
        
        /* 
         * update the flow fields p(x,i=1...nlab) and pt(x,i=1...nlab) 
         * at each layer 
         */
        
        for (id = 0; id < iLab; id++){
            
            indd = id*S3D;
            
            /* 
             * update the spatial flow field p(x,id) = (bx(x,id),by(x,id),bz(x,id)) 
             * bx, by, bz define the three components of the vector field p(x,id)
             */
            
            for (iz=0; iz< iNz; iz++){
                indz = iz*S2D;
                for (ix=0; ix< iNx; ix++){
                    indy = ix*iNy+indz;
                    for (iy = 0; iy < iNy; iy++){
                        index = indy + iy;
                        index1 = indd + index;
                        
                        pfgk[index] = pfdiv[index1] - (pfps[index] 
                                - pfpt[index1] + pfu[index1]/cc);
                    }
                }
            }
            
            for (iz=0; iz< iNz; iz++){
                indz = iz*S2D;
                for (ix=0; ix< iNx-1; ix++){
                    indy = (ix+1)*iNy+indz; 
                    for (iy=0; iy< iNy; iy++){
                        index = indy + iy; 
                        index1 = index + indd;
                        
                        pfbx[index1] = steps*(pfgk[index] - pfgk[index-iNy]) + 
                                pfbx[index1];
                    }
                }
            }

            
            for (iz=0; iz< iNz; iz++){
                indz = iz*S2D;
                for (ix=0; ix< iNx; ix++){
                    indy = ix*iNy+indz; 
                    for(iy = 0; iy < iNy-1; iy ++){
                        index = (iy+1) + indy;
                        index1 = index + indd;
                        
                        pfby[index1] = steps*(pfgk[index] - pfgk[index-1]) + 
                                pfby[index1];
                    }
                }
           }
            
            for (iz=0; iz< iNz-1; iz++){
                indz = (iz+1)*S2D;
                for (ix=0; ix< iNx; ix++){
                    indy = ix*iNy + indz; 
                    for(iy = 0; iy < iNy; iy ++){
                        index = indy + iy;
                        index1 = index + indd;
                        
                        pfbz[index1] = steps*(pfgk[index] - pfgk[index-S2D]) + 
                                pfbz[index1];
                    }
                }
           }

            /* projection step */
            
            for (iz=0; iz< iNz; iz++){
                indz = iz*S2D;
                for (ix=0; ix< iNx; ix++){
                    indy = ix*iNy + indz;
                    for (iy = 0; iy < iNy; iy++){
                        index = iy + indy;
                        index1 = index + indd;
                        
                        fpt = SQRT((SQR(pfbx[index1+iNy]) + SQR(pfbx[index1]) + 
                            SQR(pfby[index1+1]) + SQR(pfby[index1]) + SQR(pfbz[index1+S2D]) +
                            SQR(pfbz[index1]))*0.5);

                        if (fpt > pfpenalty[index])
                            pfgk[index] = pfpenalty[index] / fpt;
                        else
                            pfgk[index] = 1;
                    }
                }
            }
            
            /* update the component bx(x,id) */
            
            for (iz=0; iz< iNz; iz++){
                indz = iz*S2D;
                for (ix=0; ix< iNx-1; ix++){
                    indy = (ix+1)*iNy + indz;        
                    for (iy=0; iy< iNy; iy++){
                        index = iy + indy;
                        index1 = index + indd;
                        
                        pfbx[index1] = (pfgk[index] + pfgk[index-iNy])*0.5*pfbx[index1];
                    }
                }
            }

            /* update the component by(x,id) */
            
            for (iz=0; iz< iNz; iz++){
                indz = iz*S2D;
                for (ix=0; ix< iNx; ix++){
                    indy = ix*iNy + indz; 
                    for (iy=0; iy<iNy-1; iy++){
                       index = (iy+1) + indy;
                       index1 = index + indd;
                       
                       pfby[index1] = 0.5*(pfgk[index-1] + pfgk[index])*pfby[index1];
                    }
                }
            }      
            
            /* update the component bz(x,id) */
            
            for (iz=0; iz< iNz-1; iz++){
                indz = (iz+1)*S2D;
                for (ix=0; ix< iNx; ix++){
                    indy = ix*iNy + indz; 
                    for (iy=0; iy<iNy; iy++){
                       index = iy + indy;
                       index1 = index + indd;
                       
                       pfbz[index1] = 0.5*(pfgk[index-S2D] + pfgk[index])*pfbz[index1];
                    }
                }
            } 
            
            /* update the sink flow field pt(x,id)  */
            
            for (iz=0; iz< iNz; iz++){
                indz = iz*S2D;
                for (ix=0; ix< iNx; ix++){
                    indy = ix*iNy+indz;
                    for (iy = 0; iy < iNy; iy++){
                        index = iy + indy;
                        index1 = index + indd;
                        
                        /* update the divergence field div(x,id)  */
                        
                        pfdiv[index1] = pfbx[index1+iNy] - pfbx[index1] + pfby[index1+1] - 
                                pfby[index1] + pfbz[index1+S2D] - pfbz[index1];
                        
                        fpt = pfps[index] + pfu[index1]/cc - pfdiv[index1];
                        
                        pfpt[index1] = MIN(fpt , pfCt[index1]);
                    }
                }
            }
        }
        
        /* update the source flow field ps(x)  */
        
        fps = 0;
        
        for (iz=0; iz< iNz; iz++){
            indz = iz*S2D;
            for (ix=0; ix< iNx; ix++){
                indy = ix*iNy + indz;
                for (iy = 0; iy < iNy; iy++){
                    index = iy + indy;

                    fpt = 0;
                    for (id = 0; id < iLab; id++){
                        index1 = index + id*S3D;

                        pfft[id] = pfdiv[index1] + pfpt[index1];

                        fpt += (pfft[id] -pfu[index1]/cc);  

                    }

                    pfps[index] = fpt/iLab + 1/(cc*iLab);

                    /* update the multipliers */
                    
                    for (id = 0; id < iLab; id++){
                        fpt = cc*(pfft[id] - pfps[index]);
                        pfu[index+id*S3D] -= fpt;
                        fps += ABS(fpt);
                     }
                }
            }
        }

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
    free( (float *) pfbz );
    free( (float *) pfdiv );
    free( (float *) pfps );
    free( (float *) pfpt );
    free( (float *) pfgk );
    free( (float *) pfft );
    
    tt[0] = difftime(end_time, start_time)/1000000;
    
    mexPrintf("\nComputing Time for max-flow = %.4f sec\n \n",tt[0]);
    
    
}

/**********************************************/
/************** END MAIN FUNCTION *************/
/**********************************************/
