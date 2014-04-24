/***************************************************************************/
/*      Name:       CMF_GPU.c                               
                                                
        Performing the continuous max-flow algorithm to solve the
        2D continuous min-cut problem over GPU (Nvidia CUDA based)
 
        Usage: [u, erriter, i, timet] = CMF_GPU(penalty, C_s, C_t, para);
 
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
                    para[0,1]: rows, cols of the given image
                    para[2]: the maximum iteration number
                    para[3]: the error bound for convergence
                    para[4]: cc for the step-size of augmented Lagrangian method
                    para[5]: the step-size for the graident-projection step to the
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
 
               >> nvmex -f nvmexopts.bat CMF_GPU.cu -IC:\cuda\v4.0\include -LC:\cuda\v4.0\lib\x64 -lcufft -lcudart
                
                Note: The compilation process depends on your system configuration!
                      You should change the path to your own cuda installation path.
 
           Example:
 
               >> [u, erriter, i, timet] = CMF_GPU(single(penalty), single(Cs), single(Ct), single(para));

               >> us = max(u, beta);  % where beta in (0,1)

               >> imagesc(us), colormap gray, axis image, axis off;figure(gcf)

               >> figure, loglog(erriter,'DisplayName','erriterN');figure(gcf)

 
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
#include "cuda.h"
#include "cuda_runtime_api.h"
#include "CMF_kernels.cu"

#define YES 0
#define NO 1

#define PI 3.1415926

#define MAX(a,b) ( a > b ? a : b )
#define MIN(a,b) ( a <= b ? a : b )
#define SIGN(x) ( x >= 0.0 ? 1.0 : -1.0 )
#define ABS(x) ( (x) > 0.0 ? x : -(x) )

#ifndef HAVE_RINT 
#define rint(A) floor((A)+(((A) < 0)? -0.5 : 0.5)) 
#endif


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
    float   *pfbx, *pfby, *pfps, *pfpt, *pfgk, *tt, *pfdv;
    float   fError, cc, steps, fps;
    int     *punum, iNy, iNx, iNdim, iDim[3], ix, iy, iNI;
    int iDev;
    int     iNbIters, szImg, idx, index;
    time_t  start_time, end_time;

    //    GPU Variables
    float   *pfbx_GPU, *pfby_GPU, *pfpenalty_GPU, *pfdv_GPU;
    float   *pfps_GPU, *pfpt_GPU, *pfgk_GPU, *pfu_GPU, *pfCs_GPU, *pfCt_GPU;
    float   *FPS, *FPS_GPU;
    
    
    cudaDeviceProp prop;

    cudaGetDeviceCount(&iDev);

    if ((unsigned int)iDev == 0){
        printf("There is no CUDA device found!");
        return;
    }
    else{
        printf("There are %d CUDA devices in your computer. \n", iDev);
        for(int ii = 0; ii < iDev; ii ++){
            cudaGetDeviceProperties(&prop, ii);
            printf("------ General Information for CUDA device %d ------ \n", ii);
            printf("Name:  %s \n", prop.name);
            printf("Multiprocessor count:  %d \n", prop.multiProcessorCount);
            printf("Total global memory: %ld \n", prop.totalGlobalMem);
            printf("---------------------------------------------------- \n\n");
         }
    }

    /* Inputs */
    pfpenalty = (float*) mxGetData(pmxIn[0]); /* Given penalty */
    pfCs = (float*) mxGetData(pmxIn[1]); /* bound of source flows */
    pfCt = (float*) mxGetData(pmxIn[2]); /* bound of sink flows */
    pfVecParameters = (float*) mxGetData(pmxIn[3]); /* Vector of parameters */
    
     /* 
     *pfVecParameters Setting
     * [0] : number of columns 
     * [1] : number of rows
     * [2] : the maximum iteration number
     * [3] : error criterion
     * [4] : cc for the step-size of ALM
     * [5] : steps for the step-size of projected-gradient of p
     */
   
    /* Size */
    iNy = (int) pfVecParameters[0];
    iNx = (int) pfVecParameters[1];
    szImg = iNy * iNx;
    
    /* Choice of region segmentation model */
    iNbIters = (int) pfVecParameters[2]; /* the maximum iteration number */
    fError = (float) pfVecParameters[3]; /* error bound for convergence */
    cc = (float) pfVecParameters[4]; /* the step-size for ALM */
    steps = (float) pfVecParameters[5]; /* the step-size for each projected-gradient step */
    
   printf("Initializing ................................................ \n\n");

    /* Outputs */
    /* outputs the computed u(x)  */
    iNdim = 2;
    iDim[0] = iNy;
    iDim[1] = iNx;
    
    pmxOut[0] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    pfu = (float*) mxGetData(pmxOut[0]);
    
   
    /* outputs the convergence rate  */
    iNdim = 2;
    iDim[0] = 1;
    iDim[1] = iNbIters;
    pmxOut[1] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    pfcvg = (float*) mxGetData(pmxOut[1]);
    
    /* outputs the iteration number  */
    iNdim = 2;
    iDim[0] = 1;
    iDim[1] = 1;
    pmxOut[2] = mxCreateNumericArray(iNdim,(const int*)iDim,mxUINT16_CLASS,mxREAL);
    punum = (int*) mxGetData(pmxOut[2]);
    
    /* outputs the computation time  */
    iNdim = 2;
    iDim[0] = 1;
    iDim[1] = 1;
    pmxOut[3] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    tt = (float*)mxGetData(pmxOut[3]);
    
    /* Memory allocation */
    
    /* allocate the memory for px */
    pfbx = (float *) calloc( (unsigned)(iNy*(iNx+1)), sizeof(float) );
    if (!pfbx)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for py */
    pfby = (float *) calloc( (unsigned)((iNy+1)*iNx), sizeof(float) );
    if (!pfby)
        mexPrintf("Memory allocation failure\n");

    /* allocate the memory for ps */
    pfps = (float *) calloc( (unsigned)(iNy*iNx), sizeof(float) );
    if (!pfps)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for pt */
    pfpt = (float *) calloc( (unsigned)(iNy*iNx), sizeof(float) );
    if (!pfpt)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for gk */
    pfgk = (float *) calloc( (unsigned)(iNy*iNx), sizeof(float) );
    if (!pfgk)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for div */
    pfdv = (float *) calloc( (unsigned)(iNy*iNx), sizeof(float) );
    if (!pfdv)
        mexPrintf("Memory allocation failure\n");

    /* allocate the memory for FPS */
    FPS = (float *) calloc( (unsigned)(iNy*iNx), sizeof(float) );
    if (!FPS)
        mexPrintf("Memory allocation failure\n");


     //    GPU Memory Allocation
    
    cudaMalloc( (void**) &pfbx_GPU, sizeof(float)*(unsigned)((iNx+1)*iNy));
    cudaMalloc( (void**) &pfby_GPU, sizeof(float)*(unsigned)(iNx*(iNy+1)));
    cudaMalloc( (void**) &pfpenalty_GPU, sizeof(float)*(unsigned)(iNy*iNx));
    cudaMalloc( (void**) &pfdv_GPU, sizeof(float)*(unsigned)(iNy*iNx));
    cudaMalloc( (void**) &pfps_GPU, sizeof(float)*(unsigned)(iNy*iNx));
    cudaMalloc( (void**) &pfpt_GPU, sizeof(float)*(unsigned)(iNy*iNx));
    cudaMalloc( (void**) &pfgk_GPU, sizeof(float)*(unsigned)(iNy*iNx));
    cudaMalloc( (void**) &pfu_GPU, sizeof(float)*(unsigned)(iNy*iNx));
    cudaMalloc( (void**) &pfCs_GPU, sizeof(float)*(unsigned)(iNy*iNx));
    cudaMalloc( (void**) &pfCt_GPU, sizeof(float)*(unsigned)(iNy*iNx));
    cudaMalloc( (void**) &FPS_GPU, sizeof(float)*(unsigned)(iNy*iNx));

    /* Preprocessing initial values */
    for (ix=0; ix< iNx; ix++){
        idx = ix*iNy;
        for (iy=0; iy< iNy; iy++){
           index = idx + iy; 
            if (pfCs[index] < pfCt[index]){
                pfps[index] = pfCs[index];
                pfpt[index] = pfCs[index];
                pfdv[index] = pfbx[index+iNy] - pfbx[index] 
                         + pfby[index+1] - pfby[index];
            }
            else{
                pfu[index] = 1;
                pfps[index] = pfCt[index];
                pfpt[index] = pfCt[index];
                pfdv[index] = pfbx[index+iNy] - pfbx[index] 
                         + pfby[index+1] - pfby[index];
            }
        }
    }
    
    //    Copy Parameters from Host to Device

    cudaMemcpy( pfbx_GPU, pfbx, sizeof(float)*(unsigned)(iNy*(iNx+1)), cudaMemcpyHostToDevice);
    cudaMemcpy( pfby_GPU, pfby, sizeof(float)*(unsigned)((iNy+1)*iNx), cudaMemcpyHostToDevice);
    cudaMemcpy( pfpenalty_GPU, pfpenalty, sizeof(float)*(unsigned)(iNy*iNx), cudaMemcpyHostToDevice);
    cudaMemcpy( pfdv_GPU, pfdv, sizeof(float)*(unsigned)(iNy*iNx), cudaMemcpyHostToDevice);
    cudaMemcpy( pfps_GPU, pfps, sizeof(float)*(unsigned)(iNy*iNx), cudaMemcpyHostToDevice);
    cudaMemcpy( pfpt_GPU, pfpt, sizeof(float)*(unsigned)(iNy*iNx), cudaMemcpyHostToDevice);
    cudaMemcpy( pfgk_GPU, pfgk, sizeof(float)*(unsigned)(iNy*iNx), cudaMemcpyHostToDevice);
    cudaMemcpy( pfu_GPU, pfu, sizeof(float)*(unsigned)(iNy*iNx), cudaMemcpyHostToDevice);
    cudaMemcpy( pfCs_GPU, pfCs, sizeof(float)*(unsigned)(iNy*iNx), cudaMemcpyHostToDevice);
    cudaMemcpy( pfCt_GPU, pfCt, sizeof(float)*(unsigned)(iNy*iNx), cudaMemcpyHostToDevice);
    cudaMemcpy( FPS_GPU, FPS, sizeof(float)*(unsigned)(iNy*iNx), cudaMemcpyHostToDevice);

    /*  Main iterations */
    
    iNI = 0;

     /* Compute the execution configuration */

   dim3 dimBlock(BLOCK_SIZE,BLOCK_SIZE);
  
   dim3 dimGrid ( (iNy/dimBlock.x) + (!(iNy%dimBlock.x)?0:1) ,
                  (iNx/dimBlock.y) + (!(iNx%dimBlock.y)?0:1) );

   dim3 dimGrid_x ( ((iNy-1)/dimBlock.x) + (!((iNy-1)%dimBlock.x)?0:1) ,
                  (iNx/dimBlock.y) + (!(iNx%dimBlock.y)?0:1) );

   dim3 dimGrid_y ( (iNy/dimBlock.x) + (!(iNy%dimBlock.x)?0:1) ,
                  ((iNx-1)/dimBlock.y) + (!((iNx-1)%dimBlock.y)?0:1) );
    
    start_time = clock();

    printf("Start computing ......................................... \n\n");

    while( iNI<iNbIters ) 
    { 

        /* update px */
        krnl_1<<< dimGrid, dimBlock>>>(pfpt_GPU, pfps_GPU, pfu_GPU, 
                    pfgk_GPU, pfdv_GPU, cc, iNx, iNy);
       
        krnl_2<<< dimGrid_y, dimBlock>>>(pfbx_GPU, pfgk_GPU, steps, iNx, iNy);

        krnl_3<<< dimGrid_x, dimBlock>>>(pfby_GPU, pfgk_GPU, steps, iNx, iNy);
      
        /* projection step */
        krnl_4<<< dimGrid, dimBlock>>>(pfbx_GPU, pfby_GPU, pfgk_GPU, pfpenalty_GPU, iNx, iNy);
      
        krnl_5<<< dimGrid_y, dimBlock >>>(pfbx_GPU, pfgk_GPU, iNx, iNy);
    
        krnl_6<<< dimGrid_x, dimBlock >>>(pfby_GPU, pfgk_GPU, iNx, iNy);

        /* compute the divergence  */
        krnl_7<<< dimGrid, dimBlock>>>(pfbx_GPU, pfby_GPU, pfdv_GPU, iNx, iNy);

        /* update ps  */
        krnl_8<<< dimGrid, dimBlock>>>(pfps_GPU, pfpt_GPU, pfu_GPU, pfdv_GPU, pfCs_GPU, cc, iNx, iNy);
        
        /* update pt  */
        krnl_9<<< dimGrid, dimBlock>>>(pfps_GPU, pfpt_GPU, pfu_GPU, pfdv_GPU, pfCt_GPU, cc, iNx, iNy);
   
        /* update multipliers */
        krnl_10<<< dimGrid, dimBlock>>>(pfpt_GPU, pfdv_GPU, pfps_GPU, pfu_GPU, FPS_GPU, cc, iNx, iNy);
        cudaMemcpy( FPS, FPS_GPU, sizeof(float)*(unsigned)(szImg), cudaMemcpyDeviceToHost);

        fps = 0;
        for (int ii=0; ii< szImg; ii++){
                fps += FPS[ii];
        }

        pfcvg[iNI] = fps / szImg;
        
        if ((pfcvg[iNI] <= fError) && (iNI >= 2) ){
            break;
        }
        
        iNI ++;
     }   

    cudaMemcpy( pfu, pfu_GPU, sizeof(float)*(unsigned)(szImg), cudaMemcpyDeviceToHost);

    mexPrintf("Total iteration number = %i\n",iNI);
    end_time = clock();

    /* Outputs (see above) */
    punum[0] = iNI;
    
    /* Free memory */
    free( (float *) pfbx );
    free( (float *) pfby );
    free( (float *) pfps );
    free( (float *) pfpt );
    free( (float *) pfgk );
    free( (float *) pfdv );
    free( (float *) FPS );

    //    Free GPU Memory
    cudaFree(pfbx_GPU);
    cudaFree(pfby_GPU);
    cudaFree(pfpenalty_GPU);
    cudaFree(pfps_GPU);
    cudaFree(pfpt_GPU);
    cudaFree(pfgk_GPU);
    cudaFree(pfdv_GPU);
    cudaFree(pfu_GPU);
    cudaFree(pfCs_GPU);
    cudaFree(pfCt_GPU);
    cudaFree(FPS_GPU);

    tt[0] = difftime(end_time,start_time)/1000;
    mexPrintf("\nComputing Time for max-flow = %.4f sec\n \n",tt[0]);
    
    
}
/****************************************/






/**********************************************/
/************** END MAIN FUNCTION *************/
/**********************************************/
