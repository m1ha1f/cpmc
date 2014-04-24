/***************************************************************************/
/*      Name:       CMF_ML_GPU.c                               
                                                
        GPU program to perform the continuous max-flow algorithm to solve the
        2D continuous-cut problem with multiple labels (Potts Model)
 
        Usage: [u, erriter, i, timet] = CMF_ML_GPU(penalty, C_t, para);
 
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
 
               >> nvmex -f nvmexopts.bat CMF_ML_GPU.cu -IC:\cuda\v4.0\include -LC:\cuda\v4.0\lib\x64 -lcufft -lcudart
                
                Note: The compilation process depends on your system configuration!
                      You should change the path to your own cuda installation path.
 
           Example:
 
               >> [u, erriter, i, timet] = CMF_ML_GPU(single(penalty), single(Ct), single(para));

               >> To demonstrate the final labeled result, see the related matlab file.
 
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

/*  
    Compilation commands (under matlab) depends on the system: 
    Windows, Linux, MacOS.

    For details, please check the website: 
        http://developer.nvidia.com/matlab-cuda
*/

#include <stdio.h>
#include <stdlib.h>
#include <mex.h>
#include <math.h>
#include <time.h>

#include "cuda.h"
#include "device_functions.h"
#include "CMF_ML_kernels.cu"

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
    
    
    float   *pfu, *pfCt, *pfcvg, *pfVecParameters, *pfpenalty;
    float   *pfbx, *pfby, *pfps, *pfpt, *pfdv;
    float   fError, cc, steps;
    float   fpt, fps, *tt;
    int     *punum, iNy, iNx, iLab, iNdim, iDim[4], ix, iy, iz, id, iNI;
    int     iNbIters, szImg, SZF, idx, index, idz, iDev;
    time_t  start_time, end_time;
    
    /*    GPU Variables */
    float   *pfbx_GPU, *pfby_GPU, *pfpenalty_GPU, *pfdv_GPU;
    float   *pfps_GPU, *pfpt_GPU, *pfgk_GPU, *pfu_GPU, *pfCt_GPU;
    float   *FPS, *FPS_GPU, FCP, FINV, FINVX, FINVY;

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
    pfpenalty = (float*)mxGetData(pmxIn[0]);        /* penalty parameters */
    pfCt = (float*)mxGetData(pmxIn[1]);             /* bound of sink flows */
    pfVecParameters = (float*)mxGetData(pmxIn[2]);  /* Vector of parameters */
    
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
    SZF = iNy*iNx;

    iNbIters = (int) pfVecParameters[3]; /* total number of iterations */
    fError = (float) pfVecParameters[4]; /* error criterion */
    cc = (float) pfVecParameters[5]; /* cc for ALM */
    steps = (float) pfVecParameters[6]; /* steps for each iteration */
   
    FCP = 1/(cc*iLab);

    printf("Initializing ................................................ \n\n");

    /* Outputs */
    /* outputs the computed u(x)  */
    iNdim = 3;
    iDim[0] = iNy;
    iDim[1] = iNx;
    iDim[2] = iLab; 
    
    pmxOut[0] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    pfu = (float*)mxGetData(pmxOut[0]);
   
    /* outputs the convergence rate  */
    iNdim = 2;
    iDim[0] = 1;
    iDim[1] = iNbIters;
    pmxOut[1] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    pfcvg = (float*)mxGetData(pmxOut[1]);
    
    /* outputs the iteration number  */
    iNdim = 2;
    iDim[0] = 1;
    iDim[1] = 1;
    pmxOut[2] = mxCreateNumericArray(iNdim,(const int*)iDim,mxUINT16_CLASS,mxREAL);
    punum = (int*)mxGetData(pmxOut[2]);

    /* outputs the computation time  */
    iNdim = 2;
    iDim[0] = 1;
    iDim[1] = 1;
    pmxOut[3] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
    tt = (float *)mxGetData(pmxOut[3]);
    
    /* Memory allocation */
    /* allocate the memory for px1 */
    pfbx = (float *) calloc( (unsigned)(iNy*(iNx+1)*iLab), sizeof(float) );
    if (!pfbx)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for py1 */
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
    
    /* allocate the memory for divergence */
    pfdv = (float *) calloc( (unsigned)(iNy*iNx*iLab), sizeof(float) );
    if (!pfdv)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for FPS */
    FPS = (float *) calloc( (unsigned)(iNy*iNx), sizeof(float) );
    if (!FPS)
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
                idz = index + iz*SZF;

                pfdv[idz] = pfby[idz+1] - pfby[idz] 
                          + pfbx[idz+iNy] - pfbx[idz];
                
                if (fpt >= pfCt[idz]){
                    fpt = pfCt[idz];
                    id = iz;
                }
            }

            pfps[index] = fpt;
            pfu[index + id*SZF] = 1/cc;

            for (iz = 0; iz < iLab; iz++){
                pfpt[index + iz*SZF] = fpt;
            }
        }
    }
    
/*    GPU Memory Allocation */
    
    cudaMalloc( (void**) &pfbx_GPU, sizeof(float)*(unsigned)((iNx+1)*iNy*iLab));
    cudaMalloc( (void**) &pfby_GPU, sizeof(float)*(unsigned)(iNx*(iNy+1)*iLab));
    cudaMalloc( (void**) &pfps_GPU, sizeof(float)*(unsigned)(iNy*iNx));
    cudaMalloc( (void**) &pfpenalty_GPU, sizeof(float)*(unsigned)(iNy*iNx));
    cudaMalloc( (void**) &pfpt_GPU, sizeof(float)*(unsigned)(iNy*iNx*iLab));
    cudaMalloc( (void**) &pfgk_GPU, sizeof(float)*(unsigned)(iNy*iNx*iLab));
    cudaMalloc( (void**) &pfdv_GPU, sizeof(float)*(unsigned)(iNy*iNx*iLab));
    cudaMalloc( (void**) &pfu_GPU, sizeof(float)*(unsigned)(iNy*iNx*iLab));
    cudaMalloc( (void**) &pfCt_GPU, sizeof(float)*(unsigned)(iNy*iNx*iLab));
    cudaMalloc( (void**) &FPS_GPU, sizeof(float)*(unsigned)(iNy*iNx));

    /*    Copy Parameters from Host to Device */

    cudaMemcpy( pfbx_GPU, pfbx, sizeof(float)*(unsigned)(iNy*(iNx+1)*iLab), cudaMemcpyHostToDevice);
    cudaMemcpy( pfby_GPU, pfby, sizeof(float)*(unsigned)((iNy+1)*iNx*iLab), cudaMemcpyHostToDevice);
    cudaMemcpy( pfps_GPU, pfps, sizeof(float)*(unsigned)(iNy*iNx), cudaMemcpyHostToDevice);
    cudaMemcpy( pfpenalty_GPU, pfpenalty, sizeof(float)*(unsigned)(iNy*iNx), cudaMemcpyHostToDevice);
    cudaMemcpy( pfpt_GPU, pfpt, sizeof(float)*(unsigned)(iNy*iNx*iLab), cudaMemcpyHostToDevice);
    cudaMemcpy( pfdv_GPU, pfdv, sizeof(float)*(unsigned)(iNy*iNx*iLab), cudaMemcpyHostToDevice);
    cudaMemcpy( pfu_GPU, pfu, sizeof(float)*(unsigned)(iNy*iNx*iLab), cudaMemcpyHostToDevice);
    cudaMemcpy( pfCt_GPU, pfCt, sizeof(float)*(unsigned)(iNy*iNx*iLab), cudaMemcpyHostToDevice);
    cudaMemcpy( FPS_GPU, FPS, sizeof(float)*(unsigned)(iNy*iNx), cudaMemcpyHostToDevice);
    
    /* Compute the execution configuration */

   dim3 dimBlock(BLOCK_SIZE,BLOCK_SIZE,BLOCK_SIZE);
  
   int blocksInX = (iNy/dimBlock.x) + (!(iNy%dimBlock.x)?0:1);
   int blocksInY = (iNx/dimBlock.y) + (!(iNx%dimBlock.y)?0:1);
   int blocksInZ = (iLab/dimBlock.z) + (!(iLab%dimBlock.z)?0:1);

   dim3 dimGrid ( blocksInX, blocksInY*blocksInZ);

   blocksInX = ((iNy-1)/dimBlock.x) + (!((iNy-1)%dimBlock.x)?0:1);
   int blocksInY_x = (iNx/dimBlock.y) + (!(iNx%dimBlock.y)?0:1);
   blocksInZ = (iLab/dimBlock.z) + (!(iLab%dimBlock.z)?0:1);

   dim3 dimGrid_x (blocksInX, blocksInY_x*blocksInZ);

   blocksInX = (iNy/dimBlock.x) + (!(iNy%dimBlock.x)?0:1);
   int blocksInY_y = ((iNx-1)/dimBlock.y) + (!((iNx-1)%dimBlock.y)?0:1);
   blocksInZ = (iLab/dimBlock.z) + (!(iLab%dimBlock.z)?0:1);

   dim3 dimGrid_y ( blocksInX, blocksInY_y*blocksInZ);

   dim3 dimBlock2d(BLOCK2D_SIZE,BLOCK2D_SIZE);
  
   dim3 dimGrid2d ( (iNy/dimBlock2d.x) + (!(iNy%dimBlock2d.x)?0:1) ,
                  (iNx/dimBlock2d.y) + (!(iNx%dimBlock2d.y)?0:1) );

    /* define the constants */

    FINV = 1.0f/(float)blocksInY;
    FINVY = 1.0f/(float)blocksInY_y;
    FINVX = 1.0f/(float)blocksInY_x;

    /*  Main iterations */
    
    iNI = 0;
    start_time = clock();

    printf("Start computing ......................................... \n\n");

    while( iNI<iNbIters ) 
    { 
        
        /* update p */
        krnl_1<<<dimGrid, dimBlock>>>(pfdv_GPU,
                        pfpt_GPU, pfps_GPU, pfu_GPU, pfgk_GPU, cc, iNx, 
                        iNy, iLab, SZF, blocksInY, FINV);

        krnl_2<<<dimGrid_y, dimBlock>>>(pfbx_GPU, pfgk_GPU, steps, 
                  iNx, iNy, iLab, SZF, blocksInY_y, FINVY);
       
        krnl_3<<<dimGrid_x, dimBlock>>>(pfby_GPU, pfgk_GPU, steps, 
                  iNx, iNy, iLab, SZF, blocksInY_x, FINVX);

        krnl_4<<<dimGrid, dimBlock>>>(pfbx_GPU, pfby_GPU,  
                pfgk_GPU, pfpenalty_GPU, iNx, iNy, iLab, SZF, blocksInY, FINV);
 
        krnl_5<<<dimGrid_y, dimBlock>>>(pfbx_GPU, pfgk_GPU, 
                iNx, iNy, iLab, SZF, blocksInY_y, FINVY);

        krnl_6<<< dimGrid_x, dimBlock>>>(pfby_GPU, pfgk_GPU, iNx, 
                iNy, iLab, SZF, blocksInY_x, FINVX);

        krnl_7<<<dimGrid, dimBlock>>>(pfbx_GPU, pfby_GPU, pfdv_GPU,
                        pfps_GPU, pfpt_GPU, pfu_GPU,
                        pfCt_GPU, pfgk_GPU, cc, iNx, iNy, iLab, SZF, 
                        blocksInY, FINV);

        krnl_8<<< dimGrid2d, dimBlock2d>>>(pfps_GPU, pfu_GPU, pfgk_GPU, FPS_GPU, cc,
                        iNx, iNy, iLab, SZF, FCP);    
        
        cudaMemcpy( FPS, FPS_GPU, sizeof(float)*(unsigned)(SZF), cudaMemcpyDeviceToHost);
        
        fps = 0;
        for (int ii=0; ii< SZF; ii++){
                fps += FPS[ii];
        }

        pfcvg[iNI] = cc * fps / szImg;
        
        if (pfcvg[iNI] <= fError)
            break;
                
        iNI ++;
     }   
    
    cudaMemcpy( pfu, pfu_GPU, sizeof(float)*(unsigned)(szImg), cudaMemcpyDeviceToHost);

    end_time = clock();

    /* Outputs (see above) */
    punum[0] = ++iNI;
    mexPrintf("number of iterations = %i\n",iNI);
    
    /* Free memory */
    free( (float *) pfbx );
    free( (float *) pfby );
    free( (float *) pfps );
    free( (float *) pfpt );
    free( (float *) pfdv );
    free( (float *) FPS );

    //    Free GPU Memory
    cudaFree(pfbx_GPU);
    cudaFree(pfby_GPU);
    cudaFree(pfps_GPU);
    cudaFree(pfpenalty_GPU);
    cudaFree(pfpt_GPU);
    cudaFree(pfgk_GPU);
    cudaFree(pfdv_GPU);
    cudaFree(pfu_GPU);
    cudaFree(pfCt_GPU);
    cudaFree(FPS_GPU);
    
    tt[0] = difftime(end_time, start_time)/1000000;
   
    mexPrintf("\nComputing Time for max-flow = %.4f sec\n \n", tt[0]);
    
}
/****************************************/






/**********************************************/
/************** END MAIN FUNCTION *************/
/**********************************************/
