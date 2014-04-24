/***************************************************************************/
/*      Name:       CMF3D_GPU.c                               
                                                
        Performing the continuous max-flow algorithm to solve the
        3D continuous min-cut problem over GPU (Nvidia CUDA based)
 
        Usage: [u, erriter, i, timet] = CMF3D_GPU(penalty, C_s, C_t, para);
 
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
                    para[0,1,2]: rows, cols, heights of the given image
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
 
               >> nvmex -f nvmexopts.bat CMF3D_GPU.cu -IC:\cuda\v4.0\include -LC:\cuda\v4.0\lib\x64 -lcufft -lcudart
                
                Note: The compilation process depends on your system configuration!
                      You should change the path to your own cuda installation path.
 
           Example:
 
               >> [u, erriter, i, timet] = CMF3D_GPU(single(penalty), single(Cs), single(Ct), single(para));

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

#include "cuda.h"
#include "device_functions.h"
#include "CMF3D_kernels.cu"

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
    float   *pfbx, *pfby, *pfbz, *pfps, *pfpt, *pfgk, *tt, *pfdv;
    float   fError, cc, steps, fps;
    int     *punum, iNy, iNx, iNz, iNdim, iDim[3], ix, iy, iz, iNI;
    int     iNbIters, szImg, SZF, idx, index, idz, iDev;
    time_t  start_time, end_time;

    /*    GPU Variables */
    float   *pfbx_GPU, *pfby_GPU, *pfbz_GPU, *pfpenalty_GPU, *pfdv_GPU;
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
    pfpenalty = (float *)mxGetData(pmxIn[0]); /* Given penalty */
    pfCs = (float *)mxGetData(pmxIn[1]); /* bound of source flows */
    pfCt = (float *)mxGetData(pmxIn[2]); /* bound of sink flows */
    pfVecParameters = (float *)mxGetData(pmxIn[3]); /* Vector of parameters */
    
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
    szImg = iNy*iNx*iNz;
    SZF = iNy*iNx;
    
    /* Choice of region segmentation model */
    iNbIters = (int) pfVecParameters[3]; /* total number of iterations */
    fError = (float) pfVecParameters[4]; /* error criterion */
    cc = (float) pfVecParameters[5]; /* cc for ALM */
    steps = (float) pfVecParameters[6]; /* steps for each iteration */
 
   printf("Initializing ................................................ \n\n");

    /* Outputs */
    /* outputs the computed u(x)  */
    iNdim = 3;
    iDim[0] = iNy;
    iDim[1] = iNx;
    iDim[2] = iNz;
    
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
    tt = (float*)mxGetData(pmxOut[3]);
    
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

    /* allocate the memory for FPS */
    FPS = (float *) calloc( (unsigned)(iNy*iNx*iNz), sizeof(float) );
    if (!FPS)
        mexPrintf("Memory allocation failure\n");

    /* Preprocessing initial values */
    for (iz=0; iz < iNz; iz++){
        idz = iz*SZF;
        for (ix=0; ix < iNx; ix++){
            idx = idz + ix*iNy;
            for (iy=0; iy < iNy; iy++){
                index = idx + iy;
                
                if (pfCs[index] < pfCt[index]){
                    pfps[index] = pfCs[index];
                    pfpt[index] = pfCs[index];
                    pfdv[index] = pfbx[index+iNy] - pfbx[index] 
                                + pfby[index+1] - pfby[index]
                                + pfbz[index+SZF] - pfbz[index];
                }
                else{
                    pfu[index] = 1/cc;
                    pfps[index] = pfCt[index];
                    pfpt[index] = pfCt[index];
                    pfdv[index] = pfbx[index+iNy] - pfbx[index] 
                                + pfby[index+1] - pfby[index]
                                + pfbz[index+SZF] - pfbz[index];
                }
            }
        }
    }
 
    /*    GPU Memory Allocation */
    
    cudaMalloc( (void**) &pfbx_GPU, sizeof(float)*(unsigned)((iNx+1)*iNy*iNz));
    cudaMalloc( (void**) &pfby_GPU, sizeof(float)*(unsigned)(iNx*(iNy+1)*iNz));
    cudaMalloc( (void**) &pfbz_GPU, sizeof(float)*(unsigned)(iNx*iNy*(iNz+1)));
    cudaMalloc( (void**) &pfps_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &pfpt_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &pfgk_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &pfdv_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &pfpenalty_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &pfu_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &pfCs_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &pfCt_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &FPS_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));

    /*    Copy Parameters from Host to Device */

    cudaMemcpy( pfbx_GPU, pfbx, sizeof(float)*(unsigned)(iNy*(iNx+1)*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfby_GPU, pfby, sizeof(float)*(unsigned)((iNy+1)*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfbz_GPU, pfbz, sizeof(float)*(unsigned)((iNz+1)*iNx*iNy), cudaMemcpyHostToDevice);
    cudaMemcpy( pfps_GPU, pfps, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfpt_GPU, pfpt, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfdv_GPU, pfdv, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfpenalty_GPU, pfpenalty, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfu_GPU, pfu, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfCs_GPU, pfCs, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfCt_GPU, pfCt, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);

     /*  Main iterations */
    
    iNI = 0;
     
    /* Compute the execution configuration */

   dim3 dimBlock(BLOCK_SIZE,BLOCK_SIZE,BLOCK_SIZE);
  
   int blocksInX = (iNy/dimBlock.x) + (!(iNy%dimBlock.x)?0:1);
   int blocksInY = (iNx/dimBlock.y) + (!(iNx%dimBlock.y)?0:1);
   int blocksInZ = (iNz/dimBlock.z) + (!(iNz%dimBlock.z)?0:1);

   dim3 dimGrid ( blocksInX, blocksInY*blocksInZ);

   blocksInX = ((iNy-1)/dimBlock.x) + (!((iNy-1)%dimBlock.x)?0:1);
   int blocksInY_x = (iNx/dimBlock.y) + (!(iNx%dimBlock.y)?0:1);
   blocksInZ = (iNz/dimBlock.z) + (!(iNz%dimBlock.z)?0:1);

   dim3 dimGrid_x (blocksInX, blocksInY_x*blocksInZ);

   blocksInX = (iNy/dimBlock.x) + (!(iNy%dimBlock.x)?0:1);
   int blocksInY_y = ((iNx-1)/dimBlock.y) + (!((iNx-1)%dimBlock.y)?0:1);
   blocksInZ = (iNz/dimBlock.z) + (!(iNz%dimBlock.z)?0:1);

   dim3 dimGrid_y ( blocksInX, blocksInY_y*blocksInZ);

   blocksInX = (iNy/dimBlock.x) + (!(iNy%dimBlock.x)?0:1);
   int blocksInY_z = (iNx/dimBlock.y) + (!(iNx%dimBlock.y)?0:1);
   blocksInZ = ((iNz-1)/dimBlock.z) + (!((iNz-1)%dimBlock.z)?0:1);
 
   dim3 dimGrid_z (blocksInX, blocksInY_z*blocksInZ);
    
   start_time = clock();

    printf("Start computing ......................................... \n\n");

    while( iNI<iNbIters ) 
    { 
        /* update p */

        krnl_1<<<dimGrid, dimBlock>>>(pfdv_GPU, pfpt_GPU, pfps_GPU, pfu_GPU, pfgk_GPU, cc, iNx, 
                        iNy, iNz, SZF, blocksInY, 1.0f/(float)blocksInY);
        
        /* update px py and pz */

        krnl_2<<<dimGrid_y, dimBlock>>>(pfbx_GPU, pfgk_GPU, steps, 
                  iNx, iNy, iNz, SZF, blocksInY_y, 1.0f/(float)blocksInY_y);
    
        krnl_3<<<dimGrid_x, dimBlock>>>(pfby_GPU, pfgk_GPU, steps, 
                  iNx, iNy, iNz, SZF, blocksInY_x, 1.0f/(float)blocksInY_x);

        krnl_z<<<dimGrid_z, dimBlock>>>(pfbz_GPU, pfgk_GPU, steps, 
                  iNx, iNy, iNz, SZF, blocksInY_z, 1.0f/(float)blocksInY_z);

        /* projection step */

        krnl_4<<<dimGrid, dimBlock>>>(pfbx_GPU, pfby_GPU, pfbz_GPU, 
                pfgk_GPU, pfpenalty_GPU, iNx, iNy, iNz, SZF, blocksInY, 1.0f/(float)blocksInY);

        krnl_5<<<dimGrid_y, dimBlock>>>(pfbx_GPU, pfgk_GPU, 
                iNx, iNy, iNz, SZF, blocksInY_y, 1.0f/(float)blocksInY_y);

        krnl_6<<< dimGrid_x, dimBlock>>>(pfby_GPU, pfgk_GPU, iNx, 
                iNy, iNz, SZF, blocksInY_x, 1.0f/(float)blocksInY_x);

        krnl_zp<<<dimGrid_z, dimBlock>>>(pfbz_GPU, pfgk_GPU, iNx, 
                iNy, iNz, SZF, blocksInY_z, 1.0f/(float)blocksInY_z);
      
        krnl_7<<<dimGrid, dimBlock>>>(pfbx_GPU, pfby_GPU, pfbz_GPU,
                        pfdv_GPU, pfps_GPU, pfpt_GPU, pfu_GPU, FPS_GPU,
                        pfCs_GPU, pfCt_GPU, 1.0f/(float)cc, iNx, iNy, iNz, SZF, 
                        blocksInY, 1.0f/(float)blocksInY);
    
//        krnl_dot<<<dimGrid1D, dimBlock1D>>>(FPS_GPU, szImg, sumFPS_GPU);

        cudaMemcpy( FPS, FPS_GPU, sizeof(float)*unsigned(szImg), cudaMemcpyDeviceToHost);

        fps = 0;
        for (int ii=0; ii< szImg; ii++){
                fps += FPS[ii];
        }

        pfcvg[iNI] = cc * fps / szImg;
        
        if (pfcvg[iNI] <= fError)
            break;
                
        iNI ++;
 
    }   

    cudaMemcpy( pfu, pfu_GPU, sizeof(float)*(unsigned)(szImg), cudaMemcpyDeviceToHost);

    mexPrintf("number of iterations = %i\n",iNI);
    end_time = clock();
    
    for (iz=0; iz < iNz; iz++){
        idz = iz*SZF;
        for (ix=0; ix < iNx; ix++){
            idx = idz + ix*iNy;
            for (iy=0; iy < iNy; iy++){
                index = idx + iy;
                
                pfu[index] = cc*pfu[index];
            }
        }
    }

    /* Outputs (see above) */
    punum[0] = iNI;

    
    /* Free memory */
    free( (float *) pfbx );
    free( (float *) pfby );
    free( (float *) pfbz );
    free( (float *) pfps );
    free( (float *) pfpt );
    free( (float *) pfdv );
    free( (float *) pfgk );

    free( (float *) FPS );

    //    Free GPU Memory
    cudaFree(pfbx_GPU);
    cudaFree(pfby_GPU);
    cudaFree(pfbz_GPU);
    cudaFree(pfps_GPU);
    cudaFree(pfpt_GPU);
    cudaFree(pfgk_GPU);
    cudaFree(pfdv_GPU);
    cudaFree(pfpenalty_GPU);
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
