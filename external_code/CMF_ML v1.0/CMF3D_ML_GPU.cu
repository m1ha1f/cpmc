/***************************************************************************/
/*      Name:       CMF3D_ML_GPU.c                               
                                                
        GPU program to perform the continuous max-flow algorithm to solve the
        3D continuous-cut problem with multiple labels (Potts Model)
 
        Usage: [u, erriter, i, timet] = CMF3D_ML_GPU(penalty, C_t, para);
 
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
                    para[0,1,2]: rows, cols, heights of the given 3D data
                    para[3]: the number of labels or segmentations
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
 
               >> nvmex -f nvmexopts.bat CMF3D_ML_GPU.cu -IC:\cuda\v4.0\include ...
               >>     -LC:\cuda\v4.0\lib\x64 -lcufft -lcudart
                
                Note: The compilation process depends on your system configuration!
                      You should change the path to your own cuda installation path.
 
           Example:
 
               >> [u, erriter, i, timet] = CMF3D_ML_GPU(single(penalty), single(Ct), single(para));

               Use the matlab command imagesc or isosurface to show the computed slice or 
               3D data respectively.
 
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
#include "CMF3D_ML_kernels.cu"

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
    float   *pfbx, *pfby, *pfbz, *pfps, *pfpt, *pfdiv;
    float   fError, cc, steps;
    float   fps, fpt, *tt;
    int     *punum, iNy, iNx, iNz, iLab, iNdim, iDim[4], iNI;
    int     iNbIters, szImg, S2D, S3D, ix, iy, iz, id, iDev;
    time_t  start_time, end_time;
    
    /*    GPU Variables */
    float   *pfbx_GPU, *pfby_GPU, *pfbz_GPU, *pfpenalty_GPU;
    float   *pfdiv_GPU, *pfpt_GPU, *pfps_GPU, *pfgk_GPU, *pfu_GPU, *pfCt_GPU;
    float   *FPS, *FPS_GPU, FCP, FINV, FINVX, FINVY, FINVZ;

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
    cc = (float) pfVecParameters[6];    /* cc for ALM */
    steps = (float) pfVecParameters[7]; /* steps for each iteration */

    FCP = 1/(cc*iLab);

    printf("Initializing ................................................ \n\n");

    /* Outputs */
    /* outputs the computed u(x)  */
    iNdim = 4;
    iDim[0] = iNy;
    iDim[1] = iNx;
    iDim[2] = iNz;
    iDim[3] = iLab; 
    
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
    pfbx = (float *) calloc( (unsigned)(iNy*(iNx+1)*iNz*iLab), sizeof(float) );
    if (!pfbx)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for py1 */
    pfby = (float *) calloc( (unsigned)((iNy+1)*iNx*iNz*iLab), sizeof(float) );
    if (!pfby)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for pz1 */
    pfbz = (float *) calloc( (unsigned)(iNy*iNx*(iNz+1)*iLab), sizeof(float) );
    if (!pfbz)
        mexPrintf("Memory allocation failure\n");

    /* allocate the memory for ps */
    pfps = (float *) calloc( (unsigned)(iNy*iNx*iNz), sizeof(float) );
    if (!pfps)
        mexPrintf("Memory allocation failure\n");

    /* allocate the memory for pt */
    pfpt = (float *) calloc( (unsigned)(iNy*iNx*iNz*iLab), sizeof(float) );
    if (!pfpt)
        mexPrintf("Memory allocation failure\n");
    
    /* allocate the memory for div */
    pfdiv = (float *) calloc( (unsigned)(iNy*iNx*iNz*iLab), sizeof(float) );
    if (!pfdiv)
        mexPrintf("Memory allocation failure\n");


    /* allocate the memory for FPS */
    FPS = (float *) calloc( (unsigned)(iNy*iNx*iNz), sizeof(float) );
    if (!FPS)
        mexPrintf("Memory allocation failure\n");
    
    /* Preprocessing initial values */

    for (iz=0; iz< iNz; iz++){
        int indz = iz*S2D;
        for (ix=0; ix < iNx; ix++){
            int indy = ix*iNy;
            for (iy=0; iy < iNy; iy++){
                int index = indz + indy + iy;

                fpt = pfCt[index];
                int ik = 0;
                
                for (id = 1; id < iLab; id++){
                    int index1 = index + id*S3D;
                    if (fpt >= pfCt[index1]){
                        fpt = pfCt[index1];
                        ik = id;
                    }
                }
                    
                pfps[index] = fpt;
                pfu[index+ik*S3D] = 1/cc;
                
                for (id = 0; id < iLab; id++){        
                    pfpt[index+id*S3D] = fpt;
                }
            }
        }          
    }

/*    GPU Memory Allocation */
    
    cudaMalloc( (void**) &pfbx_GPU, sizeof(float)*(unsigned)((iNx+1)*iNy*iNz*iLab));
    cudaMalloc( (void**) &pfby_GPU, sizeof(float)*(unsigned)(iNx*(iNy+1)*iNz*iLab));
    cudaMalloc( (void**) &pfbz_GPU, sizeof(float)*(unsigned)(iNx*iNy*(iNz+1)*iLab));
    cudaMalloc( (void**) &pfpenalty_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &pfps_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));
    cudaMalloc( (void**) &pfpt_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz*iLab));
    cudaMalloc( (void**) &pfgk_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz*iLab));
    cudaMalloc( (void**) &pfu_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz*iLab));
    cudaMalloc( (void**) &pfCt_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz*iLab));
    cudaMalloc( (void**) &pfdiv_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz*iLab));
    cudaMalloc( (void**) &FPS_GPU, sizeof(float)*(unsigned)(iNy*iNx*iNz));

    /*    Copy Parameters from Host to Device */

    cudaMemcpy( pfbx_GPU, pfbx, sizeof(float)*(unsigned)(iNy*(iNx+1)*iNz*iLab), cudaMemcpyHostToDevice);
    cudaMemcpy( pfby_GPU, pfby, sizeof(float)*(unsigned)((iNy+1)*iNx*iNz*iLab), cudaMemcpyHostToDevice);
    cudaMemcpy( pfbz_GPU, pfbz, sizeof(float)*(unsigned)(iNy*iNx*(iNz+1)*iLab), cudaMemcpyHostToDevice);
    cudaMemcpy( pfpenalty_GPU, pfpenalty, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfdiv_GPU, pfdiv, sizeof(float)*(unsigned)(iNy*iNx*iNz*iLab), cudaMemcpyHostToDevice);
    cudaMemcpy( pfps_GPU, pfps, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);
    cudaMemcpy( pfpt_GPU, pfpt, sizeof(float)*(unsigned)(iNy*iNx*iNz*iLab), cudaMemcpyHostToDevice);
    cudaMemcpy( pfu_GPU, pfu, sizeof(float)*(unsigned)(iNy*iNx*iNz*iLab), cudaMemcpyHostToDevice);
    cudaMemcpy( pfCt_GPU, pfCt, sizeof(float)*(unsigned)(iNy*iNx*iNz*iLab), cudaMemcpyHostToDevice);
    cudaMemcpy( FPS_GPU, FPS, sizeof(float)*(unsigned)(iNy*iNx*iNz), cudaMemcpyHostToDevice);

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

    /* define the constants */

    FINV = 1.0f/(float)blocksInY;
    FINVY = 1.0f/(float)blocksInY_y;
    FINVX = 1.0f/(float)blocksInY_x;
    FINVZ = 1.0f/(float)blocksInY_z;

    /*  Main iterations */
    
    iNI = 0;

    printf("Start computing ......................................... \n\n");

    start_time = clock();
    
    while( iNI<iNbIters ) 
    { 

        krnl_1<<<dimGrid, dimBlock>>>(pfps_GPU, pfpt_GPU, pfdiv_GPU, pfu_GPU, pfgk_GPU, cc, iNx, 
                    iNy, iNz, iLab, S2D, S3D, blocksInY, FINV);

        krnl_2<<<dimGrid_y, dimBlock>>>(pfbx_GPU, pfgk_GPU, steps, 
              iNx, iNy, iNz, iLab, S2D, S3D, blocksInY_y, FINVY);

        krnl_3<<<dimGrid_x, dimBlock>>>(pfby_GPU, pfgk_GPU, steps, 
                          iNx, iNy, iNz, iLab, S2D, S3D, blocksInY_x, FINVX);

        krnl_z<<<dimGrid_z, dimBlock>>>(pfbz_GPU, pfgk_GPU, steps, 
                  iNx, iNy, iNz, iLab, S2D, S3D, blocksInY_z, FINVZ);

        krnl_4<<<dimGrid, dimBlock>>>(pfbx_GPU, pfby_GPU, pfbz_GPU, 
                        pfgk_GPU, pfpenalty_GPU, iNx, iNy, iNz, iLab, S2D, S3D, blocksInY, FINV);

        krnl_5<<<dimGrid_y, dimBlock>>>(pfbx_GPU, pfgk_GPU, 
                        iNx, iNy, iNz, iLab, S2D, S3D, blocksInY_y, FINVY);

        krnl_6<<< dimGrid_x, dimBlock>>>(pfby_GPU, pfgk_GPU, iNx, 
                         iNy, iNz, iLab, S2D, S3D, blocksInY_x, FINVX);
        
        krnl_zp<<<dimGrid_z, dimBlock>>>(pfbz_GPU, pfgk_GPU, 
                iNx, iNy, iNz, iLab, S2D, S3D, blocksInY_z, FINVZ);

        krnl_7<<<dimGrid, dimBlock>>>(pfbx_GPU, pfby_GPU, pfbz_GPU, 
                        pfdiv_GPU, FPS_GPU, pfps_GPU, pfpt_GPU, pfu_GPU, pfCt_GPU, cc, iNx, iNy, iNz,
                        iLab, S2D, S3D, blocksInY, FINV, FCP);
    
        cudaMemcpy( FPS, FPS_GPU, sizeof(float)*(unsigned)(S3D), cudaMemcpyDeviceToHost);
        
        fps = 0;
        for (int ii=0; ii< S3D; ii++){
                fps += FPS[ii];
        }

        pfcvg[iNI] = cc* fps / szImg;
        
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
    free( (float *) pfbz );

    free( (float *) pfps );
    free( (float *) pfpt );

    free( (float *) pfdiv );
    free( (float *) FPS );

    //    Free GPU Memory
    cudaFree(pfbx_GPU);
    cudaFree(pfby_GPU);
    cudaFree(pfbz_GPU);
    cudaFree(pfpenalty_GPU);
    cudaFree(pfdiv_GPU);
    cudaFree(pfps_GPU);
    cudaFree(pfpt_GPU);
    cudaFree(pfgk_GPU);
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
