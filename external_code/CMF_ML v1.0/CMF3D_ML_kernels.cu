#include <stdio.h>

#define BLOCK_SIZE 8
#define BLOCK2D_SIZE 16
#define MAX(a,b) ( a > b ? a : b )
#define MIN(a,b) ( a <= b ? a : b )
#define SIGN(x) ( x >= 0.0 ? 1.0 : -1.0 )
#define ABS(x) ( (x) > 0.0 ? x : -(x) )
#define SQR(x) __fmul_rz((x), (x))

static __global__ void krnl_1( float *pfps, float *pfpt, float *pfdiv, 
                        float *pfu, float *pfgk, float cc, 
                        int iNx, int iNy, int iNz, int iLab, int S2D,
                        int S3D, int blocksInY, float invBlocksInY)
{
   int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
   int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
   int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
   int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
   int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
   int index3d, index, id;

   if( idx<iNy && idy<iNx && idz<iNz)
   {
       index3d = idx + __umul24(idy,iNy) + __umul24(idz,S2D);
      
       for(id = 0; id < iLab; id ++){
            index = index3d + __umul24(id,S3D);

            pfgk[index] = __fadd_rz(__fadd_rz(__fadd_rz(pfdiv[index], - pfps[index3d]), 
                    pfpt[index]), -pfu[index]);
        }
    }
}

static __global__ void krnl_2(float *pfbx, float *pfgk, float steps, 
                int iNx, int iNy, int iNz, int iLab, int S2D, int S3D, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    int index3d, index, id;

    if( idx<iNy && idy<(iNx-1) && idz<iNz)
    {
        index3d = idx + __umul24(idy+1,iNy) + __umul24(idz,S2D);

        for(id = 0; id < iLab; id ++){
            index = index3d + __umul24(id,S3D);

            pfbx[index] = __fadd_rz(__fmul_rz(steps, __fadd_rz(pfgk[index], -pfgk[index-iNy])), 
                            pfbx[index]);
        }
    }
}

static __global__ void krnl_3(float *pfby, float *pfgk, float steps, 
                int iNx, int iNy, int iNz, int iLab, int S2D, int S3D, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    int index3d, index, id;
    
    if( idx<(iNy-1) && idy<iNx && idz<iNz)
    {
        index3d = (idx+1) + __umul24(idy,iNy) + __umul24(idz,S2D);

        for(id = 0; id < iLab; id ++){
            index = index3d + __umul24(id, S3D);

            pfby[index] = __fadd_rz(__fmul_rz(steps, __fadd_rz(pfgk[index], -pfgk[index-1])), 
                            pfby[index]);
        }
    }
}

static __global__ void krnl_z(float *pfbz, float *pfgk, float steps, 
                int iNx, int iNy, int iNz, int iLab, int S2D, int S3D, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    int index3d, index, id;
    
    if( idx<iNy && idy<iNx && idz<(iNz-1))
    {
        index3d = idx + __umul24(idy,iNy) + __umul24(idz+1,S2D);

        for(id = 0; id < iLab; id ++){
            index = index3d + __umul24(id,S3D);

            pfbz[index] = __fadd_rz(__fmul_rz(steps, __fadd_rz(pfgk[index], -pfgk[index-S2D])), 
                            pfbz[index]);
        }
    }
}

static __global__ void krnl_4(float *pfbx, float *pfby, float *pfbz, float *pfgk, float *pfpenalty, 
                    int iNx, int iNy, int iNz, int iLab, int S2D, int S3D, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    int index3d, index, id;
    float fpt;

    if( idx<iNy && idy<iNx && idz<iNz)
    {
        index3d = idx + __umul24(idy,iNy) + __umul24(idz,S2D);

        for(id = 0; id < iLab; id++){
            index = index3d + __umul24(id,S3D);

            fpt = __fsqrt_rz(__fmul_rz(SQR(pfbx[index+iNy]) + SQR(pfbx[index]) 
                      + SQR(pfby[index]) + SQR(pfby[index+1]) 
                      + SQR(pfbz[index]) + SQR(pfbz[index+S2D]),0.5));
                
            if (fpt > pfpenalty[index3d])
                 pfgk[index] = __fdividef(pfpenalty[index3d],fpt);
            else
                 pfgk[index] = 1;
                
        }
    }
}

static __global__ void krnl_5(float *pfbx, float *pfgk, int iNx, int iNy, int iNz,
        int iLab, int S2D, int S3D, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    int index3d, index, id;
    
    if( idx<iNy && idy<(iNx-1) && idz<iNz)
    {
        index3d = idx + __umul24(idy+1, iNy) + __umul24(idz,S2D);  

        for(id = 0; id < iLab; id ++){
            index = index3d + __umul24(id,S3D);

            pfbx[index] = __fmul_rz(__fmul_rz(__fadd_rz(pfgk[index], pfgk[index-iNy]), 
                            0.5), pfbx[index]);
        }
    }
}

static __global__ void krnl_6(float *pfby, float *pfgk, int iNx, 
        int iNy, int iNz, int iLab, int S2D, int S3D, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    int index3d, index, id;
    
    if( idx<(iNy-1) && idy<iNx && idz<iNz)
    {
        index3d = (idx+1) + __umul24(idy,iNy) + __umul24(idz,S2D);

        for(id = 0; id < iLab; id ++){
            index = index3d + __umul24(id,S3D);

            pfby[index] = __fmul_rz(__fmul_rz(__fadd_rz(pfgk[index-1], pfgk[index]),
                            0.5), pfby[index]);
        }
    }
}

static __global__ void krnl_zp(float *pfbz, float *pfgk, int iNx, 
        int iNy, int iNz, int iLab, int S2D, int S3D, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    int index3d, index, id;
    
    if( idx<iNy && idy<iNx && idz<(iNz-1))
    {
        index3d = idx + __umul24(idy,iNy) + __umul24(idz+1,S2D);

        for(id = 0; id < iLab; id ++){
            index = index3d + __umul24(id,S3D);

            pfbz[index] = __fmul_rz(__fmul_rz(__fadd_rz(pfgk[index], pfgk[index-S2D]), 
                            0.5), pfbz[index]);
        }
    }
}

static __global__ void krnl_7(float *pfbx, float *pfby, float *pfbz,
                        float *pfdiv, float *FPS, float *pfps, float *pfpt, float *pfu, 
                        float *pfCt, float cc, int iNx, int iNy, int iNz,
                        int iLab, int S2D, int S3D, int blocksInY, 
                        float invBlocksInY, float FCP)
{
   int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
   int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
   int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
   int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
   int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
   int index3d, index, id;
   float fpt, fps;

    if( idx<iNy && idy<iNx && idz<iNz)
    {
        index3d = idx + __umul24(idy,iNy) + __umul24(idz,S2D);
        fps = 0;

        for(id = 0; id < iLab; id ++){

            index = index3d + __umul24(id,S3D);

            /* recompute the divergence */

            pfdiv[index] = pfbx[index+iNy] - pfbx[index] 
                         + pfby[index+1] - pfby[index]
                         + pfbz[index+S2D] - pfbz[index];

            fpt = pfps[index3d] + pfu[index] - pfdiv[index];

            pfpt[index] = MIN(fpt , pfCt[index]);

            fps += __fadd_rz(__fadd_rz(pfdiv[index], pfpt[index]), -pfu[index]); 
        }
        
        pfps[index3d] = __fadd_rz(__fdiv_rz(fps,iLab), FCP);

        FPS[index3d] = 0;

         /* update the multipliers */

        for (id = 0; id < iLab; id++){
            index = index3d + __umul24(id,S3D);

            fpt = __fadd_rz(__fadd_rz(pfdiv[index], pfpt[index]), -pfps[index3d]);

            pfu[index] -= fpt;

            FPS[index3d] += ABS(fpt);
         }
    }
}
