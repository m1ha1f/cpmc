#include <stdio.h>

#define BLOCK_SIZE 8
#define BLOCK2D_SIZE 16
#define MAX(a,b) ( a > b ? a : b )
#define MIN(a,b) ( a <= b ? a : b )
#define SIGN(x) ( x >= 0.0 ? 1.0 : -1.0 )
#define ABS(x) ( (x) > 0.0 ? x : -(x) )
#define SQR(x) __fmul_rz((x), (x))

static __global__ void krnl_1(float *pfdv, float *pfpt, float *pfps, 
                        float *pfu, float *pfgk, float cc, 
                        int iNx, int iNy, int iLab, int SZF,
                        int blocksInY, float invBlocksInY)
{
   int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
   int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
   int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
   int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
   int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
   int index2d, index;

   if( idx<iNy && idy<iNx && idz<iLab)
   {
      index2d = idx +  __mul24(idy, iNy);
      index = index2d + __mul24(idz, SZF);

      pfgk[index] = __fadd_rz(__fadd_rz(__fadd_rz(pfdv[index], -pfps[index2d]), pfpt[index]), -pfu[index]);
    }
}

static __global__ void krnl_2(float *pfbx, float *pfgk, float steps, 
                int iNx, int iNy, int iLab, int SZF, int blocksInY, 
                float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    int index;

    if( idx<iNy && idy<(iNx-1) && idz<iLab )
    {
      index = idx + __mul24(idy+1, iNy) + __mul24(idz, SZF);
    
      pfbx[index] = __fadd_rz(__fmul_rz(steps, __fadd_rz(pfgk[index], -pfgk[index-iNy])),
                            pfbx[index]);
    }
}

static __global__ void krnl_3(float *pfby, float *pfgk, float steps, 
                int iNx, int iNy, int iLab, int SZF, int blocksInY, 
                float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    int index;
    
    if( idx<(iNy-1) && idy<iNx && idz<iLab)
    {
      index = idx + __mul24(idy, iNy) + __mul24(idz, SZF) + 1;
    
      pfby[index] = __fadd_rz(__fmul_rz(steps, __fadd_rz(pfgk[index], -pfgk[index-1])), 
                        pfby[index]);
    }
}

static __global__ void krnl_4(float *pfbx, float *pfby, float *pfgk, 
                    float *pfpenalty, int iNx, int iNy, int iLab, 
                    int SZF, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    float fpt;
    int index, index2d;

    if( idx<iNy && idy<iNx && idz<iLab)
    {
        index2d = idx + __mul24(idy, iNy);
        index = index2d + __mul24(idz, SZF);

        fpt = __fsqrt_rz(__fmul_rz(
            __fadd_rz(__fadd_rz(__fadd_rz(SQR(pfbx[index+iNy]), SQR(pfby[index+1])),
            SQR(pfbx[index])), SQR(pfby[index])), 0.5));
                
                if (fpt > pfpenalty[index2d])
                    pfgk[index] = __fdividef(pfpenalty[index2d], fpt);
                    
                else
                    pfgk[index] = 1;
    }
}

static __global__ void krnl_5(float *pfbx, float *pfgk, int iNx, int iNy, 
        int iLab, int SZF, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    int index;
    
    if( idx<iNy && idy<(iNx-1) && idz<iLab)
    {
      index = idx + __mul24(idy+1, iNy) + __mul24(idz, SZF);  

      pfbx[index] = __fmul_rz(__fmul_rz(__fadd_rz(pfgk[index], pfgk[index-iNy]), 0.5), pfbx[index]);
    }
}

static __global__ void krnl_6(float *pfby, float *pfgk, int iNx, 
        int iNy, int iLab, int SZF, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    int index;
    
    if( idx<(iNy-1) && idy<iNx && idz<iLab)
    {
      index = idx + __mul24(idy, iNy) + __mul24(idz, SZF) + 1;

      pfby[index] = __fmul_rz(__fmul_rz(__fadd_rz(pfgk[index], pfgk[index-1]), 0.5), pfby[index]);
    }
}


static __global__ void krnl_7(float *pfbx, float *pfby, float *pfdv,
                        float *pfps, float *pfpt, float *pfu, 
                        float *pfCt, float *pfgk, float cc,
                        int iNx, int iNy, int iLab, int SZF, 
                        int blocksInY, float invBlocksInY)
{
   int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
   int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
   int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
   int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
   int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
   int index, index2d; 
   float fpt1, fpt2;

    if( idx<iNy && idy<iNx && idz<iLab)
    {
         index2d = idx + __mul24(idy, iNy);
         index = index2d + __mul24(idz, SZF);

         fpt1 = __fadd_rz(__fadd_rz(__fadd_rz(pfbx[index+iNy], -pfbx[index]), 
                 pfby[index+1]), -pfby[index]);

         fpt2 = __fadd_rz(__fadd_rz(pfps[index2d], pfu[index]), -fpt1);

         pfdv[index] = fpt1;

         pfpt[index] = MIN(fpt2, pfCt[index]);

         pfgk[index] = fpt1 + pfpt[index];

    }
}

static __global__ void krnl_8(float *pfps, float *pfu, 
                        float *pfgk, float *FPS, float cc,
                        int iNx, int iNy, int iLab, int SZF, float ft)
{
   int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
   int idy   = __mul24(blockIdx.y,blockDim.y)+threadIdx.y;
   float fpt;
   int iz, index, index2d;

   if( idx<iNy && idy<iNx )
   {
         index = idx + __mul24(idy, iNy);

         fpt = 0;
         for (iz = 0; iz < iLab; iz++){
            index2d = index + __mul24(iz, SZF);

            fpt += __fadd_rz(pfgk[index2d], -pfu[index2d]);
         }

         pfps[index] = __fadd_rz(__fdiv_rz(fpt, iLab), ft);

         FPS[index] = 0;

         for (iz = 0; iz < iLab; iz++){
            index2d = index + __mul24(iz, SZF);

            fpt = __fadd_rz(pfgk[index2d], -pfps[index]);
            pfu[index2d] -= fpt;
            FPS[index] += ABS(fpt);
         }

    }
} 
