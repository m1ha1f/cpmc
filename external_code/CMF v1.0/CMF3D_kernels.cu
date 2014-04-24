#define BLOCK_SIZE 8
#define MAX(a,b) ( a > b ? a : b )
#define MIN(a,b) ( a <= b ? a : b )
#define SIGN(x) ( x >= 0.0 ? 1.0 : -1.0 )
#define ABS(x) ( (x) > 0.0 ? x : -(x) )
#define SQR(x) __fmul_rz((x), (x))

static __global__ void krnl_1(float *pfdv, float *pfpt, float *pfps, 
                        float *pfu, float *pfgk, float cc, 
                        int iNx, int iNy, int iNz, int SZF,
                        int blocksInY, float invBlocksInY)
{
   int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
   int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
   int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
   int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
   int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;

   if( idx<iNy && idy<iNx && idz<iNz)
   {
      int index = idx + __mul24(idy, iNy) + __mul24(idz, SZF);

      pfgk[index] = __fadd_rz(__fadd_rz(__fadd_rz(pfdv[index], - pfps[index]), pfpt[index]), - pfu[index]);
    }
}

static __global__ void krnl_2(float *pfbx, float *pfgk, float steps, 
                int iNx, int iNy, int iNz, int SZF, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;

    if( idx<iNy && idy<(iNx-1) && idz<iNz )
    {
      int index = idx + __mul24(idy+1, iNy) + __mul24(idz, SZF);
    
      pfbx[index] = __fadd_rz(__fmul_rz(steps, __fadd_rz(pfgk[index], - pfgk[index-iNy])), pfbx[index]);
    }
}

static __global__ void krnl_3(float *pfby, float *pfgk, float steps, 
                int iNx, int iNy, int iNz, int SZF, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    
    if( idx<(iNy-1) && idy<iNx && idz<iNz)
    {
      int index =idx + __mul24(idy, iNy) + __mul24(idz, SZF) + 1;
    
      pfby[index] = __fadd_rz(__fmul_rz(steps, __fadd_rz(pfgk[index], - pfgk[index-1])), pfby[index]);
    }
}

static __global__ void krnl_z(float *pfbz, float *pfgk, float steps, 
                int iNx, int iNy, int iNz, int SZF, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    
    if( idx<iNy && idy<iNx && idz<(iNz-1))
    {
      int index = idx + __mul24(idy, iNy) + __mul24(idz+1, SZF);
    
      pfbz[index] = __fadd_rz(__fmul_rz(steps, __fadd_rz(pfgk[index], - pfgk[index-SZF])), pfbz[index]);
    }
}

static __global__ void krnl_4(float *pfbx, float *pfby, float *pfbz, float *pfgk, float *pfpenalty, 
                    int iNx, int iNy, int iNz, int SZF, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    float fpt;

    if( idx<iNy && idy<iNx && idz<iNz)
    {
      int index = idx + __mul24(idy, iNy) + __mul24(idz, SZF);

      fpt = __fadd_rz(__fadd_rz(__fadd_rz(__fadd_rz(__fadd_rz(
                SQR(pfbx[index+iNy]), SQR(pfbx[index])), SQR(pfby[index])),
                SQR(pfby[index+1])), SQR(pfbz[index])), SQR(pfbz[index+SZF])); 

      fpt = sqrtf(__fmul_rz(fpt, 0.5));
                
                if (fpt > pfpenalty[index])
                    fpt = __fdividef(fpt, pfpenalty[index]);
                else
                    fpt = 1;
                
                pfgk[index] = __frcp_rz(fpt);
    }
}

static __global__ void krnl_5(float *pfbx, float *pfgk, int iNx, int iNy, 
        int iNz, int SZF, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    
    if( idx<iNy && idy<(iNx-1) && idz<iNz)
    {
      int index = idx + __mul24(idy+1, iNy) + __mul24(idz, SZF);  

      pfbx[index] = __fmul_rz(__fmul_rz(__fadd_rz(pfgk[index], pfgk[index-iNy]), 0.5), pfbx[index]);
    }
}

static __global__ void krnl_6(float *pfby, float *pfgk, int iNx, 
        int iNy, int iNz, int SZF, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    
    if( idx<(iNy-1) && idy<iNx && idz<iNz)
    {
      int index = idx + __mul24(idy, iNy) + __mul24(idz, SZF)+1;

      pfby[index] = __fmul_rz(__fmul_rz(__fadd_rz(pfgk[index], pfgk[index-1]), 0.5), pfby[index]);
    }
}

static __global__ void krnl_zp(float *pfbz, float *pfgk, int iNx, 
        int iNy, int iNz, int SZF, int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    
    if( idx<iNy && idy<iNx && idz<(iNz-1))
    {
      int index = idx + __mul24(idy, iNy) + __mul24(idz+1, SZF);

      pfbz[index] = __fmul_rz(__fmul_rz(__fadd_rz(pfgk[index], pfgk[index-SZF]), 0.5), pfbz[index]);
    }
}

static __global__ void krnl_7(float *pfbx, float *pfby, float *pfbz, float *pfdv,
                        float *pfps, float *pfpt, float *pfu, float *FPS,
                        float *pfCs, float *pfCt, float tcc,
                        int iNx, int iNy, int iNz, int SZF, 
                        int blocksInY, float invBlocksInY)
{
    int blockIdxz = __float2uint_rd(blockIdx.y * invBlocksInY);
    int blockIdxy = blockIdx.y - __umul24(blockIdxz,blocksInY);
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdxy,blockDim.y)+threadIdx.y;
    int idz   = __mul24(blockIdxz,blockDim.z)+threadIdx.z;
    float fpt;

    if( idx<iNy && idy<iNx && idz<iNz)
    {
     int index = idx + __mul24(idy, iNy) + __mul24(idz, SZF);

     pfdv[index] = __fadd_rz(__fadd_rz(__fadd_rz(__fadd_rz(__fadd_rz(
                pfbx[index+iNy], - pfbx[index]), 
                pfby[index+1]), - pfby[index]),
                pfbz[index+SZF]), - pfbz[index]);

     fpt = __fadd_rz(__fadd_rz(__fadd_rz(pfpt[index], - pfu[index]), pfdv[index]), tcc);
     pfps[index] = MIN(fpt, pfCs[index]);

     fpt = __fadd_rz(__fadd_rz(pfps[index], pfu[index]), - pfdv[index]);
     pfpt[index] = MIN(fpt, pfCt[index]);

     fpt = __fadd_rz(__fadd_rz(pfpt[index], pfdv[index]), - pfps[index]);
     FPS[index] = fabsf(fpt);
     pfu[index] = __fadd_rz(pfu[index], -fpt);
    }
}

