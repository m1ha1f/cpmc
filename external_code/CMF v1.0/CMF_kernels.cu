#define BLOCK_SIZE 16
#define MAX(a,b) ( a > b ? a : b )
#define MIN(a,b) ( a <= b ? a : b )
#define SIGN(x) ( x >= 0.0 ? 1.0 : -1.0 )
#define ABS(x) ( (x) > 0.0 ? x : -(x) )
#define SQR(x) (x)*(x)

static __global__ void krnl_1(float *pfpt, float *pfps, float *pfu, 
        float *pfgk, float *pfdv, float cc, int iNx, int iNy)
{
   int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
   int idy   = __mul24(blockIdx.y,blockDim.y)+threadIdx.y;

   if( idx<iNy && idy<iNx )
   {
    int index = idx + idy*iNy;
    pfgk[index] = pfdv[index] - (pfps[index] - pfpt[index] + pfu[index]/cc);
   }

}

static __global__ void krnl_2(float *pfbx, float *pfgk, float steps, 
                            int iNx, int iNy)
{
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdx.y,blockDim.y)+threadIdx.y;

    if( idx<iNy && idy<(iNx-1) )
    {
      int index = idx + (idy+1)*iNy; 
      pfbx[index] = steps*(pfgk[index] - pfgk[index-iNy]) + pfbx[index];
    }
}

static __global__ void krnl_3(float *pfby, float *pfgk, float steps, 
                            int iNx, int iNy)
{
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdx.y,blockDim.y)+threadIdx.y;
    
    if( idx<(iNy-1) && idy<iNx )
    {
      int index = idx + idy*iNy + 1;
      pfby[index] = steps*(pfgk[index] - pfgk[index-1]) + pfby[index];
    }
}


static __global__ void krnl_4(float *pfbx, float *pfby, float *pfgk, float *pfpenalty, 
                            int iNx, int iNy)
{
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdx.y,blockDim.y)+threadIdx.y;
    float fpt;
   
    if( idx<iNy && idy<iNx )
    {
      int index = idx + idy*iNy;

      fpt = sqrt((SQR(pfbx[index]) + SQR(pfbx[index+iNy]) 
                      + SQR(pfby[index]) + SQR(pfby[index+1]))*0.5);
                
                if (fpt > pfpenalty[index])
                    fpt = fpt / pfpenalty[index];
                else
                    fpt = 1;
                
                pfgk[index] = 1/fpt;
    }
}

static __global__ void krnl_5(float *pfbx, float *pfgk, int iNx, int iNy)
{
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdx.y,blockDim.y)+threadIdx.y;

    if( idx<iNy && idy<(iNx-1) )
    {
      int index = idx + (idy+1)*iNy;   
      pfbx[index] = (pfgk[index] + pfgk[index-iNy])*0.5*pfbx[index];
    }
}


static __global__ void krnl_6(float *pfby, float *pfgk, 
                        int iNx, int iNy)
{
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdx.y,blockDim.y)+threadIdx.y;
    
    if( idx<(iNy-1) && idy<iNx )
    {
      int index = idx + idy*iNy+1;
      pfby[index] = 0.5*(pfgk[index] + pfgk[index-1])*pfby[index];
    }
}

static __global__ void krnl_7(float *pfbx, float *pfby, float *pfdv, 
                        int iNx, int iNy)
{
   int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
   int idy   = __mul24(blockIdx.y,blockDim.y)+threadIdx.y;

   if( idx<iNy && idy<iNx )
   {
     int index = idx + idy*iNy;
     pfdv[index] = pfbx[index+iNy] - pfbx[index] + pfby[index+1] - pfby[index];
    }
}
 

static __global__ void krnl_8(float *pfps, float *pfpt, float *pfu, float *pfdv, 
                float *pfCs, float cc, int iNx, int iNy)
{
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdx.y,blockDim.y)+threadIdx.y;
    float fpt;
    
    if( idx<iNy && idy<iNx )
    {
      int index = idx + idy*iNy;

      fpt = pfpt[index] - pfu[index]/cc + pfdv[index] + 1/cc;
      fpt = MIN(fpt, pfCs[index]);
      pfps[index] = fpt;
    }
}


static __global__ void krnl_9(float *pfps, float *pfpt, float *pfu, float *pfdv, 
                    float *pfCt, float cc, int iNx, int iNy)
{
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdx.y,blockDim.y)+threadIdx.y;
    float fpt;
    
    if( idx<iNy && idy<iNx )
    {
      int index = idx + idy*iNy;

      fpt = pfps[index] + pfu[index]/cc - pfdv[index];
      fpt = MIN(fpt, pfCt[index]);
      pfpt[index] = fpt;
    }
}


static __global__ void krnl_10(float *pfpt, float *pfdv, float *pfps, float *pfu, 
                        float *FPS, float cc, int iNx, int iNy)
{
    int idx   = __mul24(blockIdx.x,blockDim.x)+threadIdx.x;
    int idy   = __mul24(blockIdx.y,blockDim.y)+threadIdx.y;
    float fpt;

    if( idx<iNy && idy<iNx )
    {
      int index = idx + idy*iNy;

      fpt = cc * (pfpt[index] + pfdv[index] - pfps[index]);
      FPS[index] = ABS(fpt);
      pfu[index] -= fpt;
    }
}