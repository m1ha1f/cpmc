#include <stdio.h>
#include <stdexcept>
#include <cuda.h>
#include <npp.h>
#include <helper_cuda.h>

#include <ImagesCPU.h>
#include <ImagesNPP.h>

#include <mex.h>

inline int cudaDeviceInit()
{
    int deviceCount;
    checkCudaErrors(cudaGetDeviceCount(&deviceCount));

    if (deviceCount == 0)
    {
        fprintf(stderr, "CUDA error: no devices supporting CUDA.\n");
        exit(EXIT_FAILURE);
    }

    int dev = 0;

    if (dev > deviceCount-1)
    {
        fprintf(stderr, ">> %d CUDA capable GPU device(s) detected. <<\n", deviceCount);
        fprintf(stderr, ">> cudaDeviceInit (-device= %d) is not a valid GPU device. <<\n\n", dev);
        return -dev;
    }
    else
    {
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, dev);
        fprintf(stderr, "cudaSetDevice GPU %d = %s\n", dev, deviceProp.name);
    }

    return dev;
}

int printfNPPinfo(int cudaVerMajor, int cudaVerMinor)
{
    fprintf(stderr, "Getting npp version\n");
    const NppLibraryVersion *libVer   = nppGetLibVersion();

    fprintf(stderr, "NPP Library Version %d.%d.%d\n", libVer->major, libVer->minor, libVer->build);

    int driverVersion, runtimeVersion;
    cudaDriverGetVersion(&driverVersion);
    cudaRuntimeGetVersion(&runtimeVersion);

    fprintf(stderr, "  CUDA Driver  Version: %d.%d\n", driverVersion/1000, (driverVersion%100)/10);
    fprintf(stderr, "  CUDA Runtime Version: %d.%d\n", runtimeVersion/1000, (runtimeVersion%100)/10);

    int bVal = checkCudaCapabilities(cudaVerMajor, cudaVerMinor);
    return bVal;
}

void graphCut(int width, int height, Npp32f *pTerminals, Npp32f *pLeftTransposed,
Npp32f *pRightTransposed, Npp32f *pTop, Npp32f *pBottom, Npp8u *labels)
{
    fprintf(stderr, "Starting cuda graphcut computation on image of size %d x %d...\n\n", width, height);

    cudaDeviceInit();

    // Min spec is SM 1.1 devices
    if (!printfNPPinfo(1, 1))
    {
        fprintf(stderr, "Insufficient Compute Capability (must be >= 1.1)\n");
        cudaDeviceReset();
        exit(EXIT_SUCCESS);
    }


    NppiSize size;
    size.width = width;
    size.height = height;

    //Alocate memory on the device
    Npp32f *d_terminals;
    Npp32f *d_left_transposed, *d_right_transposed;
    Npp32f *d_top, *d_bottom;
    size_t step, transposed_step;
    
    checkCudaErrors(cudaMallocPitch(&d_terminals, &step, width*sizeof(Npp32f), height));
    checkCudaErrors(cudaMallocPitch(&d_top, &step, width*sizeof(Npp32f), height));
    checkCudaErrors(cudaMallocPitch(&d_bottom, &step, width*sizeof(Npp32f), height));
    checkCudaErrors(cudaMallocPitch(&d_left_transposed, &transposed_step, height*sizeof(Npp32f), width));
    checkCudaErrors(cudaMallocPitch(&d_right_transposed, &transposed_step, height*sizeof(Npp32f), width));

    fprintf(stderr, "Malloc done\n");


    //Copy capacities to device
    checkCudaErrors(cudaMemcpy2D(d_terminals, step, pTerminals, width * sizeof(Npp32f), width*sizeof(Npp32f), height, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy2D(d_top,       step, pTop,       width * sizeof(Npp32f), width*sizeof(Npp32f), height, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy2D(d_bottom,    step, pBottom,    width * sizeof(Npp32f), width*sizeof(Npp32f), height, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy2D(d_left_transposed,  transposed_step, pLeftTransposed, height * sizeof(Npp32f), height*sizeof(Npp32f), width, cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy2D(d_right_transposed, transposed_step, pRightTransposed, height * sizeof(Npp32f), height*sizeof(Npp32f), width, cudaMemcpyHostToDevice));

    fprintf(stderr, "cudaMemcpy2d done\n");


    // Allocate temp storage for graphcut computation
    Npp8u *pBuffer;
    int bufferSize;
    nppiGraphcutGetSize(size, &bufferSize);
    checkCudaErrors(cudaMalloc(&pBuffer, bufferSize));

    NppiGraphcutState *pGraphcutState;
    nppiGraphcutInitAlloc(size, &pGraphcutState, pBuffer);

     // Allocate label storage
    npp::ImageNPP_8u_C1 oDeviceDst(width, height);

    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // Compute the graphcut, result is 0 / !=0
    cudaEventRecord(start,0);

    fprintf(stderr, "Calling nppiGraphcut_32f8u\n");
    NPP_CHECK_NPP(nppiGraphcut_32f8u(d_terminals, d_left_transposed, d_right_transposed,
                       d_top, d_bottom, step, transposed_step,
                       size, oDeviceDst.data(), oDeviceDst.pitch(), pGraphcutState));
    fprintf(stderr, "nppiGraphcut_32f8u done\n");
    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);

    float time;
    cudaEventElapsedTime(&time, start, stop);
    fprintf(stderr, "Elapsed Time:  %f ms\n", time);

    // declare a host image object for an 8-bit grayscale image
    npp::ImageCPU_8u_C1 oHostAlpha(width, height);

    // convert graphcut result to 0/255 alpha image using new nppiCompareC_8u_C1R primitive (CUDA 5.0)
    npp::ImageNPP_8u_C1 oDeviceAlpha(width, height);
    NPP_CHECK_NPP(nppiCompareC_8u_C1R(oDeviceDst.data(), oDeviceDst.pitch(), 0, oDeviceAlpha.data(), oDeviceAlpha.pitch(), size,
                        NPP_CMP_GREATER));

    // and copy the result to host
    oDeviceAlpha.copyTo(oHostAlpha.data(), oHostAlpha.pitch());

    labels = new unsigned char[width*height];
    for (int i = 0; i < height; ++i)
        for (int j = 0; j < width; ++j)
            labels[i*width+j] = *oHostAlpha.data(i, j);

    checkCudaErrors(cudaFree(d_terminals));
    checkCudaErrors(cudaFree(d_top));
    checkCudaErrors(cudaFree(d_bottom));
    checkCudaErrors(cudaFree(d_left_transposed));
    checkCudaErrors(cudaFree(d_right_transposed));
    checkCudaErrors(cudaFree(pBuffer));
    nppiGraphcutFree(pGraphcutState);

    checkCudaErrors(cudaEventDestroy(start));
    checkCudaErrors(cudaEventDestroy(stop));
}

template<typename T>
T *transpose(T* mat, int rows, int cols)
{
    T* tmp = new T[rows*cols];
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            tmp[i*cols+j] = mat[j*rows+i];
    return tmp;
}

template<typename T>
void transposeInPlace(T* mat, int rows, int cols)
{
    T* tmp = new T[rows*cols];
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            tmp[i*cols+j] = mat[j*rows+i];
    memcpy(mat, tmp, sizeof(T)*rows*cols);
    delete [] tmp;
}

extern void mexFunction(int iNbOut, mxArray *pmxOut[],
    int iNbIn, const mxArray *pmxIn[])
{
    int width = (int)mxGetScalar(pmxIn[0]);
    int height = (int)mxGetScalar(pmxIn[1]);

    Npp32f *pTerminals = transpose((Npp32f*)mxGetData(pmxIn[2]), height, width);

    Npp32f* pLeftTransposed = transpose((Npp32f*)mxGetData(pmxIn[3]), width, height);
    for (int j = 0; j < height; ++j)
        if (pLeftTransposed[j] != 0)
            throw std::invalid_argument("pLeftTransposed[0][*] must be 0");

    Npp32f* pRightTransposed = transpose((Npp32f*)mxGetData(pmxIn[4]), width, height);
    for (int j = 0; j < height; ++j)
        if (pRightTransposed[(width-1)*height + j] != 0)
            throw std::invalid_argument("pRightTransposed[width-1][*] must be 0");  

    fprintf(stderr, "Assertions passed\n");  

    Npp32f* pTop = transpose((Npp32f*)mxGetData(pmxIn[5]), height, width);
    for (int j = 0 ; j < width; ++j)
        if (pTop[j] != 0)
            throw std::invalid_argument("pTop[0][*] must be 0");

    Npp32f* pBottom = transpose((Npp32f*)mxGetData(pmxIn[6]), height, width);
    for (int j = 0; j < width; ++j)
        if (pBottom[width*(height-1) + j] != 0)
            throw std::invalid_argument("pBottom[height-1][*] must be 0"); 

    pmxOut[0] = mxCreateNumericMatrix(height, width, mxINT8_CLASS, mxREAL);
    Npp8u *outmat = (Npp8u*)mxGetData(pmxOut[0]);

    graphCut(width, height, pTerminals, pLeftTransposed, pRightTransposed, pTop, pBottom, outmat);
    transposeInPlace(outmat, height, width); 

    delete [] pTerminals;
    delete [] pLeftTransposed;
    delete [] pRightTransposed;
    delete [] pTop;
    delete [] pBottom;
}
