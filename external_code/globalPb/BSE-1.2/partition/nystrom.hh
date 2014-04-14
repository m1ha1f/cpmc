#ifndef __nystrom_h__
#define __nystrom_h__

//given A and B, computes normalized A and B.
//return shoud be preallocated
void approxNormalize(const double* A, const double* B, const int nsamples, 
            const int nothers, double* NA, double* NB);

//given NA and NB, approximate embedding vectors 
//return shoud be preallocated
//
// use TwoStep technique
void approxEigenvectors(const double* NA, const double* NB, 
                        const int nsamples, const int nothers,
                        double* evec, double* eval);

//// linear algebra support

//pseudo-inverse LAPACK wrapper
// Ainv shoud be preallocated
void pinv(const double* A, int n, int m, double* Ainv); 

//svd LAPACK wrapper
// [U,S,V] shoud be preallocated
void svd(const double* A, int n, int m, double* U, double* S, double* V);

//matrix zero
void mz(double* A, const int An, const int Am);

//matrix copy 
void mc(const double* A, const int An, const int Am,
          double* B, const int Bn, const int Bm);

//matrix transpose
void mt(const double* A, const int An, const int Am,
          double* B, const int Bn, const int Bm);

//matrix matrix addition, C should be preallocated
void mma(const double* A, const int An, const int Am,
     const double* B, const int Bn, const int Bm,
           double* C, const int Cn, const int Cm);

//matrix matrix multiply BLAS wrapper, C should be preallocated
void mmm(const double* A, const int An, const int Am,
     const double* B, const int Bn, const int Bm,
     double* C, const int Cn, const int Cm);

#endif

