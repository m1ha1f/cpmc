#include <math.h>
#include <iostream.h>
#include <fstream.h>
#include <assert.h>
#include <string.h>
#include "common.h" 
#include "nystrom.h" 
#include "kmeans.h" 
#include "matrix.h" 
#include "Util.h" 

#define __USE_BLAS

//DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO)
extern "C" void dgesvd_ ( char* jobu, char* jobvt, 
       int* m, int* n, double* A, int* lda, 
       double* S, double* U, int* ldu, double* V, int* ldv, 
       double* work, int* lwork, int* info
     );

#ifdef __USE_BLAS
  //DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )
  extern "C" void dgemm_ (char* transa, char* transb, int* m, int* n, int* k,
          double* alpha, const double* A, int* lda, const double* B, 
          int* ldb, double* beta, double* C, int* ldc);
#endif


///////////////////////////////////////////////////////////////////////////////


// * normalize [A,B] 
//
//   d1 = sum([A;B'],1);
//   d2 = sum(B,1) + sum(B',1)*pinv(A)*B;
//   d = [d1 d2]';
//   v = sqrt(1./d);
//   NA = A.*(v(1:Nsamples)*v(1:Nsamples)');
//   NB = B.*(v(1:Nsamples)*v(Nsamples+(1:Nothers))');

void approxNormalize(const double* A, const double* B, const int nsamples,
                 const int nothers, double* NA, double* NB)
{
  cerr << "  [a] pseudo-inverting A...";
  double* Ainv = new double[nsamples*nsamples];
  pinv(A,nsamples,nsamples,Ainv);
  cerr << "done" << endl;

  cerr << "  [b] approximating row and column sums" << endl;
  cerr << "    (1) computing ar,br,bc...";
  double* ar = new double[nsamples]; memset(ar,0,nsamples*sizeof(double));
  double* br = new double[nsamples]; memset(br,0,nsamples*sizeof(double));
  double* bc = new double[nothers]; memset(bc,0,nothers*sizeof(double));
  for (int i = 0; i < nsamples; i++)
  {
    for (int j = 0; j < nsamples; j++)
    {
      ar[i] += A[i*nsamples + j];
    }
    for (int j = 0; j < nothers; j++)
    {
      br[i] += B[i*nothers + j];
      bc[j] += B[i*nothers + j];
    }
  }
  cerr << "done" << endl;

  cerr << "    (2) multiplying Ainv * B...";
  double* AinvB = new double[nsamples*nothers];
  mmm(Ainv,nsamples,nsamples,B,nsamples,nothers,AinvB,nsamples,nothers);
  cerr << "done" << endl;

  cerr << "    (3) multiplying br * (Ainv*B)...";
  double* brAinvB = new double[nothers];
  mmm(br,1,nsamples,AinvB,nsamples,nothers,brAinvB,1,nothers);
  cerr << "done" << endl;

  cerr << "    (4) adding ar + br...";
  double* d1 = new double[nsamples];
  mma(ar,nsamples,1,br,nsamples,1,d1,nsamples,1);
  cerr << "done" << endl;

  cerr << "    (5) adding bc + (br*Ainv*B)...";
  double* d2 = new double[nothers];
  mma(bc,nothers,1,brAinvB,nothers,1,d2,nothers,1);
  cerr << "done" << endl;

  cerr << "  [c] normalizing A and B by D^(-1/2)...";
  for (int i = 0; i < nsamples; i++)
  {
    for (int j = 0; j < nothers; j++)
    {
      if (j < nsamples)
      {
        NA[i*nsamples+j] = A[i*nsamples+j] / sqrt(d1[i]*d1[j]);
        NB[i*nothers+j] = B[i*nothers+j] / sqrt(d1[i]*d2[j]);
      }
      else
      {
        NB[i*nothers+j] = B[i*nothers+j] / sqrt(d1[i]*d2[j]);
      }
    }
  } 

  int badcount = 0;
  for (int j = 0; j < nothers; j++)
  {
    if (d2[j] <= 0) 
    {
      badcount++;
    }
  }
  if (badcount > 0)
  {
    cerr << "WARNING: approximation may fail, row sums are negative!" << endl;
    cerr << "[[BADCOUNT " << badcount << "/" << nothers << " (";
    cerr << (100.0*(double)badcount/(double)nothers) << "%) ]] ..." ;
  }

  cerr << "done" << endl;
  cerr << "  [d] cleanup ";
  cerr << " [a] ";
  delete[] Ainv;
  cerr << " [b] ";
  delete[] ar;
  delete[] br;
  delete[] bc;
  delete[] d1;
  delete[] d2;
  delete[] AinvB;
  delete[] brAinvB;
  cerr << " [c] ";
  cerr << "done" << endl;
}

// *approximate eigenvectors using 2-step technique
//
// [U,S,junk]=svd(NA);
// Z = [sqrt(S)*U' pinv(sqrt(S))*U'*NB];
// [U,L,j] = svd(Z*Z');
// Va=Z'*U*pinv(sqrt(S));
// for i = 1:Nsamples-1
//   V(:,i) = Va(:,i+1)./Va(:,1);
// end;
void approxEigenvectors(const double* NA, const double* NB,
               const int nsamples, const int nothers, double* evec, double* eval)
{

  cerr << "  [a] computing svd of normalized A...";
  double* U = new double[nsamples*nsamples];
  double* S = new double[nsamples];
  double* Junk = new double[nsamples*nsamples];
  svd(NA,nsamples,nsamples,U,S,Junk);
  cerr << "done" << endl;

  cerr << "  [b] computing Z = [sqrt(S)*U' pinv(sqrt(S))*U'*B]" << endl;
  cerr << "    (1) computing sqrt(S) and pinv(sqrt(S))...";
  double* sqrtS = new double[nsamples*nsamples];
  double* pinvsqrtS = new double[nsamples*nsamples];
  for (int i = 0; i < nsamples; i++)
  {
    for (int j = 0; j < nsamples; j++)
    {
      if (i == j)
      {
        sqrtS[i*nsamples+j] = sqrt(S[i]);
        if (S[i] != 0)
        {
          pinvsqrtS[i*nsamples+j] = 1 / sqrt(S[i]);
        }
        else
        {
          pinvsqrtS[i*nsamples+j] = 0; 
        }
      }
      else
      {
        sqrtS[i*nsamples+j] = 0;
        pinvsqrtS[i*nsamples+j] = 0; 
      }
    }
  }
  cerr << "done" << endl;
  cerr << "    (2) transposing U... ";
  double* Utrans = new double[nsamples*nsamples];
  mt(U,nsamples,nsamples,Utrans,nsamples,nsamples);
  cerr << "done" << endl;
  cerr << "    (3) multiplying sqrt(S) * U'...";
  double* sqrtSUtrans = new double[nsamples*nsamples];
  mmm(sqrtS,nsamples,nsamples,
     Utrans,nsamples,nsamples,
     sqrtSUtrans,nsamples,nsamples);
  cerr << "done" << endl;
  cerr << "    (4) multiplying pinv(sqrt(S)) * U'...";
  double* pinvsqrtSUtrans = new double[nsamples*nsamples];
  mmm(pinvsqrtS,nsamples,nsamples,
     Utrans,nsamples,nsamples,
     pinvsqrtSUtrans,nsamples,nsamples);
  cerr << "done" << endl;
  cerr << "    (5) multiplying (pinv(sqrt(S))*U')*NB...";
  double* pinvsqrtSUtransNB = new double[nsamples*nothers];
  mmm(pinvsqrtSUtrans,nsamples,nsamples,
      NB,nsamples,nothers,
      pinvsqrtSUtransNB,nsamples,nothers);
  cerr << "done" << endl;
  cerr << "    (6) copying into Z...";
  double* Z = new double[nsamples*(nsamples+nothers)];
  for (int i = 0; i < nsamples; i++)
  {
    for (int j = 0; j < nsamples; j++)
    {
      Z[i*(nsamples+nothers)+j] = sqrtSUtrans[i*nsamples+j];
    }
    for (int j = 0; j < nothers; j++)
    {
      Z[i*(nsamples+nothers)+j+nsamples] = pinvsqrtSUtransNB[i*nothers+j];
    }
  }
  cerr << "done" << endl;
 


  cerr << "  [c] computing [U,L,Junk] = svd(Z*Z')" << endl;
  cerr << "    (1) transposing Z...";
  double* Ztrans = new double[nsamples*(nsamples+nothers)];
  mt(Z,nsamples,(nsamples+nothers),Ztrans,(nsamples+nothers),nsamples);
  cerr << "done" << endl;
  cerr << "    (2) multiplying Z * Z'...";
  double* ZZtrans = new double[nsamples*nsamples];
  mmm(Z,nsamples,(nsamples+nothers),
      Ztrans,(nsamples+nothers),nsamples,
      ZZtrans,nsamples,nsamples);
  cerr << "done" << endl;
  cerr << "    (3) compute svd(Z*Z')...";
  //reuse U,S, and Junk
  svd(ZZtrans,nsamples,nsamples,U,S,Junk);
  cerr << "done" << endl;



  cerr << "  [e] computing Z'*U*pinv(sqrt(S))" << endl;
  cerr << "    (1) computing pinv(sqrt(S))...";
  //reuse pinvsqrtS
  for (int i = 0; i < nsamples; i++)
  {
    for (int j = 0; j < nsamples; j++)
    {
      if ((i == j) && (S[i] != 0))
      {
        pinvsqrtS[i*nsamples+j] = 1 / sqrt(S[i]);
      }
      else
      {
        pinvsqrtS[i*nsamples+j] = 0; 
      }
    }
  }
  cerr << "done" << endl;
  cerr << "    (2) multiplying U * (pinv(sqrt(S)))...";
  double* UpinvsqrtS = new double[nsamples*nsamples];
  mmm(U,nsamples,nsamples,
      pinvsqrtS,nsamples,nsamples,
      UpinvsqrtS,nsamples,nsamples);
  cerr << "done" << endl;
  cerr << "    (3) multiplying Z' * (U*pinv(sqrt(L)))...";
  double* ZtransUpinvsqrtS = evec;
  mmm(Ztrans,(nsamples+nothers),nsamples,
      UpinvsqrtS,nsamples,nsamples,
      ZtransUpinvsqrtS,(nsamples+nothers),nsamples);
  cerr << "done" << endl;

  //copy S into evals
  for (int i = 0; i < nsamples; i++)
  {
    eval[i] = S[i];
  }

  cerr << "  [f] cleanup...";
  cerr << "[a] ";
  delete[] U;
  delete[] S;
  delete[] Junk;
  cerr << "[b] ";
  delete[] sqrtS;
  delete[] pinvsqrtS;
  delete[] Utrans;
  delete[] pinvsqrtSUtrans;
  delete[] pinvsqrtSUtransNB;
  delete[] Z;
  cerr << "[c] ";
  delete[] Ztrans;
  delete[] ZZtrans;
  cerr << "[d] ";
  cerr << "[e] ";
  delete[] UpinvsqrtS;
  //delete[] ZtransUpinvsqrtS;  evec is pre-allocated 
  cerr << " done" << endl;
}


//svd wrapper, [U,S,V] shoud be preallocated
void svd(const double* A, const int rows, const int cols, 
           double* U, double* S, double* V)
{
  //set up job
  char jobu = 'S';                      //fill out U
  char jobvt = 'S';                     //fill out V'
  int m = rows;                          //rows of A
  int n = cols;                          // cols of A
  int lda = m;                           //leading dimension of A
  double* LA = new double[lda*n];        //LA gets overwritten
  mt(A,m,n,LA,n,m);      //LAPACK expects column major so transpose and
                         //lie about the dimensions... :)
  double* LS = new double[Util::min(m,n)];      //singular values of A
  int ldu = m;
  double* LU = new double[ldu*Util::min(m,n)];  //left singular vectors
  int ldvt = Util::min(m,n);
  double* LVt = new double[ldvt*m];       //right singular vectors
  int lwork = 2 * Util::max( 3*Util::min(m,n) + Util::max(m,n) , 5*Util::min(m,n) );
  double* work = new double[lwork];
  int info;

  //invoke LAPACK
  dgesvd_ (&jobu,&jobvt,&m,&n,LA,&lda,LS,LU,&ldu,LVt,&ldvt,work,&lwork,&info);
  if (info != 0)
  {
    cerr << endl << "[LAPACK ERROR] dgesvd returned " << info << endl;
  }

  //copy out return values, transposing so that A = U*diag(S)*V' 
  mt(LU,ldu,Util::min(m,n),U,Util::min(m,n),ldu);
  mc(LS,Util::min(m,n),1,S,Util::min(m,n),1);
  mc(LVt,ldvt,m,V,ldvt,m);

  //cleanup
  delete[] LA;
  delete[] LS;
  delete[] LU;
  delete[] LVt;
  delete[] work;
}

//pseudo-inverse wrapper, return shoud be preallocated
//
//  [m,n] = size(A);
//  [U,S,V] = svd(A,0);
//  s = diag(S);
//  tol = max(m,n) * max(s) * eps;
//  r = sum(s > tol);
//  if (r == 0)
//    X = zeros(size(A'));
//  else
//    Sinv = diag(ones(r,1)./s(1:r));
//    X = V(:,1:r)*Sinv*U(:,1:r)';
//  end
//
//  TODO: currently assumes m < n, this is foolish
//
void pinv(const double* A, const int n, const int m, double* Ainv)
{

  assert(m <= n);

  double* U = new double[n*m];
  double* S = new double[m];
  double* V = new double[m*n];
  svd(A,n,m,U,S,V);

  double maxS = 0;
  for (int i = 0; i < m; i++)
  {
    maxS = Util::max(S[i],maxS);
  }
  double tolerance = maxS * Util::max(m,n) * 2.2204e-14;

  for (int i = 0; i < m; i++)
  {
    if (S[i] > tolerance)
    {
      S[i] = 1.0 / S[i];
    }
    else
    {
      S[i] = 0;  
    }
  }

  double* SinvUtrans = new double[m*n];
  memset(SinvUtrans,0,m*n*sizeof(double));
  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
    {
      SinvUtrans[i*n + j] = S[i]*U[j*m + i];
    }
  }

  memset(Ainv,0,m*m*sizeof(double));
  mmm(V,m,n,SinvUtrans,n,m,Ainv,m,m);    // Ainv = V*(1/S)*U'
}


// A = 0
// matrix zero
void mz(double* A, const int An, const int Am)
{
  for (int i = 0; i < An; i++)
  {
    for (int j = 0; j < Am; j++)
    {
      A[i*Am + j] = 0.0;
    }
  }
  
}


// B = A
// matrix copy 
void mc(const double* A, const int An, const int Am,
              double* B, const int Bn, const int Bm)
{
  assert(An == Bn);
  assert(Am == Bm);
  for (int i = 0; i < An; i++)
  {
    for (int j = 0; j < Am; j++)
    {
      B[i*Bm + j] = A[i*Am + j];
    }
  }
}


// B = A'
// matrix transpose
void mt(const double* A, const int An, const int Am,
              double* B, const int Bn, const int Bm)
{
  assert(Am == Bn);
  assert(An == Bm);
  for (int i = 0; i < An; i++)
  {
    for (int j = 0; j < Am; j++)
    {
      B[j*An + i] = A[i*Am + j];
    }
  }
}


// C = A + B
// matrix matrix addition 
void mma(const double* A, const int An, const int Am, 
         const double* B, const int Bn, const int Bm,
         double* C, const int Cn, const int Cm)
{
  assert(An == Bn);
  assert(An == Cn);
  assert(Am == Bm);
  assert(Am == Cm);
  for (int i = 0; i < Cn; i++)
  {
    for (int j = 0; j < Cm; j++)
    {
      C[i*Cm + j] = A[i*Am + j] + B[i*Bm + j];
    }
  }
}


// C = A*B
// matrix matrix multiply BLAS wrapper 
// 
// TODO: based on the dimensions of A and B, call
// the approriate BLAS function for matrix-matrix
// or matrix-vector multiply.
void mmm(const double* A, const int An, const int Am, 
         const double* B, const int Bn, const int Bm,
         double* C, const int Cn, const int Cm)
{
  assert(Am == Bn);
  assert(Cn == An);
  assert(Cm == Bm);

#ifdef __USE_BLAS 
  //we are using row-major so multiply the transposes...
  char transa = 'N'; char transb = 'N';
  int m = Bm; int n = An; int k = Bn;
  int lda = Am; int ldb = Bm; int ldc = Bm;
  double alpha = 1.0; double beta = 0;
  dgemm_ (&transa, &transb, &m, &n, &k, &alpha, 
          B, &ldb, A, &lda, &beta, C, &ldc);
#else
  memset(C,0,Cn*Cm*sizeof(double));
  for (int i = 0; i < Cn; i++)
  {
    for (int j = 0; j < Cm; j++)
    {
      for (int k = 0; k < Am; k++)
      {
        C[i*Cm + j] += (A[i*Am + k] * B[k*Bm + j]);
      }
    }
  }
#endif
}
