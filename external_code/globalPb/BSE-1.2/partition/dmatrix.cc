// Copyright (C) 2002 Charless C. Fowlkes <fowlkes@eecs.berkeley.edu>
// Copyright (C) 2002 David R. Martin <dmartin@eecs.berkeley.edu>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
// 02111-1307, USA, or see http://www.gnu.org/copyleft/gpl.html.


#include <iostream>
#include <stdio.h>
#include <math.h>
#include "array.hh"
#include "dmatrix.hh"

//DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO)
extern "C" void dgesvd_ ( char* jobu, char* jobvt,
       int* m, int* n, double* A, int* lda,
       double* S, double* U, int* ldu, double* V, int* ldv,
       double* work, int* lwork, int* info
     );
                                                                                                                 
DMatrix::DMatrix (const Util::Array2D<float>& mat)
{
    assert(mat.size(0) == mat.size(1));
    n = mat.size(0);
    std::cerr << "dmatrix n=" << n << std::endl;
    values = new double[n*n];
    for (int j = 0; j < n; j++)
    {
      for (int i = 0; i < n; i++)
      {
        values[j*n + i] = mat(i,j);
      } 
    }
}

DMatrix::~DMatrix ()
{
  delete[] values;
}

void DMatrix::symmetrize()
{
  for (int j = 0; j < n; j++)
  {
    for (int i = j+1; i < n; i++)
    {
      values[j*n + i] = (values[j*n + i] + values[i*n + j]) / 2;
      values[i*n + j] = values[j*n + i];
    } 
  }
}

double*
DMatrix::getRowSum () const
{
    double* rowSum = new double [n];
    memset(rowSum,0,n*sizeof(double));
    for (int row = 0; row < n; row++) 
    {
      for (int col = 0; col < n; col++)
      {
        rowSum[row] += values[col*n + row];
      }
    }
    return rowSum;
}

void
DMatrix::normalizedLaplacian(const double* rowSum) 
{
    double* isrd = new double[n];
    for (int i = 0; i < n; i++) {
      isrd[i] = 1.0 / sqrt(rowSum[i]);
    }
    for (int row = 0; row < n; row++) 
    {
      for (int col = 0; col < n; col++)
      {
        values[col*n+row] = isrd[col]*values[col*n+row]*isrd[row];
      } 
    }
    delete[] isrd;
}

void
DMatrix::undoNormalizedLaplacian(const double* rowSum) 
{
    double* srd = new double[n];
    for (int i = 0; i < n; i++) {
      srd[i] = sqrt(rowSum[i]);
    }
    for (int row = 0; row < n; row++) 
    {
      for (int col = 0; col < n; col++)
      {
        values[col*n+row] = srd[col]*values[col*n+row]*srd[row];
      } 
    }
    delete[] srd;
}


void 
DMatrix::eig(int nvectors, double** evec, double** eval)
{
  //set up job
  char jobu = 'S';                      //fill out NxN elements of U since matrix is symmetric
  char jobvt = 'S';                     //dont touch V'
  int m = n;                            //rows of A
  int lda = n;                          //leading dimension of A
  double* LA = new double[n*n];         //LA gets overwritten
  for (int row = 0; row < n; row++)     //LAPACK expects column major 
  {
    for (int col = 0; col < n; col++)
    {
      LA[col*n+row] = values[col*n+row]; 
    } 
  }

  double* LS = new double[n];            //singular values 
  int ldu = n;
  double* LU = new double[ldu*n];        //left singular vectors
  int ldvt = n;
  double* LVt = new double[ldvt*n];      //right singular vectors, not computed
  int lwork = 10*n;
  double* work = new double[lwork];
  int info;
                                                                                                                 
  //invoke LAPACK
  dgesvd_ (&jobu,&jobvt,&m,&n,LA,&lda,LS,LU,&ldu,LVt,&ldvt,work,&lwork,&info);
  if (info != 0)
  {
    std::cerr << std::endl << "[LAPACK ERROR] dgesvd returned " << info << std::endl;
  }
   
  //copy out return values
  *eval = new double[nvectors];
  memset(*eval,0,nvectors*sizeof(double));

  // Allocate space for computed eigenvectors.
  *evec = new double [n * nvectors];
  memset(*evec,1,n*nvectors*sizeof(double));

  for (int vec = 0; vec < nvectors; vec++)
  {
    (*eval)[vec] = LS[vec]; 
    for (int row = 0; row < n; row++)
    {
      (*evec)[vec*n+row] = LU[vec*n+row];
    }
  }

  //cleanup
  delete[] LA;
  delete[] LS;
  delete[] LU;
  delete[] LVt;
  delete[] work;
}



double
DMatrix::computeNCut(const double* rowSum, Util::Array1D<int>& membership) const
{
  const int nsegs = membership.size(0);
  if (nsegs > 1)
  {
    double* num = new double[nsegs];
    double* den = new double[nsegs];
    memset(num,0,nsegs*sizeof(double));
    memset(den,0,nsegs*sizeof(double));
    for (int row = 0; row < n; row++)
    {
      int segi = membership(row);
      den[segi] += rowSum[row]; 
      for (int col = 0; col < n; col++)
      {
        int segj = membership(col);
        if (segi == segj)
        {
          num[segi] += values[col*n+row]; 
        }
      }
    }

    double assoc = 0;
    for (int s = 0; s < nsegs; s++)
    {
      if (den[s] > 0)
      {
        assoc += (num[s] / den[s]); 
      }
    }
    return (1 - ((1/(double)nsegs)*assoc));
  }
  else
  {
    std::cerr << "only 1 segment" << std::endl;
    return 0;
  }
}

