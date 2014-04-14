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
#include "smatrix.hh"
#include "util.hh"
#include "string.hh"
#include "exception.hh"
#include "message.hh"

SMatrix::SMatrix (int n, int* nz, int** col, double** values)
{
    this->n = n;
    this->nz = nz;
    this->col = col;
    this->values = values;
    int nnz = 0;
    for (int i = 0; i < n; i++)
      nnz += nz[i];
    Util::Message::debug(Util::String("creating sparse matrix with %d nonzero entries",nnz));
}

SMatrix::SMatrix(FILE* fp) 
{
  fread(&n,sizeof(int),1,fp);
  int nnz;
  fread(&nnz,sizeof(int),1,fp);
  nz = new int[n];
  col = new int*[n];
  values = new double*[n];

  for (int row = 0; row < n; row++) {
    fread(&nz[row],sizeof(int),1,fp);
    col[row] = new int[nz[row]]; 
    values[row] = new double[nz[row]]; 
    fread(values[row],sizeof(double),nz[row],fp);
    fread(col[row],sizeof(int),nz[row],fp);
  }
}

SMatrix::~SMatrix ()
{
  for (int i = 0; i < n; i++)
  {
    delete[] col[i];
    delete[] values[i];
  }
  delete col;
  delete values;
  delete nz;
}

void 
SMatrix::dump(FILE* fp) 
{
  fwrite(&n,sizeof(int),1,fp);
  int nnz = getNNZ();
  fwrite(&nnz,sizeof(int),1,fp);
  for (int row = 0; row < n; row++) {
    fwrite(&nz[row],sizeof(int),1,fp);
    fwrite(values[row],sizeof(double),nz[row],fp);
    fwrite(col[row],sizeof(int),nz[row],fp);
  }
}

int SMatrix::getNNZ() const
{
  int nnz = 0;
  for (int i = 0; i < n; i++)
    nnz += nz[i];
  return nnz;
}

//older, slower symmetrize
/*
void SMatrix::symmetrize()
{
  for (int r = 0; r < n; r++) 
  {
    for (int i = 0; i < nz[r]; i++) 
    {
      int c = col[r][i];
      int j = -1;
      for (int k = 0; k < nz[c]; k++)
      {
        if (col[c][k] == r)
        {
          j = k;   
        }
      }
      assert(j != -1);
      double v_rc = values[r][i];
      double v_cr = values[c][j];
      values[r][i] = 0.5*(v_rc+v_cr);
      values[c][j] = 0.5*(v_rc+v_cr);
    }
  }  
}
*/


void SMatrix::symmetrize()
{
  int* tail = new int[n];  
  memset(tail,0,n*sizeof(int));
  for (int r = 0; r < n; r++) 
  {
    int offset = 0;
    while ((offset < nz[r]) && (col[r][offset] < r+1))
    {
      offset++;
    }
    for (int i = offset; i < nz[r]; i++) 
    {
      int c = col[r][i];
      assert( col[c][tail[c]] == r ); 
      double v_rc = values[r][i];
      double v_cr = values[c][tail[c]];
      values[r][i] = 0.5*(v_rc+v_cr);
      values[c][tail[c]] = 0.5*(v_rc+v_cr);
      tail[c]++;
    }
  }  
}


double*
SMatrix::getRowSum () const
{
    double* rowSum = new double[n];
    memset(rowSum,0,n*sizeof(double));
    Util::Message::debug("computing row sum");
    for (int row = 0; row < n; row++) {
      double sum = 0;;
      for (int j = 0; j < nz[row]; j++) {
        sum += values[row][j];
      }
      rowSum[row] = sum;
    }
    return rowSum;
}

/*
double
SMatrix::computeNCut(const double* rowSum, const int* membership, const int nsegs) const
{
  if (nsegs > 1)
  {
    double* num = new double[nsegs];
    double* den = new double[nsegs];
    memset(num,0,nsegs*sizeof(double));
    memset(den,0,nsegs*sizeof(double));
    for (int row = 0; row < n; row++)
    {
      int segi = membership[row];
      den[segi] += rowSum[row]; 
      for (int j = 0; j < nz[row]; j++)
      {
        int segj = membership[col[row][j]];
        if (segi == segj)
        {
          num[segi] += values[row][j]; 
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
*/


double
SMatrix::computeNCut(const double* rowSum, const Util::Array1D<int> membership, const int nsegs) const
{
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
      for (int j = 0; j < nz[row]; j++)
      {
        int segj = membership(col[row][j]);
        if (segi == segj)
        {
          num[segi] += values[row][j]; 
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
    delete num;
    delete den;
    return (1 - ((1/(double)nsegs)*assoc));
  }
  else
  {
    Util::Message::debug("only 1 segment!!");
    return 0;
  }
}


void
SMatrix::normalizedLaplacian(const double* rowSum) 
{
    double* isrd = new double[n];
    for (int i = 0; i < n; i++) {
      isrd[i] = 1.0 / sqrt(rowSum[i]);
    }
    Util::Message::debug("scaling rows");
    for (int row = 0; row < n; row++) {
      double isrdrow = isrd[row];
      for (int j = 0; j < nz[row]; j++) {
        values[row][j] = isrdrow * values[row][j] * isrd[col[row][j]];
      }
    }
    delete[] isrd;
}

void
SMatrix::undoNormalizedLaplacian(const double* rowSum) 
{
    double* isrd = new double[n];
    for (int i = 0; i < n; i++) {
      isrd[i] = sqrt(rowSum[i]);
    }
    Util::Message::debug("scaling rows");
    for (int row = 0; row < n; row++) {
      double isrdrow = isrd[row];
      for (int j = 0; j < nz[row]; j++) {
        values[row][j] = isrdrow * values[row][j] * isrd[col[row][j]];
      }
    }
    delete[] isrd;
}

void 
SMatrix::mvmul (const double* a, double* b) const
{
    for (int row = 0; row < n; row++) {
      double bval = 0;
      for (int j = 0; j < nz[row]; j++) {
          bval += a[col[row][j]] * values[row][j];
      }
      b[row] = bval;
    }
}

void 
SMatrix::mvmul (const double* a1, const double* a2, 
			double* b1, double* b2) const
{
    for (int row = 0; row < n; row++) {
      double bval1 = 0;
      double bval2 = 0;
      for (int j = 0; j < nz[row]; j++) {
          double v = values[row][j];
          bval1 += a1[col[row][j]] * v;
          bval2 += a2[col[row][j]] * v;
      }
      b1[row] = bval1;
      b2[row] = bval2;
    }
}

void 
SMatrix::mvmul (const double* a1, const double* a2, 
			const double* a3, const double* a4, 
			double* b1, double* b2,
			double* b3, double* b4) const
{
    for (int row = 0; row < n; row++) {
      double bval1 = 0;
      double bval2 = 0;
      double bval3 = 0;
      double bval4 = 0;
      for (int j = 0; j < nz[row]; j++) {
          double v = values[row][j];
          bval1 += a1[col[row][j]] * v;
          bval2 += a2[col[row][j]] * v;
          bval3 += a3[col[row][j]] * v;
          bval4 += a4[col[row][j]] * v;
      }
      b1[row] = bval1;
      b2[row] = bval2;
      b3[row] = bval3;
      b4[row] = bval4;
    }
}


typedef void (*OpFunc) (int* nrow, int* ncol, double* xin, 
			int* ldx, double* yout, int* ldy);

extern "C" void trlan77_ (OpFunc op, int ipar[32], int* nrow1, int* mev,
			  double* eval, double* evec,
			  int* nrow2, double* wrk, int* lwrk);

static SMatrix* opMatrix = NULL;

static void
op (int* nrow, int* ncol, double* xin, int* ldx, double* yout, int* ldy)
{
    Util::Message::stepBlock();
    assert (opMatrix != NULL);
    assert (*nrow == (int)opMatrix->n);
    int col = 0;
    for (; col < *ncol - 3; col += 4) {
       opMatrix->mvmul (&xin[col*(*ldx)], &xin[(col+1)*(*ldx)], &xin[(col+2)*(*ldx)], &xin[(col+3)*(*ldx)],
                         &yout[col*(*ldy)], &yout[(col+1)*(*ldy)], &yout[(col+2)*(*ldy)], &yout[(col+3)*(*ldy)]);
    }
    for (; col < *ncol - 1; col += 2) {
       opMatrix->mvmul (&xin[col*(*ldx)], &xin[(col+1)*(*ldx)], &yout[col*(*ldy)], &yout[(col+1)*(*ldy)]);
    }
    for (; col < *ncol; col++) {
       opMatrix->mvmul (&xin[col*(*ldx)], &yout[col*(*ldy)]);
    }

}

void 
SMatrix::eigs(int nvectors, double** eval, double** evec)
{
    opMatrix = this;

    int stat = 0;	                  // Error status.
    int nrow = this->n;             // Matrix size.
    int lohi = 1;                   // -1 == small , 1 == large eigenvalues?
    int ned = nvectors;             // Number of eigenpairs to compute.
    int maxlan = 0;	                // Maximum Lanczos basis size, set below.
    int mev = ned;                  // Number of eigenpairs that can be returned.
    int maxmv = 1000*nvectors;	      // Max number of MV muls allowed.
    int restart = 1;	              // Restarting scheme.
    int nec = 0;	                  // Number of eigenpairs already converged.
    int iguess = -1;	              // Initial guess (0=all ones; -1=ones + perturb
    int cpflag = 0;	                // Checkpoint flag.
    int cpio = 98;	                // Fortran IO unit number for checkpoint files.
    int mpicom = 0;	                // MPI communicator.
    int logio = 99;	                // Fortran IO unit number for debug messages.
    int verbose = 1;	              // Verbose level.
    double tol = 1e-8;	            // Relative tolerance on residual norms.
    int mvop = 2*this->getNNZ();      // FLOPs per MV mul.

    int ipar[32] = {stat, lohi, ned, nec, maxlan, restart, maxmv, mpicom,
		    verbose, logio, iguess, cpflag, cpio, mvop};

    // Allocate space for computed eigenvalues.
    *eval = new double[mev];
    memset(*eval,0,mev*sizeof(double));

    // Allocate space for computed eigenvectors.
    *evec = new double [nrow * mev];
    memset(*evec,1,nrow*mev*sizeof(double));

    // Allocate workspace for TRLan.  Rule of thumb given 
    // in Section 6.1 of TRLan manual for minimum size of maxlan.  
    maxlan = Util::max (maxlan, ned + Util::min (6, ned));
    int lwrk = maxlan * (maxlan + 10);
    double* wrk = new double [lwrk];

    // Put tolerance at beginning of workspace!
    wrk[0] = tol;

    Util::Message::startBlock(maxmv,"trlan");
    trlan77_ (op, ipar, &nrow, &mev, *eval, *evec, &nrow, wrk, &lwrk);
    Util::Message::endBlock();

    // cleanup.
    delete [] wrk;

    stat = ipar[0];
    nec = ipar[3];
    if (stat == 0 && nec >= ned) 
    {
      std::cerr << "  trlan exited normally" << std::endl;
    } 
    else 
    {
        throw Util::Exception (Util::String (
        "TRLan error code %d.  See the TRLan User Guide for a "
        "description of this error condition and tips for how to "
        "fix it.\n", stat));
    }

    //reverse ordering of the eigenvectors so the
    //first eigenvector has the largest eigenvalue
    //and make them unit norm
    double* teval = new double[mev];
    double* tevec = new double[nrow*mev];
    for (int i = 0; i < mev; i++)
    {
      teval[i] = (*eval)[mev-i-1];
      float norm = 0.0f;
      for (int j = 0; j < nrow; j++)
      {
        tevec[i*nrow+j] = (*evec)[(mev-i-1)*nrow+j];
        norm += tevec[i*nrow+j]*tevec[i*nrow+j];
      }
      norm = sqrt(norm);
      for (int j = 0; j < nrow; j++)
      {
        tevec[i*nrow+j] /= norm;
      }
    }
    delete *evec;
    delete *eval;
    *evec = tevec;
    *eval = teval;

    std::cerr << "  smallest eigenvalue: " << (*eval)[mev-1] << std::endl;
}

