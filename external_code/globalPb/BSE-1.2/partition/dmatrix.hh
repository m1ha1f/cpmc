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


#ifndef DMATRIX_H
#define DMATRIX_H

#include <stdio.h>
#include "array.hh"

class DMatrix
{
  public:
    DMatrix(const Util::Array2D<float>& mat);
    ~DMatrix();
    void symmetrize();
    double* getRowSum() const;
    void normalizedLaplacian(const double* rowSum);     //in place, converts into norm laplacian 
    void undoNormalizedLaplacian(const double* rowSum); //in place, converts back into original matrix
    void eig(int nvectors, double** evec, double** eval);
    double computeNCut(const double* rowSum, Util::Array1D<int>& membership) const;

    int n;
    double* values;
};

#endif

