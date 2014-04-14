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


#ifndef HISTOGRAM_HH 
#define HISTOGRAM_HH 

#include "array.hh"

//
// histogram code.  common to both color.cc and texture.cc
//

namespace Group
{
  typedef Util::Array3D<float> Histogram;  // (x,y,binnumber)

  //
  // L1 normalized a histogram so that the sum of bins at
  // each location is 1.
  //
  void normalizeHistogram(Histogram& histogram);

  //
  // given two points and the histograms, compute
  // chisquared distance between them.  the histogram
  // should already be L1 normalized
  //
  void chiSquared(const Histogram& histogram, 
                  const int x1, const int y1,
                  const int x2, const int y2, 
                  float& similarity);

}; //namespace Group

#endif
