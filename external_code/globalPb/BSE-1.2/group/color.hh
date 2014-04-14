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


#ifndef COLOR_H
#define COLOR_H

#include "image.hh"
#include "texture.hh"

//
// code for computing color gradients and comparing
// color patches
//

namespace Group
{
  //
  // Compute color gradient (CG) for multiple orientations.
  // Uses chi-squared distance for histogram comparision.
  // Color is specified in 3 channels; values must be in [0,1]. 
  // Binning granularity can be different for each channel.
  // Setting the number of bins in a channel to 1 effectivly removes
  // that channel from the computation.
  //
  void  computeCG (const Util::Image& channel,
                   const int nbins,     // number of bins in the resulting histogram
                   const float scale,   // disc radius
                   const int norient,   // number of orientations
                   const float sigma,   // kernel sigma
                   const float support, // kernel support (units of sigma)
                   const float zoom,    // kernel sampling (samples per sigma)
                   Util::ImageStack& cg);
  //
  // precompute a histogram for patch comparison
  //
  void computeColorHistograms (const Util::Image& channel, const int nbins,
                                 const float scale, const float sigma,
                                 const float support, const float zoom,
                                  Histogram& histogram);
  //
  // precompute a histogram for patch comparison
  //
  void computeColorHistograms (const Util::Image& channel, const int nbins, const float sigma,
                               const float support, const float zoom, const SupportMap& supportMap,
                               Histogram& histogram);


  //
  // given two points and the color histograms, compute
  // chisquared distance between them
  //
  void colorSimilarity(const Histogram& histogram, 
                         const int x1, const int y1,
                         const int x2, const int y2, float& similarity);

} //namespace Group

#endif
