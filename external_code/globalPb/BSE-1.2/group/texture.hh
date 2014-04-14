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


#ifndef TEXTURE_H
#define TEXTURE_H

#include "array.hh"
#include "image.hh"
#include "ic.hh"
#include "histogram.hh"
#include "filterbank.hh"
#include "string.hh"

namespace Group
{
  //
  // code for finding textons, texture gradients and 
  // comparing texture patches 
  //
  struct TextonMap
  {
    Util::Array2D<int> map;         
    int numChannels;
  };

  //
  // given an image, filterbank and a universal texton file, this function
  // returns a texton map
  //
  void computeTextons(const Util::Image& im, const Group::FilterBank& fbank, 
                      const Util::String& textonFile, TextonMap& textons);

  //
  // given an image and filterbank, this function returns a texton map
  //
  void computeTextons(const Util::Image& im, const Group::FilterBank& fbank,
                       const int numChannels, TextonMap& textons);

  //
  // given a set of image filter responses, this function returns a texton map 
  //
  void computeTextons(const Util::ImageStack& hypercolumn, const int numChannels, TextonMap& textons);

  //
  //compute texture gradient using chi-squared
  //
  void computeTG(const TextonMap& textons, const float radius, 
                 const int norient, Util::ImageStack& tg);

  //
  // given the texton channel membership image and a scale, compute the texton
  // histogram associated with each point 
  //
  void computeTextonHistograms (const TextonMap& textons, const float scale, Histogram& hist);

  //
  // given the texton channel membership image and a scale, compute the texton
  // histogram associated with each point 
  //
  void computeTextonHistograms (const TextonMap& textons, const SupportMap& supportMap, Histogram& hist);

  //
  // given two points and the stack of texture histograms, this 
  // method computes the chisquared distance between them
  //
  void textureSimilarity(const Histogram& textureHistogram,
                          const int x1, const int y1,
                          const int x2, const int y2, float& similarity);

} //namespace group

#endif
