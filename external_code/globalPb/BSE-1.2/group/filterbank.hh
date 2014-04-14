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


#ifndef FILTERBANK_HH
#define FILTERBANK_HH

#include "array.hh"
#include "image.hh"

namespace Group
{
  //
  // gaussian filterbank class.  even/odd quadrature pairs
  // at some set of scales and orientations.
  //

  class FilterBank
  {
    public:

      FilterBank ();

      //create a gaussian filterbank with the given number of orientations/scales
      FilterBank (const int numOrientations, const int numScales,
                  const float startScale, const float scaleFactor);
      ~FilterBank ();


      //run the filterbank on an image and return
      //a stack of images
      void filter(const Util::Image& im, Util::ImageStack& filtered) const;

      // filtered is the result of running the filterbank
      // returns OE at all scales.
      void orientationEnergy (const Util::ImageStack& filtered, Util::ImageStack& energy) const;

      //create a single 2nd derivative gaussian kernel
      static void createGaussianKernel (const float sigmaX, const float sigmaY,
                                   const float support, const float orient,
                                   const int deriv, const bool hilbert,
                                   Util::Image& kernelImage);

      //create a nice image of the filterbank for presentation or debugging
      void getFilterBankImage(Util::Image& fbim) const;

      //return a const reference to a particular filter
      const Util::Image& getFilter(const int i) const;

      //retrieve properies of the filter bank
      int getNumFilters() const;
      int getNscales() const;
      int getNorient() const;
      float getSigmaX(int index) const;
      float getSigmaY(int index) const;
      float getOrient(int index) const;
      int getIndexElong (int scale, int orient, int evenodd) const;
      int getIndexCenter (int scale) const;

      // abstract indexing based on scale,orientation,evenodd
      static int calcNumFilters(const int nscales, const int norient);
      static int calcIndexElong(const int nscales, const int norient,
                                 const int scale, const int orient,
                                 const int evenodd);
      static int calcIndexCenter(const int nscales, const int norient,
                                 const int scale);

    protected:
      int m_nfilters;
      int m_nscales;
      int m_norient;

      //an array of filters
      Util::Array1D<Util::Image> m_filters;

      //filter parameters
      Util::Array1D<float> m_sigmasX;
      Util::Array1D<float> m_sigmasY;
      Util::Array1D<float> m_orients;
  };

} //namespace Group

#endif
