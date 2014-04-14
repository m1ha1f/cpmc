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


#include <math.h>
#include "exception.hh"
#include "message.hh"
#include "util.hh"
#include "array.hh"
#include "image.hh"
#include "color.hh"

namespace Group
{

  // compute high-quality discrete 1D gaussian kernel
  // for histogram smoothing
  static void gaussian1D (const float sigma,      // sigma of gaussian
                          const float support,    // units of sigma
                          const float res,        // bin width
                          Util::Array1D<float>& kernel)    // result
  {
      const int zoom = 100;       // must be even
      const int r = (int) ceil (sigma * support / res);
      const int d = 2 * r + 1;
      kernel.resize (d);
      kernel.init (0);
      float sum = 0;
      for (int i = 0; i < r * zoom + zoom / 2; i++)
      {
        const float x = ((float) i + 0.5) / zoom * res;
        const int index = (int) rint ((float) x / res);
        const float val = exp (-x * x / (2 * sigma * sigma));
        assert (index >= 0 && index <= r);
        kernel (r + index) += val;
        kernel (r - index) += val;
        sum += 2 * val;
      }

      // Unit L1 norm.
      for (int i = 0; i < d; i++)
      {
        kernel (i) /= sum;
      }
  }


  void computeCG(const Util::Image& channel, const int nbins, const float scale, const int norient,
                 const float sigma, const float support, const float zoom, Util::ImageStack& cg)
  {
      // check arguments
      assert (nbins > 0);
      assert (scale > 0);
      assert (norient > 0);
      assert (sigma > 0);
      assert (support > 0);
      assert (zoom > 0);

      const int width = channel.size(0);
      const int height = channel.size(1);

      // all color values must be in [0,1]
      for (int x = 0; x < width; x++)
      {
        for (int y = 0; y < height; y++)
        {
          const float a = channel(x,y);
          assert (a >= 0 && a <= 1);
        }
      }

      // allocate cg
      cg.resize(norient,width,height);

      // construct our gaussian blob that we will add to the 
      // histogram for each sample
      Util::Array1D<float> kernel;
      gaussian1D (sigma, support, sigma / zoom, kernel);
      assert ((kernel.size () % 2) == 1);
      const int kd = kernel.size ();
      const int kr = kd / 2;
      assert (kd == 2 * kr + 1);

      // precompute spatial offsets for the  kernel
      Util::Array1D<float> offsets(kd);
      for (int i = -kr; i <= kr; i++)
      {
        offsets(kr+i) = (float) i / kr / zoom * sigma * support;
      }

      // compute pie slice membership for each window location
      // count the number of pixels inside the disc
      const int windowRadius = (int) ceil (scale);
      const int windowDiameter = 2 * windowRadius;
      Util::Array2D<int> pie (windowDiameter, windowDiameter);
      for (int u = -windowRadius; u < windowRadius; u++)
      {
        for (int v = -windowRadius; v < windowRadius; v++)
        {
          const int i = u + windowRadius;
          const int j = v + windowRadius;
          float theta = atan2(v+0.5,u+0.5);
          if (theta < 0)
          {
            theta += 2 * M_PI;
          }
          assert (theta >= 0 && theta < 2 * M_PI);
          const int slice = (int) floor (theta / (M_PI / norient));
          assert (slice >= 0 && slice < 2 * norient);
          pie(i,j) = slice;
        }
      }

      // left and right disc histograms
      Util::Array1D<float> histL(nbins);
      Util::Array1D<float> histR(nbins);
      Util::Array2D<float> histPie(2*norient,nbins);

      Util::Message::startBlock(width,"CG computation");
      for (int x = 0; x < width; x++)
      {
        Util::Message::stepBlock();
        for (int y = 0; y < height; y++)
        {
          // comute the pie slice histograms
          histPie.init(0);
          for (int u = -windowRadius; u < windowRadius; u++)
          {
            for (int v = -windowRadius; v < windowRadius; v++)
            {
              if ((u+0.5)*(u+0.5) + (v+0.5)*(v+0.5) > scale*scale) { continue; }
              const int xi = x + u;
              const int yi = y + v;
              if (xi < 0 || xi >= width) { continue; }
              if (yi < 0 || yi >= height) { continue; }

              // figure out which pie slice we're in
              const int slice = pie(windowRadius+u,windowRadius+v);

              // get this point's color
              const float a = channel(xi,yi);
              for (int i = 0; i < kd; i++)
              {
                const float ai = a + offsets(i);
                int bin1 = (int) floor (ai * nbins);
                bin1 = Util::minmax(0, bin1, nbins-1);
                histPie(slice, bin1) += kernel(i);
              }
            }
          }

  #define ACC(h,OP,slice) { \
    const int __slice = slice; \
    for (int i = 0; i < nbins; i++) \
      h(i) OP histPie(__slice,i); \
    }
          // set up L histogram (minus last slice)
          histL.init (0);
          for (int orient = 0; orient < norient - 1; orient++)
          {
            ACC(histL,+=,orient);
          }

          // set up R histogram (minus last slice)
          histR.init (0);
          for (int orient = norient; orient < 2 * norient - 2; orient++)
          {
            ACC(histR,+=,orient);
          }

          // spin the disc
          for (int orient = 0; orient < norient; orient++)
          {
            // add next slice into L and R histograms
            ACC(histL,+=,orient+norient-1);
            ACC(histR,+=,Util::wrap(orient-1,2*norient));
            // normalize L,R histograms
            float sumL = 0;
            float sumR = 0;
            for (int i = 0; i < nbins; i++)
            {
              sumL += histL(i);
              sumR += histR(i);
            }
            if (sumL <= 0)
            {
              sumL = 1;
            }
            if (sumR <= 0)
            {
              sumR = 1;
            }
            assert(sumL > 0);
            assert(sumR > 0);

            // compute chi-squared distance
            float chi2 = 0;
            for (int i = 0; i < nbins; i++)
            {
              float LL = (histL(i) / sumL);
              float RR = (histR(i) / sumR);
              float sum = LL + RR;
              if (sum > 0)
              {
                float dif = LL - RR;
                chi2 += dif*dif / sum;
              }
            }
            chi2 *= 0.5f;
            chi2 = Util::minmax (0.0f,chi2,1.0f);
            assert (chi2 >= 0 && chi2 <= 1);
            cg(orient,x,y) = chi2;
            
            //subtract last slice
            if (orient < norient - 1)
            {
              ACC(histL,-=,orient);
              ACC(histR,-=,orient+norient);
            }
          }
        }
  #undef ACC
      }
      Util::Message::endBlock();
  }

  // precompute histogram for a given scale 
  void computeColorHistograms (const Util::Image& channel, const int nbins,
                               const float scale, const float sigma,
                               const float support, const float zoom,
                               Histogram& histogram)
  {
      // check arguments
      assert (nbins > 0);
      assert (scale> 0);
      assert (sigma > 0);
      assert (support > 0);
      assert (zoom > 0);

      const int width = channel.size(0);
      const int height = channel.size(1);

      // all color values must be in [0,1]
      for (int x = 0; x < width; x++)
      {
        for (int y = 0; y < height; y++)
        {
          const float a = channel(x,y);
          assert (a >= 0 && a <= 1);
        }
      }

      // construct our gaussian blob that we will add to the 
      // histogram for each sample
      Util::Array1D<float> kernel;
      gaussian1D (sigma, support, sigma / zoom, kernel);
      assert ((kernel.size () % 2) == 1);
      const int kd = kernel.size ();
      const int kr = kd / 2;
      assert (kd == 2 * kr + 1);

      // precompute spatial offsets for the  kernel
      Util::Array1D<float> offsets(kd);
      for (int i = -kr; i <= kr; i++)
      {
        offsets(kr+i) = (float) i / kr / zoom * sigma * support;
      }

      // allocate histogram 
      histogram.resize(width,height,nbins);
      histogram.init(0);

      // minimum square window radius to include disc
      const int windowRadius = (int)ceil(scale);

      Util::Message::startBlock(width,"color histogram patch computation");
      for (int x = 0; x < width; x++)
      {
        Util::Message::stepBlock();
        for (int y = 0; y < height; y++)
        {
          for (int u = -windowRadius; u <= windowRadius; u++)
          {
            for (int v = -windowRadius; v <= windowRadius; v++)
            {
              if (u*u + v*v > scale*scale) { continue; }
              const int xi = x + u;
              const int yi = y + v;
              if (xi < 0 || xi >= width) { continue; }
              if (yi < 0 || yi >= height) { continue; }

              // get this point's color
              const float a = channel(xi,yi);
              for (int i = 0; i < kd; i++)
              {
                const float ai = a + offsets(i);
                int bin = (int) floor (ai * nbins);
                bin = Util::minmax(0, bin, nbins-1);
                histogram(x,y,bin) += kernel(i);
              }
            }
          }
        }
      }
      Util::Message::endBlock();
      Util::Message::startBlock("histogram normalization");
      normalizeHistogram(histogram);
      Util::Message::endBlock();
  }


  // precompute histogram for a given scale 
  void computeColorHistograms (const Util::Image& channel, const int nbins, const float sigma,
                               const float support, const float zoom, const SupportMap& supportMap,
                               Histogram& histogram)
  {
      // check arguments
      assert (nbins > 0);
      assert (sigma > 0);
      assert (support > 0);
      assert (zoom > 0);

      const int width = channel.size(0);
      const int height = channel.size(1);

      // all color values must be in [0,1]
      for (int x = 0; x < width; x++)
      {
        for (int y = 0; y < height; y++)
        {
          assert (channel(x,y) >= 0 && channel(x,y) <= 1);
        }
      }

      // construct our gaussian blob that we will add to the 
      // histogram for each sample
      Util::Array1D<float> kernel;
      gaussian1D (sigma, support, sigma / zoom, kernel);
      assert ((kernel.size () % 2) == 1);
      const int kd = kernel.size ();
      const int kr = kd / 2;
      assert (kd == 2 * kr + 1);

      // precompute spatial offsets for the  kernel
      Util::Array1D<float> offsets(kd);
      for (int i = -kr; i <= kr; i++)
      {
        offsets(kr+i) = (float) i / kr / zoom * sigma * support;
      }

      // allocate histogram 
      histogram.resize(width,height,nbins);
      histogram.init(0);

      Util::Message::startBlock(width,"weighted color histogram patch computation");
      for (int x = 0; x < width; x++)
      {
        Util::Message::stepBlock();
        for (int y = 0; y < height; y++)
        {

          for (int i = 0; i < supportMap(x,y).size(); i++)
          {
            // compute coordinates of point in the window
            const int ix = supportMap(x,y)(i).x;
            const int iy = supportMap(x,y)(i).y;
            assert(ix >= 0);
            assert(ix < width);
            assert(iy >= 0);
            assert(iy < height);

            // get this point's color and weighting
            const float a = channel(ix,iy);
            const float weight = supportMap(x,y)(i).sim;

            for (int i = 0; i < kd; i++)
            {
              const float ai = a + offsets(i);
              int bin = (int) floor (ai * nbins);
              bin = Util::minmax(0, bin, nbins-1);
              histogram(x,y,bin) += weight*kernel(i);
            }
          }
        }
      }
      Util::Message::endBlock();

      Util::Message::startBlock("histogram normalization");
      normalizeHistogram(histogram);
      Util::Message::endBlock();
  }

  //
  // given two points and the color histograms, compute
  // chisquared distance between them
  //
  void colorSimilarity(const Histogram& histogram,
                         const int x1, const int y1,
                         const int x2, const int y2, float& similarity)
  {
    chiSquared(histogram,x1,y1,x2,y2,similarity);
  }

} //namespace Group


