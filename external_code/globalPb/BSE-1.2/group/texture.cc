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


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <values.h>

#include "string.hh"
#include "exception.hh"
#include "array.hh"
#include "image.hh"
#include "util.hh"
#include "message.hh"
#include "kmeans.hh"

#include "texture.hh"


namespace Group
{

  //
  // given an image, filterbank and a universal texton file, this function 
  // returns a texton map
  //
  void computeTextons(const Util::Image& im, const Group::FilterBank& fbank, 
                       const Util::String& textonFile, TextonMap& textons)
  {
      const int numFilters = fbank.getNumFilters();
      textons.map.resize(im.size(0),im.size(1));
                                                                                                                 
      // run the filterbank
      Util::Message::debug("filtering image");
      Util::ImageStack filtered;
      fbank.filter(im,filtered);
      Util::Message::debug("done filtering image");
                                                                                                                 
      // use the universal textons
      Util::Array2D<float> universalTextons;
      Util::Message::debug(Util::String("reading universal textons from %s",textonFile.text()));
                                                                                                               
      std::ifstream is(textonFile);
      if(!is.good())
      {
        throw Util::Exception(Util::String("could not open '%s' for reading.",textonFile.text()));
      }
      is >> universalTextons;
      is.close();
                                                                                                               
      textons.numChannels = universalTextons.size(0);
      assert(universalTextons.size(1) == numFilters);

      // compute the texton map
      for(int x = 0; x < im.size(0); x++)
      {
        for(int y = 0; y < im.size(1); y++)
        {
          int best = -1;
          float minDist = FLT_MAX;
          for(int i = 0; i < textons.numChannels; i++)
          {
            float dist = 0;
            for(int f = 0; f < numFilters; f++)
            {
              const float v = filtered(f,x,y) - universalTextons(i,f);
              dist += v*v;
            }
            if(dist < minDist)
            {
              minDist = dist;
              best = i;
            }
          }
          assert(best != -1);
          textons.map(x,y) =  best;
        }
      }
}                                                                                                               
                                                                                                                 
  //
  // given an image, filterbank and desired number of channels,
  // this function returns a texton map
  //
  void computeTextons(const Util::Image& im, const Group::FilterBank& fbank, 
                       const int numChannels, TextonMap& textons)
  {
    // run the filterbank
    Util::Message::debug("filtering image");
    Util::ImageStack filtered;
    fbank.filter(im,filtered);

    // cluster
    Util::Message::debug(Util::String("clustering filter responses into %d channels",numChannels));
    computeTextons(filtered,numChannels,textons);
  }
 

  // given a set of image filter responses, this function returns a texton map
  void computeTextons(const Util::ImageStack& hypercolumn, const int numChannels, TextonMap& textons)
  {
      int numDim = hypercolumn.size(0);
      int width = hypercolumn.size(1);
      int height = hypercolumn.size(2);
      int numPoints = width*height;

      textons.numChannels = numChannels;
      textons.map.resize(width,height);

      // Transpose the data so that each hypercolumn is contiguous
      Util::Array2D<float> data(width*height, numDim);
      data.init(0);
      for (int x = 0; x < width; x++)
      {
        for (int y = 0; y < height; y++)
        {
          for (int k = 0; k < numDim; k++)
          {
            data(x*height + y,k) = hypercolumn(k,x,y);
          }
        }
      }

      // Run kmeans to get the cluster means and memberships.
      Util::Image means(numChannels, numDim);
      Util::Array1D<int> membership(numPoints);
      Util::kmeans(data, numChannels, 1e-4, 30, true, means, membership);

      //check that membership is well behaved
      for (int x = 0; x < width; x++)
      {
        for (int y = 0; y < height; y++)
        {
          assert(membership(x*height + y) >= 0);    
          assert(membership(x*height + y) < numChannels);    
          textons.map(x,y) = membership(x*height+y);
        }
      }
  }


  // compute TG at all orientations at a single scale
  void computeTG (const TextonMap& textons, const float radius, 
                   const int norient, Util::ImageStack& tg)
  {
      assert(textons.numChannels > 0);
      assert(radius > 0);
      assert(norient > 0);

      const int width = textons.map.size(0);
      const int height = textons.map.size(1);

      tg.resize(norient,width,height);
      tg.init(0);

      // pre-compute pie slice membership for each window location
      // also count the number of pixels inside the disc
      int pixelCount = 0;
      const int windowRadius = (int)ceil(radius);
      const int windowDiam = 2 * windowRadius;
      Util::Array2D<int> pie (windowDiam, windowDiam);
      for (int u = -windowRadius; u < windowRadius; u++)
      {
        for (int v = -windowRadius; v < windowRadius; v++)
        {
          const int i = u + windowRadius;
          const int j = v + windowRadius;
          float theta = atan2 (v+0.5, u+0.5);
          if (theta < 0)
          {
            theta += 2 * M_PI;
          }
          assert (theta >= 0 && theta < 2 * M_PI);
          const int slice = (int) floor (theta / (M_PI / norient));
          assert (slice >= 0 && slice < 2 * norient);
          pie(i,j) = slice;
          if ((u+0.5) * (u+0.5) + (v+0.5) * (v+0.5) <= radius * radius)
          {
            pixelCount++;
          }
        }
      }

      Util::Array1D<float> lHist (textons.numChannels);
      Util::Array1D<float> rHist (textons.numChannels);
      Util::Array2D<float> phist (norient*2, textons.numChannels);

      Util::Message::startBlock(width,"TG computation");
      for (int x = 0; x < width; x++)
      {
        Util::Message::stepBlock();
        for (int y = 0; y < height; y++)
        {
          // compute the histograms for windows over this pixel at
          // each orientation
          // compute the pie slice histograms
          phist.init (0);
          for (int u = -windowRadius; u < windowRadius; u++)
          {
            for (int v = -windowRadius; v < windowRadius; v++)
            {
              if ((u+0.5) * (u+0.5) + (v+0.5) * (v+0.5) > radius * radius)
              {
                continue;
              }
              const int xi = x + u;
              const int yi = y + v;
              if (xi < 0 || xi >= width)
              {
                continue;
              }
              if (yi < 0 || yi >= height)
              {
                continue;
              }
              // get the texton
              const int texton = textons.map(xi,yi);
              assert (texton >= 0 && texton < textons.numChannels);

              // compute which pie slice we're in
              const int i = u + windowRadius;
              const int j = v + windowRadius;
              const int slice = pie (i, j);

              // increment histogram
              phist (slice, texton) += 2.0 / pixelCount;
            }
          }

          // set up L histogram (minus last slice)
          lHist.init (0);
          for (int orient = 0; orient < norient - 1; orient++)
          {
            for (int texton = 0; texton < textons.numChannels; texton++)
            {
              lHist (texton) += phist (orient, texton);
            }
          }

          // set up R histogram (minus last slice)
          rHist.init (0);
          for (int orient = norient; orient < 2 * norient - 2; orient++)
          {
            for (int texton = 0; texton < textons.numChannels; texton++)
            {
              rHist (texton) += phist (orient, texton);
            }
          }

          // spin the disc
          for (int orient = 0; orient < norient; orient++)
          {
            // add next slice into L histogram
            int slice = orient + norient - 1;
            for (int texton = 0; texton < textons.numChannels; texton++)
            {
                lHist (texton) += phist (slice, texton);
            }
            // add next slice into R histogram
            slice = Util::wrap(orient - 1, 2 * norient);
            for (int texton = 0; texton < textons.numChannels; texton++)
            {
                rHist (texton) += phist (slice, texton);
            }
            // compute chi-squared distance between L and R histograms
            float dist = 0;
            for (int i = 0; i < textons.numChannels; i++)
            {
              const float den = lHist (i) + rHist (i);
              if (den != 0)
              {
                const float num = lHist (i) - rHist (i);
                dist += num*num / den;
              }
            }
            dist *= 0.5;
            dist = Util::max (0.0f, Util::min (1.0f, dist));
            assert (dist >= 0 && dist <= 1);
            tg(orient,x,y) = dist;

            // subtract first slice from L histogram
            slice = orient;
            for (int texton = 0; texton < textons.numChannels; texton++)
            {
                lHist (texton) -= phist (slice, texton);
            }

            // subtract first slice from R histogram
            slice = orient + norient;
            for (int texton = 0; texton < textons.numChannels; texton++)
            {
                rHist (texton) -= phist (slice, texton);
            }
          }
        }
      }
      Util::Message::endBlock();

      // verify TG values
      for (int orient = 0; orient < norient; orient++)
      {
        for (int x = 0; x < width; x++)
        {
          for (int y = 0; y < height; y++)
          {
            assert (tg(orient,x,y) >= 0);
            assert (tg(orient,x,y) <= 1);
          }
        }
      }
  }


  // given the texton channel membership image and a scale, compute the texton
  // histogram associated with each point 
  void computeTextonHistograms (const TextonMap& textons, const float scale, Histogram& histogram)
  {
      int width = textons.map.size(0);
      int height = textons.map.size(1);

      histogram.resize(width,height,textons.numChannels);
      histogram.init(0);

      int windowRadius = (int)ceil(scale);
       
      // count the number of pixels in a circular window of size scale
      int pixelCount = 0;
      for (int v = -windowRadius; v <= windowRadius; v++)
      {
        for (int u = -windowRadius; u <= windowRadius; u++)
        {
          // make sure we are inside the circular window 
          if (u * u + v * v > scale * scale)
          {
            continue;
          }
          pixelCount++;
        }
      }

      if (pixelCount > 0)
      {
        //for each pixel
        for (int x = 0; x < width; x++)
        {
          for (int y = 0; y < height; y++)
          {
            //look at the histogram over the window    
            for (int v = -windowRadius; v <= windowRadius; v++)
            {
              for (int u = -windowRadius; u <= windowRadius; u++)
              {
                // make sure we are inside the circular window 
                if (u * u + v * v > scale * scale)
                {
                  continue;
                }

                // compute coordinates of point in the window
                int ix = x + u;
                int iy = y + v;

                // if point is out of bounds, skip it
                if (ix < 0 || ix >= width)
                {
                  continue;
                }
                if (iy < 0 || iy >= height)
                {
                  continue;
                }

                assert (ix >= 0);
                assert (ix < width);
                assert (iy >= 0);
                assert (iy < height);

                int texton = textons.map(ix,iy);
                histogram(x,y,texton) += (1.0 / pixelCount);
              }
            }
          }
        }
      }
  }


  //compute texton histogram at given scale, weighted
  //by values in support map.
  void computeTextonHistograms(const TextonMap& textons, const SupportMap& supportMap, Histogram& histogram)
  {
      int width = textons.map.size(0);
      int height = textons.map.size(1);
      histogram.resize(width,height,textons.numChannels);
      histogram.init(0);

      //for each pixel
      Util::Message::startBlock(width,"weighted texton histogram patch computation");
      for (int x = 0; x < width; x++)
      {
        Util::Message::stepBlock();
        for (int y = 0; y < height; y++)
        {
          //look at the histogram over the window    
          for (int i = 0; i < supportMap(x,y).size(); i++)
          {
              // compute coordinates of point in the window
              const int ix = supportMap(x,y)(i).x;
              const int iy = supportMap(x,y)(i).y;
              assert(ix >= 0);
              assert(ix < width);
              assert(iy >= 0);
              assert(iy < height);

              int texton = textons.map(ix,iy);
              histogram(x,y,texton) += supportMap(x,y)(i).sim;
          }
        }
      }
      Util::Message::endBlock();

      Util::Message::startBlock("histogram normalization");
      normalizeHistogram(histogram);
      Util::Message::endBlock();
  }

  // given two points and the stack of texture histograms, this 
  // method computes the chisquared distance between them
  void textureSimilarity (const Histogram& histogram,
                            const int x1, const int y1,
                            const int x2, const int y2, float& similarity)
  {
    chiSquared(histogram,x1,y1,x2,y2,similarity);
  }

} //namespace Group
