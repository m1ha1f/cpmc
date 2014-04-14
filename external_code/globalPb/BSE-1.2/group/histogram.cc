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


#include <assert.h>
#include "histogram.hh"
#include "array.hh"

namespace Group
{

  //
  // L1 normalized a histogram so that the sum of bins at
  // each location is 1.
  //
  void normalizeHistogram(Histogram& histogram)
  {
    const int width = histogram.size(0);
    const int height = histogram.size(1);
    const int nbins = histogram.size(2);

    for (int x = 0; x < width; x++)
    {
      for (int y = 0; y < height; y++)
      {
        // L1-normalize the histogram at (x,y)
        float histNorm = 0;
        for (int bin = 0; bin < nbins; bin++)
        {
          assert(histogram(x,y,bin) >= 0.0f);
          histNorm += histogram(x,y,bin);
        }

        if (histNorm > 0.0f)
        {
          for (int bin = 0; bin < nbins; bin++)
          {
            histogram(x,y,bin) /= histNorm;
          }
        }
        else
        {
          for (int bin = 0; bin < nbins; bin++)
          {
            histogram(x,y,bin) = 0;
          }
        }
      } 
    }
  }


  //
  // given two points and the color histograms, compute
  // chisquared distance between them.  the histogram
  // should already be normalized.
  //
  void chiSquared(const Histogram& histogram,
                  const int x1, const int y1,
                  const int x2, const int y2, float& similarity)
  {
    assert(x1 >= 0 && x1 < histogram.size(0));
    assert(y1 >= 0 && y1 < histogram.size(1));
    assert(x2 >= 0 && x2 < histogram.size(0));
    assert(y2 >= 0 && y2 < histogram.size(1));

    const int nbins = histogram.size(2);
    similarity = 0.0f;
    for (int bin = 0; bin < nbins; bin++)
    {
      float h1 = histogram(x1,y1,bin);
      float h2 = histogram(x2,y2,bin);
      if ((h1 + h2) != 0)
      {
          similarity += ((h1 - h2) * (h1 - h2)) / (2 * (h1 + h2));
      }
    }
    if (similarity > 1)  //something fp ??
    {
      similarity = 1.0f;
    }
    assert(similarity >= 0.0f && similarity <= 1.0f); 
  }

} //namespace Group


