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
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <math.h>

#include "smatrix.hh"
#include "region.hh"
#include "ic.hh"
#include "affinity.hh"

namespace Group
{

  //////////////////////////////////////////////////////////////////////////////
  //
  // probability-of-same-segment logistic fits for color and grayscale
  // 
  // optimized for dthresh = 20
  //
  //

  float C_P_SS(float bp, float tp, float cp)
  {
    float val = 2.0404;

    val += (bp / 0.1989)*-0.4301;
    val += (tp / 0.1639)*-1.0353;
    val += (cp / 0.3554)*-0.5053;

    const float post = 1.0 / (1.0 + exp(-val));
    if (!finite(post)) { return 0; }
    if (post < 0) { return 0; }
    if (post > 1) { return 1; }
    return post;
  }

  float C_IC_SS(float ic)
  {
    float val = 1.8682;

    val += (ic / 0.3130)*-1.3113;

    const float post = 1.0 / (1.0 + exp(-val));
    if (!finite(post)) { return 0; }
    if (post < 0) { return 0; }
    if (post > 1) { return 1; }
    return post;
    
  }

  float C_B_SS(float ic, float bp, float tp, float cp)
  {
    float val = 2.1585;

    val += (ic / 0.3130)*-0.6615;
    val += (bp / 0.1989)*-0.3093;
    val += (tp / 0.1639)*-0.6605;
    val += (cp / 0.3554)*-0.2920;

    const float post = 1.0 / (1.0 + exp(-val));
    if (!finite(post)) { return 0; }
    if (post < 0) { return 0; }
    if (post > 1) { return 1; }
    return post;
  }

  //////////////////////////////////////////////////////////////////////////////

  float G_P_SS(float bp, float tp)
  {
    float val = 2.1339;

    val += (bp / 0.1932)*-0.6303;
    val += (tp / 0.1348)*-1.1480;

    const float post = 1.0 / (1.0 + exp(-val));
    if (!finite(post)) { return 0; }
    if (post < 0) { return 0; }
    if (post > 1) { return 1; }
    return post;
    
  }

  float G_IC_SS(float ic)
  {
    float val = 1.9755;

    val += (ic / 0.2841)*-1.2532;

    const float post = 1.0 / (1.0 + exp(-val));
    if (!finite(post)) { return 0; }
    if (post < 0) { return 0; }
    if (post > 1) { return 1; }
    return post;
  }

  float G_B_SS(float ic, float bp, float tp)
  {
    float val = 2.3738;

    val += (ic / 0.2841)*-0.6735;
    val += (bp / 0.1932)*-0.3661;
    val += (tp / 0.1348)*-0.8013;

    const float post = 1.0 / (1.0 + exp(-val));
    if (!finite(post)) { return 0; }
    if (post < 0) { return 0; }
    if (post > 1) { return 1; }
    return post;
  }

  //////////////////////////////////////////////////////////////////////////////
  
  //
  // compute similarities for the entire image in scanline order
  //
  void computeAffinities(const SupportMap& icmap, const Region& region, const cueModeType cueMode,
                           const float sigma, const float dthresh, const bool useColor,
                            SMatrix** affinities)
  {
    Util::Array2D<bool> full(icmap.size(0),icmap.size(1));
    full.init(true);
    computeAffinities(icmap,region,cueMode,sigma,dthresh,useColor,full,affinities);
  }

  //
  // compute similarities for the set of "true" pixels in region.  
  // affinity matrix is ordered in scanline order 
  //
  void computeAffinities(const SupportMap& icmap, const Region& region, const cueModeType cueMode,
                           const float sigma, const float dthresh, const bool useColor,
                           const Util::Array2D<bool> mask, SMatrix** affinities)
  {
    int width = icmap.size(0);
    int height = icmap.size(1);
    assert(icmap.size(0) == mask.size(0));
    assert(icmap.size(1) == mask.size(1));

    //make sure region object has color info before we request it
    if ((cueMode == patch) || (cueMode == both))
    {
      if (useColor)
      {
        assert(useColor == region.useColor());
      }
    }

    //build a scanline order index 
    int numPixels = 0;
    Util::Array2D<int> index(width,height);
    index.init(-1);
    for (int y = 0; y < height; y++)
    {
      for (int x = 0; x < width; x++)
      {
        if (mask(x,y))
        {
          index(x,y) = numPixels;
          numPixels++;
        }
      }
    }

    //sparse matrix data
    int* nz = new int[numPixels];           //number of non-zero entries in each row
    double** vals = new double*[numPixels];   //the values in each row
    int** col = new int*[numPixels];        //the column number for each value in the row  
   
    int dthreshi = (int)ceil(dthresh);
    int wd = (2*dthreshi+1);                 //window diameter 
    Util::Array1D<PointIC> connections(wd*wd);
    
    int row = 0;
    for (int x = 0; x < width; x++)
    {
        for (int y = 0; y < height; y++)
        {
            if (mask(x,y))
            {
                row = index(x,y);  //the row we are working on
                nz[row] = 0;       //connection count for row i
                int icIndex = 0;        //index into sparse supportMap
                for (int u = -dthreshi; u <= dthreshi; u++)
                {
                  int yy = y + u;
                  for (int v = -dthreshi; v <= dthreshi; v++)
                  {
                    int xx = x + v;
                    if (xx < 0 || xx >= width) {continue;}
                    if (yy < 0 || yy >= height) {continue;}
                    if (!mask(xx,yy)) {continue;}
                    if (u*u+v*v > dthresh*dthresh) {continue;}
                  
                    //increment our index into the support map
                    if ((cueMode == ic) || (cueMode == both))
                    {
                      while( icIndex < icmap(x,y).size() && 
                             icmap(x,y)(icIndex).y < yy) 
                      {
                        icIndex++;
                      }
                      while( icIndex < icmap(x,y).size() && 
                              icmap(x,y)(icIndex).x < xx) 
                      {
                        icIndex++;
                      }
                    }


                    float pss = 0.0;     //connection strength

                    if ((u == 0) && (v == 0))
                    {
                      pss = 1.0f;
                    } 
                    else
                    {
                      switch(cueMode)
                      {
                        case patch:
                        {
                          const float bp = region.getBSim(0,x,y,xx,yy);
                          const float tp = region.getTSim(0,x,y,xx,yy);
                          if (useColor)
                          {
                            float cp = region.getCSim(0,x,y,xx,yy);
                            pss = C_P_SS(bp,tp,cp);
                          } 
                          else 
                          { 
                            pss = G_P_SS(bp,tp); 
                          }
                        }
                        break;

                        case ic:
                        {
                          float icsim = 0.0f;
                          if (icIndex < icmap(x,y).size() &&
                               icmap(x,y)(icIndex).x == xx &&
                                icmap(x,y)(icIndex).y == yy)
                          {
                            icsim = icmap(x,y)(icIndex).sim;
                            icIndex++;
                          }
                          if (useColor)
                          {
                            pss = C_IC_SS(1-icsim);
                          }
                          else 
                          { 
                            pss = G_IC_SS(1-icsim);
                          }
                        }
                        break;

                        case both:
                        {
                          float icsim = 0.0f;
                          if (icIndex < icmap(x,y).size() &&
                               icmap(x,y)(icIndex).x == xx &&
                                icmap(x,y)(icIndex).y == yy)
                          {
                            icsim = icmap(x,y)(icIndex).sim;
                            icIndex++;
                          }
                          const float bp = region.getBSim(0,x,y,xx,yy);
                          const float tp = region.getTSim(0,x,y,xx,yy);
                          if (useColor)
                          {
                            float cp = region.getCSim(0,x,y,xx,yy);
                            pss = C_B_SS(1-icsim,bp,tp,cp); 
                          } 
                          else 
                          { 
                            pss = G_B_SS(1-icsim,bp,tp); 
                          }
                        }
                        break;
                      }//switch(cueMode)
                    }//if (u==0) & (v==0)

                    connections(nz[row]).x = xx;
                    connections(nz[row]).y = yy;
                    connections(nz[row]).sim = pss;
                    nz[row]++;
                  }//for v
                }//for u

                //fill in entries of sparse matrix
                vals[row] = new double[nz[row]];
                col[row] = new int[nz[row]];
                for (int j = 0; j < nz[row]; j++)
                {
                  float val = exp( -(1-connections(j).sim) / sigma);
                  assert((val >= 0.0) && (val <= 1.0));
                  vals[row][j] = val;
                  col[row][j] = index(connections(j).x,connections(j).y);
                }
            }//if mask(x,y)
        }//for y
    }//for x

    *affinities = new SMatrix(numPixels,nz,col,vals);
    (*affinities)->symmetrize();
  }

} //namespace Group



