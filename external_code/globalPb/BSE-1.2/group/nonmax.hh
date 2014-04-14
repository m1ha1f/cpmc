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

#ifndef BOUNDARY_H 
#define BOUNDARY_H 

#include "image.hh"

//
// non-maximal suppression and subpixel fitting code
//

// TODO: fix parabolaset class to be more reasonable

namespace Group
{
  //
  // structure holding a set of subpixel edgels
  // orientations are in [0,PI).
  //
  class ParabolaSet
  {
    public:
      ParabolaSet () {};
      ParabolaSet (const int width, const int height);
      ParabolaSet (const ParabolaSet& that);
      ~ParabolaSet () {};

      //resize array of parabolas
      void resize(const int width, const int height);

      //zero out all parameters
      void init(); 

      // interpolated maximum orientation energy, valid for all pixels
      Util::Image maxo;
      Util::Image maxe;

      // mask: parabolas are fit only to pixels with edgemap==1
      Util::Image edgemap;

      // parabola parameters; valid only when edgemap==1
      Util::Image X;    
      Util::Image Y;
      Util::Image orientation;
      Util::Image energy;
      Util::Image fiterror;
      Util::Image curvature;
      Util::Image dist;
      Util::Image smooth;
  };  


  //
  // store boundary information on the lattice dual to 
  // the pixel grid
  //
  struct DualLattice
  {
    Util::Image H;
    Util::Image V;
    int width;
    int height;
  };

  //
  // fit cylindrical parabola to a circular patch of radius r centered
  // at (x,y) at angle theta.  return parabola coefficients and 
  // the normalized fit error.  z = ax^2 + bx + c
  //
  void fitParabolaLeastSqr (const int x, const int y, const float theta, const float r,
                            const Util::Image& z, float& a, float& b, float& c, float& error);

  //
  // simple nonmax suppression based on linear interpolation
  // for a 3x3 neighborhood.  single orientation
  //
  void nonmaxSuppress(const Util::Image& mag, const float theta, Util::Image& suppressed);

  //
  // simple nonmax suppression based on linear interpolation
  // for a 3x3 neighborhood. orientation per pixel
  //
  void nonmaxSuppress(const Util::Image& mag, const Util::Image& theta, Util::Image& suppressed);



  //
  // interpolate the orientation matrix and energies at a single orientation.
  // fit clyndrical parabolas to each and return the resulting subpixel info
  // r = radius of square to which we fit the cylindrical parabola
  // theta = single orientation
  //
  void fitSubpixelParabolas (const Util::Image& mag, const float r, const float theta, ParabolaSet& parabs);

  //
  // fit clyndrical parabolas at given max orientation and return the 
  // resulting subpixel info 
  // r = radius of square to which we fit the cylindrical parabola
  // theta = per pixel oreintation
  //
  void fitSubpixelParabolas (const Util::Image& mag, const Util::Image& theta, const float r,
                              ParabolaSet& parabs);

  //
  // interpolate the maximal orientation.
  // fit clyndrical parabolas to each and return the resulting subpixel info
  // r = radius of square to which we fit the cylindrical parabola
  //
  void fitSubpixelParabolas (const Util::ImageStack& mag, const float r, ParabolaSet& parabs);


  //
  // given the subpixel parabolas and the orientation fits, this function
  // projects pixels down onto the dual lattice by tracing the appropriate
  // digital line
  //
  void intersectPixels (const ParabolaSet& parabolas, DualLattice& L);

  //
  // given the subpixel parabolas and the orientation fits, this function
  // intersects edges with the pixels: one pixel per parabola 
  //
  void intersectPixels (const ParabolaSet& parabolas, Util::Image& pb);
     
  //
  // estimate derivatives using LOWESS style smoothing
  //
  void lowess(const float theta, const float r, const Util::Image& a,
              Util::Image& d0, Util::Image& d1, Util::Image& d2);

  //
  // estimate smoothed 0th,1st,2nd derivatives using gaussian derivatives convolution
  // sigma = 3/r;
  //
  void derivs(const float theta, const float r, const Util::Image& a,
              Util::Image& d0, Util::Image& d1, Util::Image& d2);
  

} //namespace Group
#endif
