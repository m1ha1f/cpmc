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
#include <limits.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <ctype.h>
#include <memory.h>
#include <time.h>

#include "filterbank.hh"
#include "nonmax.hh"
#include "util.hh"
#include "image.hh"

namespace Group
{

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // data structure for a set of subpixel edgels
  //

  ParabolaSet::ParabolaSet (const int width, const int height)
  {
    resize(width,height);
  }

  ParabolaSet::ParabolaSet (const ParabolaSet& that)
  {
      maxo = that.maxo;
      maxe = that.maxe;
      edgemap = that.edgemap;
      X = that.X;
      Y = that.Y;
      orientation = that.orientation;
      energy = that.energy;
      fiterror = that.fiterror;
      curvature = that.curvature;
      dist = that.dist;
      smooth = that.smooth;
  }

  void ParabolaSet::resize(const int width, const int height)
  {
      maxo.resize(width,height);
      maxe.resize(width,height);
      edgemap.resize(width,height);
      X.resize(width,height);
      Y.resize(width,height);
      orientation.resize(width,height);
      energy.resize(width,height);
      fiterror.resize(width,height);
      curvature.resize(width,height);
      dist.resize(width,height);
      smooth.resize(width,height);
  }

  void ParabolaSet::init()
  {
      maxo.init(0);
      maxe.init(0);
      edgemap.init(0);
      X.init(0);
      Y.init(0);
      orientation.init(0);
      energy.init(0);
      fiterror.init(0);
      curvature.init(0);
      dist.init(0);
      smooth.init(0);
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // fit cylindrical parabola to a circular patch of radius r centered
  // at (x,y) at angle theta.  return parabola coefficients and
  // the normalized fit error.  z = ax^2 + bx + c
  //
  void fitParabolaLeastSqr (
      const int x, const int y,
      const float theta, const float r,
      const Util::Image& z,
      float& a, float& b, float& c, float& error)
  {
      a = b = c = 0;
      error = 0;

      const int width = z.size(0);
      const int height = z.size(1);
      const int rad = (int) ceil (r);
      const float sint = sin(theta);
      const float cost = cos(theta);

      // d := projection of the 3x3 grid onto line perpendicular to theta
      // di := Sum over u,v of d[u,v]^i
      float d0 = 0, d1 = 0, d2 = 0, d3 = 0, d4 = 0;
      float v0 = 0, v1 = 0, v2 = 0;
      for (int v = -rad; v <= rad; v++) {
          for (int u = -rad; u <= rad; u++) {
              const int xi = x + u;
              const int yi = y + v;
              if (xi < 0 || xi >= width) { continue; }
              if (yi < 0 || yi >= height) { continue; }
              if (u*u + v*v > r*r) { continue; }
              const float zi = z(xi,yi);
              const float di = -u*sint + v*cost; // distance to diameter
              d0 += 1;
              d1 += di;
              d2 += di * di;
              d3 += di * di * di;
              d4 += di * di * di * di;
              v0 += zi;
              v1 += di * zi;
              v2 += di * di * zi;
          }
      }

      // solve the linear system A [c b a] = v
      // A[3][3] = {{d0, d1, d2},
      //            {d1, d2, d3},
      //            {d2, d3, d4}}
      const float detA =
          -d2*d2*d2 + 2*d1*d2*d3 - d0*d3*d3 - d1*d1*d4 + d0*d2*d4;
      const float invA[3][3] = {
          {-d3*d3 + d2*d4,  d2*d3 - d1*d4, -d2*d2 + d1*d3},
          { d2*d3 - d1*d4, -d2*d2 + d0*d4,  d1*d2 - d0*d3},
          {-d2*d2 + d1*d3,  d1*d2 - d0*d3, -d1*d1 + d0*d2}
      };

      if (detA != 0)
      {
        c = (invA[0][0] * v0 + invA[0][1] * v1 + invA[0][2] * v2) / detA;
        b = (invA[1][0] * v0 + invA[1][1] * v1 + invA[1][2] * v2) / detA;
        a = (invA[2][0] * v0 + invA[2][1] * v1 + invA[2][2] * v2) / detA;
      }
      else
      {
        //parabola fit failed due to invertability issue.
        //this can happen for very small fit radius when
        //we are clipped by an image boundary.  set some
        //default values so that we never choose this as
        //a boundary location
        c = 0.0f;
        b = 10.0f;
        a = 1.0f;
      }

      //compute residual error
      float znorm = 0;
      float dznorm = 0;
      for (int v = -rad; v <= rad; v++) {
          for (int u = -rad; u <= rad; u++) {
              const int xi = x + u;
              const int yi = y + v;
              if (xi < 0 || xi >= width) { continue; }
              if (yi < 0 || yi >= height) { continue; }
              if (u*u + v*v > r*r) { continue; }
              const float zi = z(xi,yi);
              const float di = -u*sint + v*cost;
              znorm += zi * zi;
              const float dzi = a*di*di + b*di + c - zi;
              dznorm += dzi * dzi;
          }
      }
      if (znorm != 0) {
          error = sqrt (dznorm / znorm); // normalized RMS error
      }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // simple nonmax suppression based on linear interpolation
  // for a 3x3 neighborhood.
  //
  // single orientation 
  //
  void nonmaxSuppress(const Util::Image& mag, const float theta, Util::Image& suppressed)
  {
    const int width = mag.size(0);
    const int height = mag.size(1);
    suppressed.resize(width,height);

    float tval = Util::mod(theta + (M_PI/2),M_PI);

    for (int y = 1; y < height-1; y++)
    {
      for (int x = 1; x < width-1; x++)
      {
         float mval = mag(x,y);
         suppressed(x,y) = mval;

         if ((tval >= 0) && (tval < (M_PI / 4)))
         {
           float d = tan(tval);
           float v1 = (1-d)*mag(x+1,y) + d*mag(x+1,y+1);
           float v2 = (1-d)*mag(x-1,y) + d*mag(x-1,y-1);
           if ((v2 > mval) || (v1 > mval))
           {
             suppressed(x,y) = 0;
           }
         }
         else if ((tval >= (M_PI/4)) && (tval < (M_PI / 2)))
         {
           float d = tan((M_PI/2) - tval);
           float v1 = (1-d)*mag(x,y+1) + d*mag(x+1,y+1);
           float v2 = (1-d)*mag(x,y-1) + d*mag(x-1,y-1);
           if ((v2 > mval) || (v1 > mval))
           {
             suppressed(x,y) = 0;
           }
         }
         else if ((tval >= (M_PI/2)) && (tval < (3*M_PI / 4)))
         {
           float d = tan(tval - (M_PI/2));
           float v1 = (1-d)*mag(x,y+1) + d*mag(x-1,y+1);
           float v2 = (1-d)*mag(x,y-1) + d*mag(x+1,y-1);
           if ((v2 > mval) || (v1 > mval))
           {
             suppressed(x,y) = 0;
           }
         }
         else if ((tval >= (3*M_PI/4)) && (tval < M_PI))
         {
           float d = tan(M_PI - tval);
           float v1 = (1-d)*mag(x-1,y) + d*mag(x-1,y+1);
           float v2 = (1-d)*mag(x+1,y) + d*mag(x+1,y-1);
           if ((v2 > mval) || (v1 > mval))
           {
             suppressed(x,y) = 0;
           }
         }
      }
    }
  }

  //
  // simple nonmax suppression based on linear interpolation
  // for a 3x3 neighborhood.
  //
  // oreintation per pixel
  //
  void nonmaxSuppress(const Util::Image& mag, const Util::Image& theta, Util::Image& suppressed)
  {
    const int width = mag.size(0);
    const int height = mag.size(1);
    suppressed.resize(width,height);

    for (int y = 1; y < height-1; y++)
    {
      for (int x = 1; x < width-1; x++)
      {
         float tval = Util::mod(theta(x,y) + (M_PI/2),M_PI);
         float mval = mag(x,y);
         suppressed(x,y) = mval;

         if ((tval >= 0) && (tval < (M_PI / 4)))
         {
           float d = tan(tval);
           float v1 = (1-d)*mag(x+1,y) + d*mag(x+1,y+1);
           float v2 = (1-d)*mag(x-1,y) + d*mag(x-1,y-1);
           if ((v2 > mval) || (v1 > mval))
           {
             suppressed(x,y) = 0;
           }
         }
         else if ((tval >= (M_PI/4)) && (tval < (M_PI / 2)))
         {
           float d = tan((M_PI/2) - tval);
           float v1 = (1-d)*mag(x,y+1) + d*mag(x+1,y+1);
           float v2 = (1-d)*mag(x,y-1) + d*mag(x-1,y-1);
           if ((v2 > mval) || (v1 > mval))
           {
             suppressed(x,y) = 0;
           }
         }
         else if ((tval >= (M_PI/2)) && (tval < (3*M_PI / 4)))
         {
           float d = tan(tval - (M_PI/2));
           float v1 = (1-d)*mag(x,y+1) + d*mag(x-1,y+1);
           float v2 = (1-d)*mag(x,y-1) + d*mag(x+1,y-1);
           if ((v2 > mval) || (v1 > mval))
           {
             suppressed(x,y) = 0;
           }
         }
         else if ((tval >= (3*M_PI/4)) && (tval < M_PI))
         {
           float d = tan(M_PI - tval);
           float v1 = (1-d)*mag(x-1,y) + d*mag(x-1,y+1);
           float v2 = (1-d)*mag(x+1,y) + d*mag(x+1,y-1);
           if ((v2 > mval) || (v1 > mval))
           {
             suppressed(x,y) = 0;
           }
         }
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // interpolate max using 3 point technique
  //
  static void
  interpolateMaxLagrange (const float x1, const float y1,
                          const float x2, const float y2,
                          const float x3, const float y3,
                          float &x, float &y)
  {
      assert (x1 < x2);
      assert (x2 < x3);

      assert (y2 >= y1);
      assert (y2 >= y3);

      // Location of max.
      x = (x3 * x3 * (y1 - y2) + x1 * x1 * (y2 - y3) + x2 * x2 * (y3 - y1))
          / (2. * (x3 * (y1 - y2) + x1 * (y2 - y3) + x2 * (y3 - y1)));

      // Lagrange interpolation to find max value.
      y = y1 * ((x - x2) * (x - x3)) / ((x1 - x2) * (x1 - x3))
          + y2 * ((x - x1) * (x - x3)) / ((x2 - x1) * (x2 - x3))
          + y3 * ((x - x1) * (x - x2)) / ((x3 - x1) * (x3 - x2));

      if (!finite (x) || !finite (y))
        {
            x = x2;
            y = y2;
        }
  }


  //
  // fit clyndrical parabolas at given max orientation and return the
  // resulting subpixel info
  // r = radius of disc on which we fit the cylindrical parabola
  // theta = orientation of fits
  //
  void fitSubpixelParabolas (const Util::Image& mag, const float theta, 
                             float r, ParabolaSet& parabs)
  {
      const int width = mag.size(0);
      const int height = mag.size(1);
      parabs.resize(width,height);
      parabs.init();

      r = Util::max(r,1.5f); 
      //compute the best paraboloid fit at each point (x,y)
      for (int x = 0; x < width; x++)
      {
        for (int y = 0; y < height; y++)
        {
              float maxe = mag(x,y);
              assert (theta >= 0 && theta < M_PI);

              // save the maximal orientation and energy at all pixels
              parabs.maxo(x,y) = theta;
              parabs.maxe(x,y) = maxe;

              // these are the least squares paraboloid coeff. estimates
              // for z = a*d^2 + b*d + c, where d = -x*sint + y*cost
              // error is the normalized fit error
              float a = 0, b = 0, c = 0;
              float error = 0;

              fitParabolaLeastSqr ( x, y, theta, r, mag, a, b, c, error);

              // compute various parameters of fitted parabola
              // curvature is normalized to the window radius
              const float curvature = a * r * r;
              const float dist = -b / (2*a);
              const float energy = a*dist*dist + b*dist + c;
              const float smooth = Util::max(c,0.0f);
              const float ycoord = dist*cos(theta);
              const float xcoord = -dist*sin(theta);

              // save fit info for all putative needles
              parabs.X(x,y) = xcoord;
              parabs.Y(x,y) = ycoord;
              parabs.energy(x,y) = energy;
              parabs.orientation(x,y) = theta;
              parabs.curvature(x,y) = curvature;
              parabs.fiterror(x,y) = error;
              parabs.dist(x,y) = dist;
              parabs.smooth(x,y) = smooth;

              // if any of the values are infinite or NaN, then skip
              if (!finite(xcoord)) { continue; }
              if (!finite(ycoord)) { continue; }
              if (!finite(energy)) { continue; }
              if (!finite(theta)) { continue; }
              if (!finite(curvature)) { continue; }
              if (!finite(error)) { continue; }

              // if curvature is positive, then this is not a max; skip
              if (curvature > -1e-10) { continue; }

              // if the center is out of this pixel's neighborhood, then skip
              if (fabs(dist) > 1) { continue; }

              // only mark valid those needles that are maxima in the
              // pixel's immediate neighborhood
              parabs.edgemap(x,y) = 1;
          }
      }
  }


  //
  // fit clyndrical parabolas at given max orientation and return the
  // resulting subpixel info
  // r = radius of square to which we fit the cylindrical parabola
  //
  void fitSubpixelParabolas (const Util::Image& mag, const Util::Image& theta, 
                             float r, ParabolaSet& parabs)
  {
      const int width = mag.size(0);
      const int height = mag.size(1);
      parabs.resize(width,height);
      parabs.init();

      r = Util::max(r,1.5f); 

      //compute the best paraboloid fit at each point (x,y)
      for (int x = 0; x < width; x++)
      {
        for (int y = 0; y < height; y++)
        {
              float eta = theta(x,y); 
              float maxe = mag(x,y);
              assert (eta >= 0 && eta < M_PI);

              // save the maximal orientation and energy at all pixels
              parabs.maxo(x,y) = eta;
              parabs.maxe(x,y) = maxe;

              // these are the least squares paraboloid coeff. estimates
              // for z = a*d^2 + b*d + c, where d = -x*sint + y*cost
              // error is the normalized fit error
              float a = 0, b = 0, c = 0;
              float error = 0;

              fitParabolaLeastSqr ( x, y, eta, r, mag, a, b, c, error);

              // compute various parameters of fitted parabola
              // curvature is normalized to the window radius
              const float curvature = a * r * r;
              const float dist = -b / (2*a);
              const float energy = a*dist*dist + b*dist + c;
              const float smooth = Util::max(c,0.0f);
              const float ycoord = dist*cos(eta);
              const float xcoord = -dist*sin(eta);

              // save fit info for all putative needles
              parabs.X(x,y) = xcoord;
              parabs.Y(x,y) = ycoord;
              parabs.energy(x,y) = energy;
              parabs.orientation(x,y) = eta;
              parabs.curvature(x,y) = curvature;
              parabs.fiterror(x,y) = error;
              parabs.dist(x,y) = dist;
              parabs.smooth(x,y) = smooth;

              // if any of the values are infinite or NaN, then skip
              if (!finite(xcoord)) { continue; }
              if (!finite(ycoord)) { continue; }
              if (!finite(energy)) { continue; }
              if (!finite(eta)) { continue; }
              if (!finite(curvature)) { continue; }
              if (!finite(error)) { continue; }

              // if curvature is positive, then this is not a max; skip
              if (curvature > -1e-10) { continue; }

              // if the center is out of this pixel's neighborhood, then skip
              if (fabs(dist) > 1) { continue; }

              // only mark valid those needles that are maxima in the
              // pixel's immediate neighborhood
              parabs.edgemap(x,y) = 1;
          }
      }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // interpolate the orientation matrix and energies,
  // fit clyndrical parabolas to each and return the resulting subpixel info
  // r = radius of square to which we fit the cylindrical parabola
  //
  void fitSubpixelParabolas (const Util::ImageStack& mag, float r, ParabolaSet& parabs)
  {
    r = Util::max(r,1.5f); 
    const int norient = mag.size(0);
    const int width = mag.size(1);
    const int height = mag.size(2);

    parabs.resize(width,height);
    parabs.init();

    //compute the best paraboloid fit at each point (x,y)
    for (int x = 0; x < width; x++)
    {
      for (int y = 0; y < height; y++)
      {
        // compute the max energy orientation at this point
        int maxOrient = 0;
        float maxEnergy = mag(0,x,y);
        for (int orient = 0; orient < norient; orient++) {
            const float e = mag(orient,x,y);
            if (e > maxEnergy) {
                maxOrient = orient;
                maxEnergy = e;
            }
        }
        const float oStep = M_PI / norient;
        float theta = maxOrient * oStep;
        float maxe = mag(maxOrient,x,y);
        assert (theta >= 0 && theta < M_PI);

        if (norient > 2) {
          // interpolate around the max energy point to get a more
          // accurate theta; we will use this orientation for the
          // parabola fit
          const float x1 = theta - oStep;
          const float x2 = theta;
          const float x3 = theta + oStep;
          const float y1 = mag(Util::wrap(maxOrient-1,norient),x,y);
          const float y2 = mag(maxOrient,x,y);
          const float y3 = mag(Util::wrap(maxOrient+1,norient),x,y);
          float maxo,maxeFit;
          interpolateMaxLagrange (x1, y1, x2, y2, x3, y3, maxo, maxeFit);
          assert(finite(maxo));
          assert(finite(maxeFit));
          if (fabs(theta-maxo) < oStep) 
          {
            theta = maxo;
            maxe = maxeFit;
            if (theta < 0) { theta += M_PI; }
            if (theta >= M_PI) { theta -= M_PI; }
            if (theta < 0) { theta = 0; }
            if (theta >= M_PI) { theta = 0; }
            assert (theta >= 0 && theta < M_PI);
          } // else bad numerical failure!
        }
        assert (theta >= 0 && theta < M_PI);

        // save the interpolated maximal orientation and energy 
        // at all pixels
        parabs.maxo(x,y) = theta;
        parabs.maxe(x,y) = maxe;

        // these are the least squares paraboloid coeff. estimates
        // for z = a*d^2 + b*d + c, where d = -x*sint + y*cost
        // error is the normalized fit error
        float a = 0, b = 0, c = 0;
        float error = 0;

        fitParabolaLeastSqr ( x, y, theta, r, *mag.slice(maxOrient), a, b, c, error);

        // compute various parameters of fitted parabola
        // curvature is normalized to the window radius
        const float curvature = a * r * r;
        const float dist = -b / (2*a);
        const float energy = a*dist*dist + b*dist + c;
        const float smooth = Util::max(c,0.0f);
        const float ycoord = dist*cos(theta);
        const float xcoord = -dist*sin(theta);

        // save fit info for all putative needles
        parabs.X(x,y) = xcoord;
        parabs.Y(x,y) = ycoord;
        parabs.energy(x,y) = energy;
        parabs.orientation(x,y) = theta;
        parabs.curvature(x,y) = curvature;
        parabs.fiterror(x,y) = error;
        parabs.dist(x,y) = dist;
        parabs.smooth(x,y) = smooth;

        // if any of the values are infinite or NaN, then skip
        if (!finite(xcoord)) { continue; }
        if (!finite(ycoord)) { continue; }
        if (!finite(energy)) { continue; }
        if (!finite(theta)) { continue; }
        if (!finite(curvature)) { continue; }
        if (!finite(error)) { continue; }

        // if curvature is positive, then this is not a max; skip
        if (curvature > -1e-10) { continue; }

        // if the center is out of this pixel's neighborhood, then skip
        if (fabs(dist) > 1) { continue; }

        // only mark valid those needles that are maxima in the
        // pixel's immediate neighborhood
        parabs.edgemap(x,y) = 1;
      }
    }
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // project an edgel with energy pb onto the dual lattice by tracing the
  // breshenham line that lies underneath it.
  // 
  void traceBoundary(const float pb, const float ex1, const float ey1, 
                     const float ex2, const float ey2, DualLattice& L)
  {
      const int width = L.H.size(0);
      const int height = L.V.size(1);

      int dx = (int)rintf(ex2) - (int)rintf(ex1);
      int dy = (int)rintf(ey2) - (int)rintf(ey1);
      float x1 = 0;
      float y1 = 0;
      float x2 = 0;
      float y2 = 0;
      if (dx > 0)
      {
        x1 = ex1; 
        y1 = ey1; 
        x2 = ex2; 
        y2 = ey2; 
      }
      else
      {
        x1 = ex2;
        y1 = ey2;
        x2 = ex1;
        y2 = ey1;
        dx = -dx;
        dy = -dy;
      }

      const int steps = Util::max(abs(dx),abs(dy));
      if (steps == 0) { return; }
                                                                                                                 
      const float xincr = (x2-x1) / (float) steps;
      const float yincr = (y2-y1) / (float) steps;
                                                                                                                 
      float x = x1;
      float y = y1;
      int oldx = (int)rint(x);
      int oldy = (int)rint(y);
      for (int k = 0; k < steps; k++)
      {
        float mx = x + 0.5*xincr; 
        float lmx = mx - floor(mx);
        float my = y + 0.5*yincr; 
        float lmy = my - floor(my);

        x += xincr;
        y += yincr;
        const int xi = (int) rintf(x);
        const int yi = (int) rintf(y);

        //don't try to accumulate out of bounds
        if ((oldx < 0) || (oldx >= width)) {continue;}
        if ((xi < 0) || (xi >= width)) {continue;}
        if ((oldy < 0) || (oldy >= height)) {continue;}
        if ((yi < 0) || (yi >= height)) {continue;}

        if ((oldx != xi) && (oldy != yi))
        {
          //diagonal move 
          if (dy > 0)
          {
            if ( (lmx-lmy) > 0 )
            {
              L.H(oldx,oldy) = pb;
              L.V(xi,oldy) = pb;
            }
            else
            {
              L.V(oldx,oldy) = pb;
              L.H(oldx,yi) = pb;
            }
          }
          else
          {
            if ( (lmx+lmy) > 1 )
            {
              L.H(oldx,oldy) = pb;
              L.V(xi,yi) = pb;
            }
            else
            {
              L.V(oldx,yi) = pb;
              L.H(oldx,yi) = pb;
            }
          }
        }
        else if (oldy != yi)
        {
          //vertical move
          if (dy > 0)
          {
            L.V(oldx,oldy) = pb;    
          }
          else
          {
            L.V(oldx,yi) = pb;    
          }
        }
        else if (oldx != xi)
        {
          //horizontal move
          L.H(oldx,oldy) = pb;    
        }
        else
        {
          //TODO: this is some bug.  testcase: 188005.jpg
          fprintf(stderr,"no move!! %d %d %3.3f %3.3f\n",dx,dy,xincr,yincr);
        }
        oldx = xi;
        oldy = yi;
      }

      assert(fabs(x-x2)<0.01);
      assert(fabs(y-y2)<0.01);
  }
                                                                                                                 

  //
  // given a set of parabolas, intersect them all with the dual lattice
  //
  void intersectPixels (const ParabolaSet& parabolas, DualLattice& L)
  {
      const float edgelLength = 1.1f;
      const int width = parabolas.edgemap.size(0);
      const int height = parabolas.edgemap.size(1);
      L.width = width;
      L.height = height;
      L.H.resize(width,height+1);
      L.V.resize(width+1,height);
      L.H.init(0);
      L.V.init(0);

      for (int x = 0; x < width; x++)
      {
        for (int y = 0; y < height; y++)
        {
          if (parabolas.edgemap(x,y) != 1) 
          {
            //nearest edgel is more than 1 pixel away
            //from my center
            continue; 
          }

          // edgel parameters
          //const float dval = parabolas.dist(x,y);
          const float energy = parabolas.energy(x,y);
          const float theta = parabolas.orientation(x,y);
          const float dx = 0.5*edgelLength*cos(theta);
          const float dy = 0.5*edgelLength*sin(theta);
          float x1 = (float)x + parabolas.X(x,y) + dx;
          float x2 = (float)x + parabolas.X(x,y) - dx;
          float y1 = (float)y + parabolas.Y(x,y) + dy;
          float y2 = (float)y + parabolas.Y(x,y) - dy;
          traceBoundary(energy,x1,y1,x2,y2,L);
        }
      }
  }


  //
  // given the subpixel parabolas and the orientation fits, this function
  // intersects edges with the pixels: one pixel per parabola
  //
  void intersectPixels (const ParabolaSet& parabolas, Util::Image& pb)
  {
      const int width = parabolas.edgemap.size(0);
      const int height = parabolas.edgemap.size(1);
      pb.resize(width,height);
      pb.init(0);
      for (int x = 0; x < width; x++)
      {
          for (int y = 0; y < height; y++)
          {
              if (parabolas.edgemap(x,y) != 1) 
              {
                continue; 
              }
              float dval = parabolas.dist(x,y);
              const float dx = parabolas.X(x,y);
              const float dy = parabolas.Y(x,y);
              int xi = (int) rint (x + dx);
              int yi = (int) rint (y + dy);
              if (xi < 0 || xi >= width) 
              {
                continue; 
              }
              if (yi < 0 || yi >= height) 
              { 
                continue; 
              }

              if (dval < (sqrt(2)/2))
              {
                float energy = parabolas.energy(x,y);
                energy = Util::max(energy,0.0f);
                energy = Util::min(energy,1.0f);
                pb(xi,yi) = Util::max(pb(xi,yi),energy);
              }
          }
      }
  }
  

  //
  // estimate smoothed derivatives using LOWESS style smoothing
  //
  void lowess(const float theta, const float r, const Util::Image& signal, 
              Util::Image& d0, Util::Image& d1, Util::Image& d2)

  {
    const int width = signal.size(0);
    const int height = signal.size(1);
    d0.resize(width,height);
    d1.resize(width,height);
    d2.resize(width,height);
    for(int x = 0; x < width; x++)
    {
      for(int y = 0; y < height; y++)
      {
        float a, b, c, error;
        fitParabolaLeastSqr(x, y, theta, r, signal, a, b, c, error);
        d2(x,y) = 2*a;
        d1(x,y) = b;
        d0(x,y) = c;
      }
    }
   
  }


  //
  // estimate smoothed 0th,1st,2nd derivatives using gaussian derivative convolution
  // gaussian sigma = r / 3
  //
  void derivs(const float theta, const float r, const Util::Image& signal, 
               Util::Image& d0, Util::Image& d1, Util::Image& d2)
  {
    const float sigma = (sqrt(2)*r) / 3.0;

    Util::Image d0Filt,d1Filt,d2Filt;
    FilterBank::createGaussianKernel(sigma, sigma, 3.0, theta, 0, false, d0Filt);
    FilterBank::createGaussianKernel(sigma, sigma, 3.0, theta, 1, false, d1Filt);
    FilterBank::createGaussianKernel(sigma, sigma, 3.0, theta, 2, false, d2Filt);

    getFiltered(signal,d0Filt,d0);
    getFiltered(signal,d1Filt,d1);
    getFiltered(signal,d2Filt,d2);
  }

} //namespace Group


