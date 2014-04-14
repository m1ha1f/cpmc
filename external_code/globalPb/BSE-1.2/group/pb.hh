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


#ifndef PB_HH
#define PB_HH

#include "nonmax.hh"
#include "image.hh"
#include "string.hh"
#include "texture.hh"

namespace Group
{
  class Pb 
  {
    public:

      //
      // configuration structure and interface 
      //
      struct Config
      {
          // filterbank for tg
          int tgFiltNumOrient;
          int tgFiltNumScales;
          float tgFiltStartScale;
          float tgFiltScaleFactor;
          bool gaussDeriv;
          float tgEpsilon;

          // textons
          Util::String textonFile;
          int textonStartChannels;
          float textonChannelFactor;

          // texture gradient
          bool tgScaleLinked;
          int tgNumOrient;
          int tgNumScales;
          float tgStartScale;
          float tgScaleFactor;

          // color gradient
          int cgNumOrient;
          int cgNumScales;
          float cgStartScale;
          float cgScaleFactor;
          int    cgBinsA, cgBinsB;
          float cgSigma, cgSupport, cgZoom;

          // brightness gradient
          int bgNumOrient;
          int bgNumScales;
          float bgStartScale;
          float bgScaleFactor;
          int bgBinsL;
          float bgSigma, bgSupport, bgZoom;
      };

      static void registerConfig ();
      static void getDefaultConfig (Config & config);
      const Config& getConfig () const;

      Pb ();
      Pb (const Config & config);
      ~Pb ();

      //
      // initialize does all the computationally expensive
      // image  processing in order to compute the gradients
      //
      void initialize(const Util::ImageStack& im, const bool useColor = true);

      //
      // compute probability of a boundary at a given location and
      // orientation.  this assumes the default configuration parameters.
      // in particular, it only uses the zeroth scale.
      //
      void computeOrientedPb(const int numOrient, Util::ImageStack& pbs);
 
      //
      // compute probability of a boundary at a given location and
      // orientation.  this assumes the default configuration parameters.
      // in particular, it only uses the zeroth scale.
      //
      void computeOrientedNMSPb(const int numOrient, Util::ImageStack& pbs);

      //
      // compute the pb image from the cached gradients
      //
      void computePb(const int numOrient, DualLattice& boundaries);

      //
      // compute the pb image from the cached gradients
      //
      void computePb(const int numOrient, Util::Image& boundaries);

      //
      // access texton map
      //
      const TextonMap& getTextons(const int scale) const;

      //
      // access to raw gradient info
      // 
      const Util::Image& getTG(const int scale, const int orient) const;
      const Util::Image& getCG(const int scale, const int orient) const;
      const Util::Image& getBG(const int scale, const int orient) const;

      //
      // hints for acessing TG, CG, and BG fields
      //
      int getTGNumScales() const;
      int getTGNumOrient() const;
      int getCGNumScales() const;
      int getCGNumOrient() const;
      int getBGNumScales() const;
      int getBGNumOrient() const;

    protected:
      // Configuration
      Config config;
      void checkConfig ();
      static bool registerConfigCalled;

      //
      // core computations
      //
      void computeTG();
      void computeCG();
      void computeBG();

      //
      // helper functions for computing features and pb
      //
      float tgFitRadius (const int scale) const;
      float cgFitRadius (const int scale) const;
      float bgFitRadius (const int scale) const;

      static void computeNeedleness (const float epsilon, const Util::Image& a, const Util::Image& b,
                        const Util::Image& c, Util::Image& needles);

      //
      // member vars
      //
      Util::ImageStack m_lab;
      int m_width;
      int m_height;
      float m_idiag;
      bool m_useColor;
      bool m_initialized;

      Util::Array1D<TextonMap> m_textons;         // Texton Maps 

      Util::ImageStack m_tg;
      Util::ImageStack m_cg;
      Util::ImageStack m_bg;
  };

} //namespace Group

#endif 


