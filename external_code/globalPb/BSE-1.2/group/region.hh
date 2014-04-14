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


#ifndef REGION_H 
#define REGION_H 

#include "string.hh"
#include "image.hh"
#include "texture.hh"

namespace Group
{
  class Region
  {
    public:
      // Configuration
      struct Config
      {
          // filterbank for tp
          int tpFiltNumOrient;
          int tpFiltNumScales;
          float tpFiltStartScale;
          float tpFiltScaleFactor;

          // textons
          Util::String textonFile;
          int textonStartChannels;
          float textonChannelFactor;

          // texture patch 
          bool tpScaleLinked;
          int tpNumScales;
          float tpStartScale;
          float tpScaleFactor;

          // color patch
          int cpNumScales;
          float cpStartScale;
          float cpScaleFactor;
          int cpBinsA, cpBinsB;
          float cpSigma, cpSupport, cpZoom;

          // brightness patch 
          int bpNumScales;
          float bpStartScale;
          float bpScaleFactor;
          int bpBinsL;
          float bpSigma, bpSupport, bpZoom;
      };

      static void registerConfig ();
      static void getDefaultConfig (Config & config);
      const Config & getConfig () const;

      Region ();
      Region (const Config & config);
      ~Region ();

      //
      // initialize does all the computationally expensive
      // image processing in order to compute the region
      // similarity cues
      //
      void initialize(const Util::ImageStack& im);
      void initialize(const Util::ImageStack& im, const SupportMap& support);

      //
      // a hint about how big the support map should be
      // to match the starting patch scale given in the 
      // config parameters
      //
      int getSupportMapRadius(const float imagediag); 

      //
      // get the probability that two patches came from the
      // same segment.  fit to default parameter settings.
      //
      float getPss(const int scale, const int x1, const int y1, const int x2, const int y2) const; 

      //
      // access the histograms
      //
      const Histogram& getTextonHist(const int scale) const;
      const Histogram& getLHist(const int scale) const;
      const Histogram& getAHist(const int scale) const;
      const Histogram& getBHist(const int scale) const;

      //
      // get the similarity between two points
      //
      float getTSim(const int scale, const int x1, const int y1, const int x2, const int y2) const;
      float getBSim(const int scale, const int x1, const int y1, const int x2, const int y2) const;
      float getCSim(const int scale, const int x1, const int y1, const int x2, const int y2) const;

      int getTPNumScales() const;
      int getBPNumScales() const;
      int getCPNumScales() const;

      bool useColor() const;

    protected:
      // Configuration
      Config config;
      void checkConfig ();
      static bool registerConfigCalled;
     
      bool m_initialized;

      //
      // core computations
      //
      void computeTextonHistograms ();
      void computeColorHistograms ();
      void computeBrightnessHistograms ();
      void computeTextonHistograms(const SupportMap& support);
      void computeColorHistograms(const SupportMap& support);
      void computeBrightnessHistograms(const SupportMap& support);

      //
      // image properties
      //
      Util::ImageStack m_lab;
      int m_height;
      int m_width;
      float m_idiag;
      bool m_useColor;

      //
      // cached histograms
      //
      Util::Array1D<Histogram> m_textonHist;
      Util::Array1D<Histogram> m_lHist;
      Util::Array1D<Histogram> m_aHist;
      Util::Array1D<Histogram> m_bHist;
  };

} //namespace Group

#endif 

