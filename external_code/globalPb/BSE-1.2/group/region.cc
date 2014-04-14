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

#include "configure.hh"
#include "exception.hh"
#include "message.hh"
#include "string.hh"
#include "array.hh"
#include "util.hh"

#include "image.hh"
#include "filterbank.hh"
#include "texture.hh"
#include "color.hh"
#include "region.hh"

namespace Group
{

  //////////////////////////////////////////////////////////////////////////////
  //
  // configuration structure and interface
  //

  bool Region::registerConfigCalled = false;

  void
  Region::registerConfig ()
  {
    if (registerConfigCalled)
    {
        return;
    }

    // filterbank for tp
    Configure::registerInt ("Region::tp::filt::norient", 6,
                            "Number of orientations for filterbank.");
    Configure::registerInt ("Region::tp::filt::nscales", 1,
                            "Number of scales for filterbank.");
    Configure::registerFloat ("Region::tp::filt::scalefactor", sqrt (2),
                              "Scale factor for filterbank.");
    Configure::registerFloat ("Region::tp::filt::startscale", 0.007,
                              "Beginning scale for filterbank, as fraction of image diagonal.");

    // textons
    Configure::registerString ("Region::textons::file", NULL,
                               "If provided, textons are read from the file rather than computed.");
    Configure::registerInt ("Region::textons::startchannels", 12,
                            "Beginning number of texton channels.");
    Configure::registerFloat ("Region::textons::channelfactor", 2,
                              "Increase number of textons by this factor for each scale.");

    // texture patch 
    Configure::registerBool ("Region::tp::scalelinked", false,
                             "Use textons computed from different scale filters for different "
                             "scale texture patches?");
    Configure::registerInt ("Region::tp::nscales", 1,
                            "Number of scales for texture patch.");
    Configure::registerFloat ("Region::tp::startscale", 0.056,
                              "Starting scale for texture patch.");
    Configure::registerFloat ("Region::tp::scalefactor", sqrt (2),
                              "Ratio between successive texture scales.");

    // color patch 
    Configure::registerInt ("Region::cp::nscales", 1,
                            "Number of scales for CP.");
    Configure::registerFloat ("Region::cp::startscale", 0.02,
                              "Starting scale for CP.");
    Configure::registerFloat ("Region::cp::scalefactor", sqrt (2),
                              "Ratio between successive CP scales.");
    Configure::registerInt ("Region::cp::bins::a", 25,
                            "Number of bins for L channel.");
    Configure::registerInt ("Region::cp::bins::b", 25,
                            "Number of bins for L channel.");
    Configure::registerFloat ("Region::cp::sigma", 0.1,
                              "Sigma for CP histogram kernel.");
    Configure::registerFloat ("Region::cp::support", 2,
                              "Support of CP histogram kernel, in units of sigma.");
    Configure::registerFloat ("Region::cp::zoom", 5,
                              "Resolution of CP histogram kernel, in units of points/sigma.");

    // brightness patch 
    Configure::registerInt ("Region::bp::nscales", 1,
                            "Number of scales for BP.");
    Configure::registerFloat ("Region::bp::startscale", 0.02,
                              "Starting scale for BP.");
    Configure::registerFloat ("Region::bp::scalefactor", sqrt (2),
                              "Ratio between successive BP scales.");
    Configure::registerInt ("Region::bp::bins::L", 6,
                            "Number of bins for L channel.");
    Configure::registerFloat ("Region::bp::sigma", 0.40,
                              "Sigma for BP histogram kernel.");
    Configure::registerFloat ("Region::bp::support", 2,
                              "Support of BP histogram kernel, in units of sigma.");
    Configure::registerFloat ("Region::bp::zoom", 5,
                              "Resolution of BP histogram kernel, in units of points/sigma.");

    registerConfigCalled = true;
  }

  void
  Region::getDefaultConfig (Config & config)
  {

      // filterbank for tp
      config.tpFiltNumOrient = Configure::getInt ("Region::tp::filt::norient");
      config.tpFiltNumScales = Configure::getInt ("Region::tp::filt::nscales");
      config.tpFiltScaleFactor =
          Configure::getFloat ("Region::tp::filt::scalefactor");
      config.tpFiltStartScale =
          Configure::getFloat ("Region::tp::filt::startscale");

      // textons
      config.textonFile = Configure::getString ("Region::textons::file");
      config.textonStartChannels =
          Configure::getInt ("Region::textons::startchannels");
      config.textonChannelFactor =
          Configure::getFloat ("Region::textons::channelfactor");

      // texture patch 
      config.tpScaleLinked = Configure::getBool ("Region::tp::scalelinked");
      config.tpNumScales = Configure::getInt ("Region::tp::nscales");
      config.tpStartScale = Configure::getFloat ("Region::tp::startscale");
      config.tpScaleFactor = Configure::getFloat ("Region::tp::scalefactor");

      // color patch 
      config.cpNumScales = Configure::getInt ("Region::cp::nscales");
      config.cpStartScale = Configure::getFloat ("Region::cp::startscale");
      config.cpScaleFactor = Configure::getFloat ("Region::cp::scalefactor");
      config.cpBinsA = Configure::getInt ("Region::cp::bins::a");
      config.cpBinsB = Configure::getInt ("Region::cp::bins::b");
      config.cpSigma = Configure::getFloat ("Region::cp::sigma");
      config.cpSupport = Configure::getFloat ("Region::cp::support");
      config.cpZoom = Configure::getFloat ("Region::cp::zoom");

      // brightness patch 
      config.bpNumScales = Configure::getInt ("Region::bp::nscales");
      config.bpStartScale = Configure::getFloat ("Region::bp::startscale");
      config.bpScaleFactor = Configure::getFloat ("Region::bp::scalefactor");
      config.bpBinsL = Configure::getInt ("Region::bp::bins::L");
      config.bpSigma = Configure::getFloat ("Region::bp::sigma");
      config.bpSupport = Configure::getFloat ("Region::bp::support");
      config.bpZoom = Configure::getFloat ("Region::bp::zoom");

  }


  void
  Region::checkConfig ()
  {
    // TODO: if any of the configuration parameters are invlaid, then 
    // throw an exception describing the problem.
  }

  const Region::Config& Region::getConfig() const
  {
    return config;    
  }

  Region::Region ()
  {
    getDefaultConfig(config);
  }

  Region::Region (const Config & config)
  {
    this->config = config;
  }

  Region::~Region () 
  {
  }
  
  /////////////////////////////////////////////////////////////////////////////////////
  // initialize precomputes the histograms.  it can be called multiple
  // times with different images if desired.

  //
  // a hint about how big the support map should be
  // to match the patch scale given in the config. parameters
  //
  int Region::getSupportMapRadius(const float imagediag) 
  {
      const float tpScaleVal = 
        imagediag * config.tpStartScale * pow (config.tpScaleFactor, config.tpNumScales-1);
      const float cpScaleVal = 
        imagediag * config.cpStartScale * pow (config.cpScaleFactor, config.cpNumScales-1);
      const float bpScaleVal = 
        imagediag * config.bpStartScale * pow (config.bpScaleFactor, config.bpNumScales-1);

      const float maxscale = Util::max(Util::max(tpScaleVal,cpScaleVal),bpScaleVal); 
      int radius = ((int)ceil(maxscale)) + 1;
      return radius;
  }

  //
  // initialize with a support map
  //
  void
  Region::initialize (const Util::ImageStack& im, const SupportMap& support)
  {
      checkConfig();

      m_width = im.size(1);
      m_height = im.size(2);
      m_idiag = sqrt(m_width * m_width + m_height * m_height);

      if (im.size(0) == 3)
      {
        m_useColor = true;
        rgb2lab(im,m_lab);
        labNormalize(m_lab);
      }
      else
      {
        m_useColor = false;
        m_lab = im;
      }

      //precompute histograms
      Util::Message::startBlock("Region initialization");
      
      Util::Message::startBlock("texture histogram computation"); 
      computeTextonHistograms(support);
      Util::Message::endBlock();
     
      if (m_useColor)
      {
        Util::Message::startBlock("color histogram computation"); 
        computeColorHistograms(support);
        Util::Message::endBlock();
      }

      Util::Message::startBlock("brightness histogram computation"); 
      computeBrightnessHistograms(support);
      Util::Message::endBlock();

      Util::Message::endBlock();

      m_initialized = true;
  }

  //
  // initialize without a support map
  //
  void
  Region::initialize (const Util::ImageStack& im)
  {
      checkConfig();

      m_width = im.size(1);
      m_height = im.size(2);
      m_idiag = sqrt(m_width * m_width + m_height * m_height);

      if (im.size(0) == 3)
      {
        m_useColor = true;
        rgb2lab(im,m_lab);
        labNormalize(m_lab);
      }
      else
      {
        m_useColor = false;
        m_lab = im;
      }

      //precompute histograms
      Util::Message::startBlock("Region initialization");
      
      Util::Message::startBlock("texture histogram computation"); 
      computeTextonHistograms();
      Util::Message::endBlock();
     
      if (m_useColor)
      {
        Util::Message::startBlock("color histogram computation"); 
        computeColorHistograms();
        Util::Message::endBlock();
      }

      Util::Message::startBlock("brightness histogram computation"); 
      computeBrightnessHistograms();
      Util::Message::endBlock();

      Util::Message::endBlock();

      m_initialized = true;
  }


  //////////////////////////////////////////////////////////////////////////////
  //
  // patch classifier.  parameters are fit with default 
  //

  //color images
  static const float bp_stdev = 0.151709;
  static const float tp_stdev = 0.0875824;
  static const float cp_stdev = 0.225726;
  float cpss_model(const Util::Array1D<float> features)
  {
    float val = 0.85422;
    val += features(0) / bp_stdev * -0.105818; /* bp0 */
    val += features(1) / tp_stdev * -1.02351;  /* tp0 */
    val += features(2) / cp_stdev * -0.268088; /* cp0 */
    float post = 1.0 / (1.0 + exp(-val));
    if (!finite(post)) { return 0; }
    post = Util::minmax(0.0f,post,1.0f);
    return post;
  }

  //grayscale
  float gpss_model(const Util::Array1D<float> features)
  {
    float val = 2.14297;
    val += features(0) / bp_stdev * -1.0145;   /* bp0 */
    val += features(1) / tp_stdev * -1.34457;  /* tp0 */
    float post = 1.0 / (1.0 + exp(-val));
    if (!finite(post)) { return 0; }
    post = Util::minmax(0.0f,post,1.0f);
    return post;
  }

  float Region::getPss(const int scale, const int x1, const int y1, const int x2, const int y2) const
  {
    assert(m_initialized);
    float pval = 0;
    if (m_useColor)
    {
      Util::Array1D<float> features(3);
      features(0) = getBSim(scale,x1,y1,x2,y2);
      features(1) = getTSim(scale,x1,y1,x2,y2);
      features(2) = getCSim(scale,x1,y1,x2,y2);
      pval = cpss_model(features);
    }
    else
    {
      Util::Array1D<float> features(2);
      features(0) = getBSim(scale,x1,y1,x2,y2);
      features(1) = getTSim(scale,x1,y1,x2,y2);
      pval = gpss_model(features);
    }
    return pval;
  }


  //////////////////////////////////////////////////////////////////////////
  //
  // access histograms computed for the image.
  // 

  const Histogram& Region::getTextonHist(const int scale) const
  {
    assert(m_initialized);
    assert(scale >= 0 && scale < config.tpNumScales);
    return m_textonHist(scale);
  }
  const Histogram& Region::getLHist(const int scale) const
  {
    assert(m_initialized);
    assert(scale >= 0 && scale < config.bpNumScales);
    return m_lHist(scale);
  }
  const Histogram& Region::getAHist(const int scale) const
  {
    assert(m_initialized);
    assert(scale >= 0 && scale < config.cpNumScales);
    return m_aHist(scale);
  }
  const Histogram& Region::getBHist(const int scale) const
  {
    assert(m_initialized);
    assert(scale >= 0 && scale < config.cpNumScales);
    return m_bHist(scale);
  }


  //////////////////////////////////////////////////////////////////////////
  //
  // access similarities computed for the image.
  // 

  float Region::getTSim(const int scale, const int x1, const int y1, const int x2, const int y2) const
  {
    assert(m_initialized);
    float similarity = 0;
    textureSimilarity(m_textonHist(scale),x1,y1,x2,y2,similarity);
    return similarity;
  }

  float Region::getBSim(const int scale, const int x1, const int y1, const int x2, const int y2) const
  {
    assert(m_initialized);
    float similarity = 0;
    colorSimilarity(m_lHist(scale),x1,y1,x2,y2,similarity);    
    return similarity;
  }

  float Region::getCSim(const int scale, const int x1, const int y1, const int x2, const int y2) const
  {
    assert(m_initialized);
    float sa = 0;
    float sb = 0;    
    colorSimilarity(m_aHist(scale),x1,y1,x2,y2,sa);    
    colorSimilarity(m_bHist(scale),x1,y1,x2,y2,sb);    
    return 0.5 * (sa + sb);
  }


  //////////////////////////////////////////////////////////////////////////
  //
  // hints for accessing TSim,BSim,CSim
  //
  int Region::getTPNumScales() const
  {
    return config.tpNumScales;
  }

  int Region::getBPNumScales() const
  {
    return config.bpNumScales;
  }

  int Region::getCPNumScales() const
  {
    return config.cpNumScales;
  }

  
  //////////////////////////////////////////////////////////////////////////
 
  bool Region::useColor() const
  {
    assert(m_initialized);
    return m_useColor;    
  }
  
  //////////////////////////////////////////////////////////////////////////

  void
  Region::computeTextonHistograms()
  {
    m_textonHist.resize(config.tpNumScales);
    for (int scale = 0; scale < config.tpNumScales; scale++)
    {
        float filtStartScale = m_idiag * config.tpFiltStartScale;
        if (config.tpScaleLinked)
        {
          filtStartScale *= pow (config.tpScaleFactor, scale);
        }
     
        // create the filterbank
        FilterBank fbank(config.tpFiltNumOrient, config.tpFiltNumScales,
                         filtStartScale, config.tpFiltScaleFactor);

        Util::Message::debug(Util::String("Filterbank %d %d %f %f",
           config.tpFiltNumOrient,config.tpFiltNumScales, filtStartScale,config.tpFiltScaleFactor),2);

        //extract textons
        TextonMap textons;
        textons.map.resize(m_lab.size(1),m_lab.size(2));
        textons.numChannels = 0;

        if(config.textonFile == "")
        {
          const int nchannels = (int)rint(config.textonStartChannels*pow(config.textonChannelFactor,scale));
          computeTextons(*m_lab.slice(Util::LAB_L),fbank,nchannels,textons);
        }
        else
        {
          computeTextons(*m_lab.slice(Util::LAB_L),fbank,config.textonFile,textons);
        }
        Util::Message::debug(Util::String("Used %d textons",textons.numChannels),2);

        //now compute histograms
        const float tpScaleVal = m_idiag * config.tpStartScale * pow (config.tpScaleFactor, scale);
        m_textonHist(scale).resize(m_width,m_height,textons.numChannels);
        Group::computeTextonHistograms(textons,tpScaleVal,m_textonHist(scale)); 
    } 
  }  

  void
  Region::computeColorHistograms ()
  {
    m_aHist.resize(config.cpNumScales);
    m_bHist.resize(config.cpNumScales);

    // compute TP at each scale
    for (int scale = 0; scale < config.cpNumScales; scale++)
    {
      const float cpScaleVal = m_idiag * config.cpStartScale * pow (config.cpScaleFactor, scale);
      m_aHist(scale).resize(m_width,m_height,config.cpBinsA);
      m_bHist(scale).resize(m_width,m_height,config.cpBinsB);
      Group::computeColorHistograms(*m_lab.slice(Util::LAB_A), config.cpBinsA, cpScaleVal, config.cpSigma,
                                                config.cpSupport, config.cpZoom, m_aHist(scale));
      Group::computeColorHistograms(*m_lab.slice(Util::LAB_B), config.cpBinsB, cpScaleVal, config.cpSigma,
                                                config.cpSupport, config.cpZoom, m_bHist(scale));
    }
  }



  void
  Region::computeBrightnessHistograms ()
  {
      m_lHist.resize(config.bpNumScales);
      for (int scale = 0; scale < config.bpNumScales; scale++)
      {
        const float bpScaleVal = m_idiag * config.bpStartScale * pow (config.bpScaleFactor, scale);
        m_lHist(scale).resize(m_width,m_height,config.bpBinsL);
        Group::computeColorHistograms(*m_lab.slice(Util::LAB_L), config.bpBinsL, bpScaleVal, config.bpSigma,
                                        config.bpSupport, config.bpZoom, m_lHist(scale));
      }
  }


  ///////////////////////////////////////////////////////////////////////////
  
  void
  Region::computeTextonHistograms(const SupportMap& support)
  {
    m_textonHist.resize(1);
    float filtStartScale = m_idiag * config.tpFiltStartScale;
     
    // create the filterbank
    FilterBank fbank(config.tpFiltNumOrient, config.tpFiltNumScales,
                     filtStartScale, config.tpFiltScaleFactor);

    Util::Message::debug(Util::String("Filterbank %d %d %f %f",
       config.tpFiltNumOrient,config.tpFiltNumScales, filtStartScale,config.tpFiltScaleFactor),2);

    //extract textons
    TextonMap textons;
    textons.map.resize(m_lab.size(1),m_lab.size(2));
    textons.numChannels = 0;

    if(config.textonFile == "")
    {
      const int nchannels = config.textonStartChannels;
      computeTextons(*m_lab.slice(Util::LAB_L),fbank,nchannels,textons);
    }
    else
    {
      computeTextons(*m_lab.slice(Util::LAB_L),fbank,config.textonFile,textons);
    }
    Util::Message::debug(Util::String("Used %d textons",textons.numChannels),2);

    m_textonHist(0).resize(m_width,m_height,textons.numChannels);

    Group::computeTextonHistograms(textons,support,m_textonHist(0)); 
  }  

  void
  Region::computeColorHistograms (const SupportMap& support)
  {
    m_aHist.resize(1);
    m_bHist.resize(1);
    m_aHist(0).resize(m_width,m_height,config.cpBinsA);
    m_bHist(0).resize(m_width,m_height,config.cpBinsB);
    Group::computeColorHistograms(*m_lab.slice(Util::LAB_A), config.cpBinsA, config.cpSigma,
                                    config.cpSupport, config.cpZoom, support, m_aHist(0));
    Group::computeColorHistograms(*m_lab.slice(Util::LAB_B), config.cpBinsB, config.cpSigma,
                                    config.cpSupport, config.cpZoom, support, m_bHist(0));
  }

  void
  Region::computeBrightnessHistograms(const SupportMap& support)
  {
      m_lHist.resize(1);
      m_lHist(0).resize(m_width,m_height,config.bpBinsL);
      Group::computeColorHistograms(*m_lab.slice(Util::LAB_L), config.bpBinsL, config.bpSigma,
                                      config.bpSupport, config.bpZoom, support, m_lHist(0));
  }

} //namespace Group
