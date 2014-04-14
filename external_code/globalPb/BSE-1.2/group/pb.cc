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


#include <float.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "configure.hh"
#include "exception.hh"
#include "message.hh"
#include "string.hh"
#include "array.hh"
#include "util.hh"

#include "image.hh"
#include "filterbank.hh"
#include "nonmax.hh"
#include "texture.hh"
#include "color.hh"
#include "pb.hh"

namespace Group
{

  //////////////////////////////////////////////////////////////////////////////
  // 
  // configuration structure and interface
  //

  bool Pb::registerConfigCalled = false;

  void Pb::registerConfig()
  {
    if(registerConfigCalled)
    {
      return;
    }

    
    // filterbank for tg
    Configure::registerInt("Pb::tg::filt::norient",
                           6, "Number of orientations for filterbank.");
    Configure::registerInt("Pb::tg::filt::nscales", 1,
                           "Number of scales for filterbank.");
    Configure::registerFloat("Pb::tg::filt::scalefactor", sqrt(2),
                             "Scale factor for filterbank.");
    Configure::registerFloat("Pb::tg::filt::startscale", 0.007,
                             "Beginning scale for filterbank, as fraction of image diagonal.");

    // textons
    Configure::registerString("Pb::textons::file", NULL,
                              "If provided, textons are read from the file rather than computed.");
    Configure::registerInt("Pb::textons::startchannels", 12,
                           "Beginning number of texton channels.");
    Configure::registerFloat("Pb::textons::channelfactor", 2,
                             "Increase number of textons by this factor for each scale.");

    // texture gradient
    Configure::registerBool("Pb::tg::scalelinked", false,
                            "Use textons computed from different scale filters for different "
                            "scale texture gradients?");
    Configure::registerInt("Pb::tg::norient", 8,
                           "Number of orientations for texture gradient.");
    Configure::registerInt("Pb::tg::nscales", 1,
                           "Number of scales for texture gradient.");
    Configure::registerFloat("Pb::tg::startscale", 0.0198,
                             "Starting scale for texture gradient.");
    Configure::registerFloat("Pb::tg::scalefactor", sqrt(2),
                             "Ratio between successive texture scales.");
    Configure::registerBool("Pb::tg::gaussDeriv", false,
                            "Use gaussian derivative estimates (instead of Savitsky-Golay filtering)?");
    Configure::registerFloat("Pb::tg::epsilon", 0.01,
                             "Epsilon for needleness function.");

    // color gradient
    Configure::registerInt("Pb::cg::norient", 8,
                           "Number of orientations for CG.");
    Configure::registerInt("Pb::cg::nscales", 1,
                           "Number of scales for CG.");
    Configure::registerFloat("Pb::cg::startscale", 0.0198,
                             "Starting scale for CG.");
    Configure::registerFloat("Pb::cg::scalefactor", sqrt(2),
                             "Ratio between successive CG scales.");
    Configure::registerInt("Pb::cg::bins::a", 25,
                           "Number of bins for a channel.");
    Configure::registerInt("Pb::cg::bins::b", 25,
                           "Number of bins for b channel.");
    Configure::registerFloat("Pb::cg::sigma", 0.10,
                             "Sigma for CG histogram kernel.");
    Configure::registerFloat("Pb::cg::support", 2,
                             "Support of CG histogram kernel, in units of sigma.");
    Configure::registerFloat("Pb::cg::zoom", 5,
                             "Resolution of CG histogram kernel, in units of points/sigma.");

    // brightness gradient
    Configure::registerInt("Pb::bg::norient", 8,
                           "Number of orientations for BG.");
    Configure::registerInt("Pb::bg::nscales", 1,
                           "Number of scales for BG.");
    Configure::registerFloat("Pb::bg::startscale", 0.0106,
                             "Starting scale for BG.");
    Configure::registerFloat("Pb::bg::scalefactor", sqrt(2),
                             "Ratio between successive BG scales.");
    Configure::registerInt("Pb::bg::bins::L", 12,
                           "Number of bins for L channel.");
    Configure::registerFloat("Pb::bg::sigma", 0.20,
                             "Sigma for BG histogram kernel.");
    Configure::registerFloat("Pb::bg::support", 2,
                             "Support of BG histogram kernel, in units of sigma.");
    Configure::registerFloat("Pb::bg::zoom", 5,
                             "Resolution of BG histogram kernel, in units of points/sigma.");


    registerConfigCalled = true;
  }


  void Pb::getDefaultConfig(Config & config)
  {

    // filterbank for tg
    config.tgFiltNumOrient = Configure::getInt("Pb::tg::filt::norient");
    config.tgFiltNumScales = Configure::getInt("Pb::tg::filt::nscales");
    config.tgFiltScaleFactor =
      Configure::getFloat("Pb::tg::filt::scalefactor");
    config.tgFiltStartScale =
      Configure::getFloat("Pb::tg::filt::startscale");
    config.gaussDeriv = Configure::getBool("Pb::tg::gaussDeriv");
    config.tgEpsilon = Configure::getFloat("Pb::tg::epsilon");

    // textons
    config.textonFile = Configure::getString("Pb::textons::file");
    config.textonStartChannels =
      Configure::getInt("Pb::textons::startchannels");
    config.textonChannelFactor =
      Configure::getFloat("Pb::textons::channelfactor");

    // texture gradient
    config.tgScaleLinked = Configure::getBool("Pb::tg::scalelinked");
    config.tgNumOrient = Configure::getInt("Pb::tg::norient");
    config.tgNumScales = Configure::getInt("Pb::tg::nscales");
    config.tgStartScale = Configure::getFloat("Pb::tg::startscale");
    config.tgScaleFactor = Configure::getFloat("Pb::tg::scalefactor");

    // color gradient
    config.cgNumOrient = Configure::getInt("Pb::cg::norient");
    config.cgNumScales = Configure::getInt("Pb::cg::nscales");
    config.cgStartScale = Configure::getFloat("Pb::cg::startscale");
    config.cgScaleFactor = Configure::getFloat("Pb::cg::scalefactor");
    config.cgBinsA = Configure::getInt("Pb::cg::bins::a");
    config.cgBinsB = Configure::getInt("Pb::cg::bins::b");
    config.cgSigma = Configure::getFloat("Pb::cg::sigma");
    config.cgSupport = Configure::getFloat("Pb::cg::support");
    config.cgZoom = Configure::getFloat("Pb::cg::zoom");

    // brightness gradient
    config.bgNumOrient = Configure::getInt("Pb::bg::norient");
    config.bgNumScales = Configure::getInt("Pb::bg::nscales");
    config.bgStartScale = Configure::getFloat("Pb::bg::startscale");
    config.bgScaleFactor = Configure::getFloat("Pb::bg::scalefactor");
    config.bgBinsL = Configure::getInt("Pb::bg::bins::L");
    config.bgSigma = Configure::getFloat("Pb::bg::sigma");
    config.bgSupport = Configure::getFloat("Pb::bg::support");
    config.bgZoom = Configure::getFloat("Pb::bg::zoom");

  }

  void Pb::checkConfig()
  {
    // TODO: if any of the configuration parameters are invalid, then 
    // throw an exception describing the problem.
  }

  const Pb::Config& Pb::getConfig () const
  {
    return config;
  }

  Pb::Pb()
  {
    getDefaultConfig(config);
  }

  Pb::Pb(const Config & config)
  {
    this->config = config;
  }

  Pb::~Pb()
  {
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // initialize does all the computationally expensive
  // image processing in order to compute the various
  // gradient features.  it can be called multiple times
  // with different images if desired....
  //
  
  void Pb::initialize(const Util::ImageStack& im, const bool useColor)
  {
    checkConfig();

    m_useColor = useColor;
    m_width = im.size(1);
    m_height = im.size(2);
    m_idiag = sqrt(m_width * m_width + m_height * m_height);
    
    if (im.size(0) == 3)
    {
      rgb2lab(im,m_lab);
      labNormalize(m_lab);
    } 
    else
    {
      if (m_useColor == true)
      {
        Util::Message::debug("WARNING: cowardly refusing to compute color pb for grayscale image");  
      }
      m_lab = im;
      m_useColor = false;
    }

    Util::Message::startBlock("Pb initialization");

    Util::Message::startBlock("brightness gradient computation");
    computeBG();
    Util::Message::endBlock();
    if (m_useColor)
    {
      Util::Message::startBlock("color gradient computation");
      computeCG();
      Util::Message::endBlock();
    }

    Util::Message::startBlock("texture gradient computation");
    computeTG();
    Util::Message::endBlock();

    Util::Message::endBlock();

    m_initialized = true;
  }

  //////////////////////////////////////////////////////////////////////////////
  //
  // compute the pb image from the cached gradients
  //

  //
  // logistic fits.  this is of course only valid for the default 
  // parameter settings!!  changing relative scales etc. can make
  // these fits "non-optimal"
  //
  static const float bg_stdev = 0.13;
  static const float cg_stdev = 0.077;
  static const float tg_stdev = 0.063;

  static float cpb_model(const Util::Array1D<float> features)
  {
    float val = -3.08;
    val += features(0) / bg_stdev * 0.31;  //bg
    val += features(1) / tg_stdev * 0.53;  //cg
    val += features(2) / cg_stdev * 0.44;  //tg
    float post = 1.0 / (1.0 + exp(-val));
    if(!finite(post)) { return 0; }
    post = Util::minmax(0.0f,post,1.0f);
    return post;
  }

  static float gpb_model(const Util::Array1D<float> features)
  {
    float val = -2.81;
    val += features(0) / bg_stdev * 0.50;  // bg 
    val += features(1) / tg_stdev * 0.52;  // tg 
    float post = 1.0 / (1.0 + exp(-val));
    if(!finite(post)) { return 0; }
    post = Util::minmax(0.0f,post,1.0f);
    return post;
  }
  
  static int closestOrient(const float theta, const int norient)
  {
    assert(theta >= 0);
    int orient = (int)rint((theta / M_PI) * norient);
    assert(orient >= 0);
    return (orient % norient);
  }

  //
  // compute probability of a boundary at a given location and
  // orientation.  this assumes the default configuration parameters.  
  // in particular, it only uses the zeroth scale.
  //
  void Pb::computeOrientedPb(const int numOrient, Util::ImageStack& pbs)
  {
    assert(m_initialized);
    pbs.resize(numOrient,m_width,m_height);
    Util::Array1D<float> features;
    if(m_useColor)
    {
      features.resize(3);
    }
    else
    {
      features.resize(2);
    }

    features.init(0);
    for(int orient = 0; orient < numOrient; orient++)
    {
      const float theta = (M_PI*orient / numOrient);
      const int bgorient = closestOrient(theta, config.bgNumOrient);
      const int cgorient = closestOrient(theta, config.cgNumOrient);
      const int tgorient = closestOrient(theta, config.tgNumOrient);
      for(int x = 0; x < m_width; x++)
      {
        for(int y = 0; y < m_height; y++)
        {
          float p = 0;
          if (m_useColor)
          {
            features(0) = m_bg(bgorient,x,y);
            features(1) = m_tg(tgorient,x,y);
            features(2) = m_cg(cgorient,x,y);
            p = cpb_model(features);
          }
          else
          {
            features(0) = m_bg(bgorient,x,y);
            features(1) = m_tg(tgorient,x,y);
            p = gpb_model(features);
          }
          assert(p >= 0 && p <= 1);
          pbs(orient,x,y) = p;
        }
      }
    }
  }

  //
  // compute probability of a boundary at a given location and
  // orientation.  this assumes the default configuration parameters.  
  // in particular, it only uses the zeroth scale.
  //
  void Pb::computeOrientedNMSPb(const int numOrient, Util::ImageStack& pbs)
  {
    assert(m_initialized);
    computeOrientedPb(numOrient,pbs);

    const float oStep = M_PI / numOrient;
    const float r = 2.5;
    for (int orient = 0; orient < numOrient; orient++)
    {
      float theta = orient*oStep;
      assert(theta >= 0 && theta < M_PI);
      ParabolaSet parabs;
      fitSubpixelParabolas(*pbs.slice(orient), theta, r, parabs);
      Util::Image mag = parabs.smooth;

      Util::Image nms(pbs.size(1),pbs.size(2));
      nonmaxSuppress(mag,theta,nms);
      for (int x = 0; x < pbs.size(1); x++)
      {
        for (int y = 0; y < pbs.size(2); y++)
        {
          pbs(orient,x,y) = nms(x,y); 
        }
      }
    }                 
  }

  //
  // compute pb maxed over orientations
  //
  void Pb::computePb(const int numOrient, DualLattice& boundaries)
  {
    assert(m_initialized);
    Util::ImageStack pbs(numOrient,m_width,m_height);
    computeOrientedPb(numOrient,pbs);

    // non-maximal suppression and sub-pixel localization
    Util::Message::debug("fitting sub-pixel parabolas");
    //const float fitRad = 0.0104*m_idiag;
    const float fitRad = 0.01*m_idiag;
    ParabolaSet parabs;
    fitSubpixelParabolas(pbs,fitRad,parabs);

    // intersect needles to get final pb
    Util::Message::debug("projecting parabolas down to pixel grid");
    intersectPixels(parabs, boundaries);
  }

  void Pb::computePb(const int numOrient, Util::Image& boundaries)
  {
    assert(m_initialized);
    Util::ImageStack pbs(numOrient,m_width,m_height);
    computeOrientedPb(numOrient,pbs);

    // non-maximal suppression and sub-pixel localization
    Util::Message::debug("fitting sub-pixel parabolas");
    //const float fitRad = 0.0104*m_idiag;
    const float fitRad = 0.01*m_idiag;
    ParabolaSet parabs;
    fitSubpixelParabolas(pbs,fitRad,parabs);

    // intersect needles to get final pb
    Util::Message::debug("projecting parabolas down to pixel grid");
    intersectPixels(parabs, boundaries);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // access to raw gradient info and texton maps
  //

  const TextonMap& Pb::getTextons(const int scale) const
  {
    assert(m_initialized);
    assert(scale >=0 && scale < config.tgNumScales);
    return m_textons(scale);
  }
  
  const Util::Image& Pb::getTG(const int scale, const int orient) const
  {
    assert(m_initialized);
    assert(scale >= 0 && scale < config.tgNumScales);
    assert(orient >= 0 && orient < config.tgNumOrient);
    return *m_tg.slice(scale*config.tgNumOrient + orient);
  }

  const Util::Image& Pb::getCG(const int scale, const int orient) const
  {
    assert(m_initialized);
    assert(scale >= 0 && scale < config.cgNumScales);
    assert(orient >= 0 && orient < config.cgNumOrient);
    if (m_useColor)
    {
      return *m_cg.slice(scale*config.cgNumOrient + orient);
    }
    else
    {
      throw Util::Exception("no color gradient available");
    }
  }

  const Util::Image& Pb::getBG(const int scale, const int orient) const
  {
    assert(m_initialized);
    assert(scale >= 0 && scale < config.bgNumScales);
    assert(orient >= 0 && orient < config.bgNumOrient);
    return *m_bg.slice(scale*config.bgNumOrient + orient);
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // hints for accessing TG,CG,BG
  //
  
  inline int Pb::getTGNumScales() const
  {
    return config.tgNumScales;
  }
  inline int Pb::getTGNumOrient() const
  {
    return config.tgNumOrient;
  }
  inline int Pb::getCGNumScales() const
  {
    return config.cgNumScales;
  }
  inline int Pb::getCGNumOrient() const
  {
    return config.cgNumOrient;
  }
  inline int Pb::getBGNumScales() const
  {
    return config.bgNumScales;
  }
  inline int Pb::getBGNumOrient() const
  {
    return config.bgNumOrient;
  }


  //////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // helper functions for computing features and pb
  //

  //
  // compute "needleness" which is a function of the derivatives at a point
  //
  void Pb::computeNeedleness(const float epsilon, const Util::Image& a, const Util::Image& b, 
                                 const Util::Image& c, Util::Image& needleness)
  {
    const int width = a.size(0);
    const int height = a.size(1);
    needleness.resize(width,height);
    for(int x = 0; x < width; x++)
    {
      for(int y = 0; y < height; y++)
      {
        needleness(x,y) = Util::max(0.0f,a(x,y)) * Util::max(0.0f,-c(x,y)) / (fabs(b(x,y))+epsilon);
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // core computation
  //

  //
  // compute texture gradient
  //
  void Pb::computeTG()
  {
    m_tg.resize(config.tgNumScales*config.tgNumOrient,m_width,m_height);
    m_textons.resize(config.tgNumScales);
    Util::Image tgSmooth(m_width,m_height);
    Util::Image tgSlope(m_width,m_height);
    Util::Image tgCurv(m_width,m_height);

    // compute TG at each scale
    for(int tgScale = 0; tgScale < config.tgNumScales; tgScale++)
    {
      const float tgScaleVal = m_idiag * config.tgStartScale * pow(config.tgScaleFactor, tgScale);
      float filtStartScale = m_idiag * config.tgFiltStartScale;
      if(config.tgScaleLinked)
      {
        filtStartScale *= pow(config.tgScaleFactor, tgScale);
      }

      // create the filterbank
      FilterBank fbank(config.tgFiltNumOrient, config.tgFiltNumScales,
                       filtStartScale, config.tgFiltScaleFactor);

      Util::Message::debug(Util::String("Filterbank %d %d %f %f",
         config.tgFiltNumOrient,config.tgFiltNumScales, filtStartScale,config.tgFiltScaleFactor),2); 

      //extract textons
      m_textons(tgScale).map.resize(m_lab.size(1),m_lab.size(2));
      m_textons(tgScale).numChannels = 0;

      if(config.textonFile == "")
      {
        const int nchannels = (int)rint(config.textonStartChannels * pow(config.textonChannelFactor,tgScale));
        computeTextons(*m_lab.slice(Util::LAB_L),fbank,nchannels,m_textons(tgScale));
      }
      else
      {
        computeTextons(*m_lab.slice(Util::LAB_L),fbank,config.textonFile,m_textons(tgScale));
      }
      Util::Message::debug(Util::String("Used %d textons",m_textons(tgScale).numChannels),2);

      // compute the TG
      Util::ImageStack tgraw;
      Group::computeTG(m_textons(tgScale),tgScaleVal,config.tgNumOrient,tgraw);

      // compute localized TG
      Util::Message::startBlock(config.tgNumOrient,Util::String("TG %d smoothing",tgScale).text());
      for(int orient = 0; orient < config.tgNumOrient; orient++)
      {
        Util::Message::stepBlock();
        int index = tgScale * config.tgNumOrient + orient;
        const float r = tgFitRadius(tgScale);
        const float theta = M_PI * orient / config.tgNumOrient;
        if (config.gaussDeriv)
        {
          derivs(theta, r, *tgraw.slice(index), tgSmooth, tgSlope, tgCurv);
        }
        else
        {
          lowess(theta, r, *tgraw.slice(index), tgSmooth, tgSlope, tgCurv);
        }

        Util::Image needle;
        computeNeedleness(config.tgEpsilon, tgSmooth, tgSlope, tgCurv, needle);
        for (int x = 0; x < m_width; x++)
        {
          for (int y = 0; y < m_height; y++)
          {
            m_tg(index,x,y) = needle(x,y);
          }   
        }
      }
      Util::Message::endBlock();
    }
  }


  // 
  // compute color gradient
  //
  void Pb::computeCG()
  {
    assert(m_useColor);

    // allocate
    m_cg.resize(config.cgNumScales*config.cgNumOrient,m_width,m_height);
    Util::ImageStack cgA;
    Util::ImageStack cgB;

    // compute CG at each scale
    for(int cgScale = 0; cgScale < config.cgNumScales; cgScale++)
    {
      const float cgScaleVal = m_idiag * config.cgStartScale * pow(config.cgScaleFactor, cgScale);
      Group::computeCG(*m_lab.slice(Util::LAB_A), config.cgBinsA, cgScaleVal, config.cgNumOrient,
                config.cgSigma, config.cgSupport, config.cgZoom, cgA);
      Group::computeCG(*m_lab.slice(Util::LAB_B), config.cgBinsB, cgScaleVal, config.cgNumOrient,
                config.cgSigma, config.cgSupport, config.cgZoom, cgB);

      for(int orient = 0; orient < config.cgNumOrient; orient++)
      {
        const int index = cgScale * config.cgNumOrient + orient;
        for(int x = 0; x < m_width; x++)
        {
          for(int y = 0; y < m_height; y++)
          {
            m_cg(index,x,y) = 0.5 * (cgA(orient,x,y) + cgB(orient,x,y));
          }
        }
      }
    }

  }


  //
  // compute brightness gradient
  //
  void Pb::computeBG()
  {
    // allocate
    m_bg.resize(config.bgNumScales * config.bgNumOrient,m_width,m_height);
    Util::ImageStack bg;

    // compute BG at each scale
    for(int bgScale = 0; bgScale < config.bgNumScales; bgScale++)
    {
      const float bgScaleVal = m_idiag * config.bgStartScale * pow(config.bgScaleFactor, bgScale);
      Group::computeCG(*m_lab.slice(Util::LAB_L),config.bgBinsL, bgScaleVal, config.bgNumOrient,
                config.bgSigma, config.bgSupport, config.bgZoom, bg);

      for(int orient = 0; orient < config.bgNumOrient; orient++)
      {
        const int index = bgScale * config.bgNumOrient + orient;
        for(int x = 0; x < m_width; x++)
        {
          for(int y = 0; y < m_height; y++)
          {
            m_bg(index,x,y) = bg(orient,x,y);
          }
        }
      }
    }
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // helper functions for computing features and pb
  //

  float Pb::tgFitRadius(const int scale) const
  {
    const float alpha = 1.0;
    const float beta = 1.0;
    assert(scale >= 0 && scale <= config.tgNumScales);

    // compute disc radius used to compute TG at this scale 
    const float rDisc = config.tgStartScale
      * pow(config.tgScaleFactor, scale);

    // get the max sigma over filterbank used for TG
    const float sigmaFilt = config.tgFiltStartScale
      * pow(config.tgFiltScaleFactor, config.tgFiltNumScales);

    // return radius for fit, in pixels
    const float r = m_idiag * (alpha * sigmaFilt + beta * rDisc);
    return Util::max(2.1f, r);
  }

  float Pb::cgFitRadius(const int scale) const
  {
    const float r = m_idiag * config.cgStartScale
      * pow(config.cgScaleFactor, scale);
     return Util::max(2.1f, r);
  }

  float Pb::bgFitRadius(const int scale) const
  {
    const float r = m_idiag * config.bgStartScale
      * pow(config.bgScaleFactor, scale);
     return Util::max(2.1f, r);
  }

}  //namespace Group


