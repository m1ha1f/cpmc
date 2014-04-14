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
#include <assert.h>
#include "util.hh"
#include "string.hh"
#include "message.hh"
#include "image.hh"
#include "array.hh"
#include "filterbank.hh"

namespace Group
{

  ////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // some static functions useful for building filters
  //

  //
  //return a new image that is zero mean
  //
  static void makeZeroMean(Util::Image& im)
  {
    int width = im.size(0);
    int height = im.size(1);

    float sum = 0;
    for (int i = 0; i < width; i++)
    {
      for (int j = 0; j < height; j++)
      {
        sum += im(i,j);
      }
    }
                                            
    if (width*height > 0)
    {
      float mean = (sum / ((float)width*(float)height));
      im -= mean;
    }
  }

  //
  //return an image that is unit L1norm
  //
  static void makeL1Normalized(Util::Image& im)
  {
    float norm = 0;
    int width = im.size(0);
    int height = im.size(1);
    for (int i = 0; i < width; i++)
    {
      for (int j = 0; j < height; j++)
      {
          norm += fabs(im(i,j));
      }
    }
    if (norm > 0)
    {
      im /= norm;
    }
  }


  //
  // compute 1 dimensional DFT for odd length signal
  //
  static void DFT(const Util::Array1D<float> re, const Util::Array1D<float> im,
            Util::Array1D<float>& dft_re, Util::Array1D<float>& dft_im)
  {
    int N = re.size(0);
    assert((N % 2) == 1);
    int halfN = (int)(N/2);
    float Ninv = 1/(float)N;

    dft_re.resize(N);
    dft_im.resize(N);
    dft_re.init(0);
    dft_im.init(0);

    for(int k = -halfN; k <= halfN; k++)
    {
      float d = -k * 2 * M_PI * Ninv;
      for(int j = 0; j < N; j++)
      {
        float c = cos(d * (j+1));
        float s = sin(d * (j+1));
        dft_re(k + halfN) += (re(j)*c - im(j)*s);
        dft_im(k + halfN) += (re(j)*s + im(j)*c);
      }
    }
    dft_re *= Ninv;
    dft_im *= Ninv;
  }

  //
  // compute inverse of 1 dimensional DFT for odd length signal
  //
  static void IDFT(const Util::Array1D<float> re, const Util::Array1D<float> im,
            Util::Array1D<float>& dft_re, Util::Array1D<float>& dft_im)
  {
    int N = re.size(0);
    assert((N % 2) == 1);
    int halfN = (int)(N/2);
    float Ninv = 1/(float)N;

    dft_re.resize(N);
    dft_im.resize(N);
    dft_re.init(0);
    dft_im.init(0);

    for(int j = 0; j < N; j++)
    {
      float d = (j+1) * 2.0 * M_PI * Ninv;
      for(int k = -halfN; k <= halfN; k++)
      {
        float c = cos(d * k);
        float s = sin(d * k);
        dft_re(j) += (re(k+halfN)*c - im(k+halfN)*s);
        dft_im(j) += (re(k+halfN)*s + im(k+halfN)*c);
      }
    }
  }

  //
  // compute the hilbert transform of a discrete sequence of odd length. 
  // operates in place on yval
  //
  static void computeHilbert (Util::Array1D<float>& yval)
  {
      int N = yval.size(0);
      assert((N % 2) == 1);

      Util::Array1D<float> re = yval;
      Util::Array1D<float> im(N);
      im.init(0);

      Util::Array1D<float> dft_re(N);
      Util::Array1D<float> dft_im(N);
      
      DFT(re,im,dft_re,dft_im);

      //double positive frequencies
      //zero out negative frequencies
      //leave the dc component
      for (int i = 0; i < (N/2); i++)
      {
        dft_re(i) *= 2.0;    
        dft_im(i) *= 2.0;    
      }
      for (int i = (N/2)+1; i < N; i++)
      {
        dft_re(i) = 0.0;    
        dft_im(i) = 0.0;    
      }

      IDFT(dft_re,dft_im,im,yval);
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // create an empty filterbank
  //
  FilterBank::FilterBank ()
  {
      Util::Message::debug("creating empty filterbank");
      m_nfilters = 0;
      m_filters.resize(m_nfilters);
      m_sigmasX.resize(m_nfilters);
      m_sigmasY.resize(m_nfilters);
      m_orients.resize(m_nfilters);
  }


  //
  // create a filterbank
  //
  FilterBank::FilterBank(const int numOrientations, const int numScales,
                          const float startScale, const float scaleFactor)
  {
      m_nscales = numScales;
      m_norient = numOrientations;
      m_nfilters = (2*numOrientations+1) * numScales;
      m_filters.resize(m_nfilters);
      m_sigmasX.resize(m_nfilters);
      m_sigmasY.resize(m_nfilters);
      m_orients.resize(m_nfilters);
      int filterIndex = 0;

      float sigmaRatio = 3.0;
      float sigma = startScale;
      for (int scale = 0; scale < numScales; scale++)
      {
        float sigmaX = sigma;
        float sigmaY = sigma / sigmaRatio;
        Util::Message::debug(Util::String("sx = %f, sy = %f",sigmaX,sigmaY),1);
        Util::Message::startBlock(m_norient,"filterbank construction");

        //create oriented derivative of gaussian filter pairs
        for (int i = 0; i < m_norient; i++)
        {
            Util::Message::stepBlock();
            float orientation = ((i * M_PI) / (float) m_norient);

            createGaussianKernel(sigmaX, sigmaY, 3, orientation, 2, true, m_filters(filterIndex));
            m_sigmasX(filterIndex) = sigmaX;
            m_sigmasY(filterIndex) = sigmaY;
            m_orients(filterIndex) = orientation;
            filterIndex++;

            createGaussianKernel(sigmaX, sigmaY, 3, orientation, 2, false, m_filters(filterIndex));
            m_sigmasX(filterIndex) = sigmaX;
            m_sigmasY(filterIndex) = sigmaY;
            m_orients(filterIndex) = orientation;
            filterIndex++;
        }
        Util::Message::endBlock();

        //create the center-surround filter
        //make sure it's zero-mean and L1 normalized
        float support = 4;
        float sigmaC = sigma / sigmaRatio;
        float sigmaS = sigmaC * scaleFactor;
        float supportC = support * sigmaS / sigmaC;
        float supportS = support;

        Util::Image center;
        createGaussianKernel(sigmaC, sigmaC, supportC, 0, 0, false,center);
        Util::Image surround;
        createGaussianKernel (sigmaS, sigmaS, supportS, 0, 0, false,surround);
        int hC = center.size(1);
        int wC = center.size(0);
        int hS = surround.size(1);
        int wS = surround.size(0);
        assert(hC == hS);
        assert(wC == wS);

        m_filters(filterIndex) = surround-center;
        makeZeroMean(m_filters(filterIndex));
        makeL1Normalized(m_filters(filterIndex));

        m_sigmasX(filterIndex) = sigmaS;
        m_sigmasY(filterIndex) = sigmaS;
        m_orients(filterIndex) = 0;
        filterIndex++;

        sigma = sigma * scaleFactor;
      }
      assert (filterIndex == m_nfilters);
  }

  FilterBank::~FilterBank ()
  {
    Util::Message::debug("destroying filterbank",1);
  }


  ////////////////////////////////////////////////////////////////////////////////////////////////

  const Util::Image& FilterBank::getFilter (int i) const
  {
    assert (i >= 0 && i < m_nfilters);
    return m_filters(i);
  }

  int FilterBank::getNumFilters () const
  {
    return m_nfilters;
  }

  float
  FilterBank::getSigmaX (int index) const
  {
    assert (index >= 0 && index < m_nfilters);
    return m_sigmasX(index);
  }

  float
  FilterBank::getSigmaY (int index) const
  {
    assert (index >= 0 && index < m_nfilters);
    return m_sigmasY(index);
  }

  float
  FilterBank::getOrient (int index) const
  {
    assert (index >= 0 && index < m_nfilters);
    return m_orients(index);
  }

  int FilterBank::getNscales () const
  {
    return m_nscales;
  }
                                                                                                       
  int FilterBank::getNorient () const
  {
    return m_norient;
  }
                                                                                                       
  int FilterBank::getIndexElong (int scale, int orient, int evenodd) const
  {
    return calcIndexElong(m_nscales, m_norient, scale, orient, evenodd);
  }
                                                                                                       
  int FilterBank::getIndexCenter (int scale) const
  {
    return calcIndexCenter(m_nscales, m_norient, scale);
  }


  int FilterBank::calcNumFilters (const int nscales, const int norient)
  {
    return nscales * (2 * norient + 1);
  }

  int FilterBank::calcIndexElong (const int nscales, const int norient,
                              const int scale, const int orient,
                              const int evenodd)
  {
    assert (scale >= 0 && scale < nscales);
    assert (orient >= 0 && orient < norient);
    assert (evenodd == 0 || evenodd == 1);
    return scale * (2 * norient + 1) + 2 * orient + evenodd;
  }

  int FilterBank::calcIndexCenter (const int nscales, const int norient,
                               const int scale)
  {
    assert (scale >= 0 && scale < nscales);
    return scale * (2 * norient + 1) + 2 * norient;
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////

  //
  //run the filterbank on an image and return
  //an array of images containing numFilters elements.
  //
  void FilterBank::filter (const Util::Image& im, Util::ImageStack& filtered) const
  {
    filtered.resize(m_nfilters,im.size(0),im.size(1));
    filtered.init(0);
    Util::Image f;
    Util::Message::startBlock(m_nfilters,"filtering",1);
    for (int i = 0; i < m_nfilters; i++)
    {
      Util::Message::stepBlock(1);
      getFiltered(im,getFilter(i),f);
      assert(f.size(0) == filtered.size(1));
      assert(f.size(1) == filtered.size(2));
      for (int x = 0; x < f.size(0); x++)
      {
        for (int y = 0; y < f.size(1); y++)
        {
          filtered(i,x,y) = f(x,y);
        } 
      }
    }
    Util::Message::endBlock(1);
  }

  // given a set of quadrature pair images, returns 
  // half as many new images which give the orientation
  // energy at each point
  void FilterBank::orientationEnergy (const Util::ImageStack& filtered, Util::ImageStack& energy) const
  {
    int width = filtered.size(1);
    int height = filtered.size(2);
    energy.resize(m_nscales*m_norient,width,height);
    for (int scale = 0; scale < m_nscales; scale++)
    {
      for (int orient = 0; orient < m_norient; orient++)
      {
        for (int x = 0; x < width; x++)
        {
          for (int y = 0; y < height; y++)
          {

            const float e1 = filtered(getIndexElong (scale, orient,0), x, y);
            const float e2 = filtered(getIndexElong (scale, orient,1), x, y);
            energy(scale*m_norient + orient,x,y) = sqrt(e1*e1 + e2*e2);
          }
        }
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // create a 2D gaussian kernel with sigmaX and sigmaY where
  //  orient - the orientation in radians.
  //  deriv - 2 or 0 depending on whether we want an edge filter
  //          or just a gaussian bump
  //  hilbert - whether to take the hilbert transform in the Y direction
  //  
  void FilterBank::createGaussianKernel (const float sigmaX, const float sigmaY,
                                    const float support, const float orient,
                                    const int deriv, const bool hilbert, 
                                    Util::Image& kernelImage)
  {
      assert (deriv >= 0 && deriv <= 2);

      //make an odd size kernel
      int kernelSize = (int) rint (2.0 * support * Util::max (sigmaX, sigmaY));
      kernelSize = kernelSize + (1 - (kernelSize % 2));
      assert ((kernelSize % 2) == 1);
      int halfSize = (kernelSize - 1) / 2;

      //create a grayscale image
      kernelImage.resize(kernelSize, kernelSize);

      //compute the kernel value along the x and y directions at a fine
      //resolution
      int maxSize = 2 * kernelSize;
      int zoom = 64;
      int nsamples = zoom * maxSize + 1;
      int halfnsamples = zoom * maxSize / 2;
      Util::Array1D<float> xval(nsamples);
      Util::Array1D<float> yval(nsamples);
      for (int y = -halfnsamples; y <= halfnsamples; y++)
      {
        const float sigmaYsqr = (sigmaY * sigmaY);
        const float sigmaXsqr = (sigmaX * sigmaX);
        const float u = (((float) y) / ((float) zoom));
        xval(y + halfnsamples) = exp (-u * u / (2 * sigmaXsqr));
        if (deriv == 0)
        {
            yval(y + halfnsamples) = exp (-u * u / (2 * sigmaYsqr));
        }
        else if (deriv == 1)
        {
            yval(y + halfnsamples) = exp (-u * u / (2 * sigmaYsqr)) * (-u / sigmaYsqr);
        }
        else if (deriv == 2)
        {
            yval(y + halfnsamples) = exp (-u * u / (2 * sigmaYsqr)) * ((u * u / sigmaYsqr) - 1);
        }
      }

      //take the hilbert transform in the Y direction if necessary
      if (hilbert == true)
      {
        computeHilbert (yval);
      }

      //now create the filter by averaging values out of the sampled
      //kernels.
      for (int x = -halfSize; x <= halfSize; x++)
      {
        for (int y = -halfSize; y <= halfSize; y++)
        {
          float u = x * cos(orient) + y * sin(orient);
          float v = x * sin(orient) - y * cos(orient);

          // compute the x value
          int u1 = (int) rint (((u - 0.5) * zoom) + halfnsamples);
          int u2 = (int) rint (((u + 0.5) * zoom) + halfnsamples);
          u1 = Util::max (u1, 0);
          u2 = Util::min (u2, nsamples - 1);
          assert (u2 >= u1);
          float xvalue = 0;
          for (int i = u1; i <= u2; i++)
          {
              xvalue += xval(i);
          }
          xvalue /= (u2 - u1 + 1);

          // compute the y value
          int v1 = (int) rint (((v - 0.5) * zoom) + halfnsamples);
          int v2 = (int) rint (((v + 0.5) * zoom) + halfnsamples);
          v1 = Util::max (v1, 0);
          v2 = Util::min (v2, nsamples - 1);
          assert (v2 >= v1);
          float yvalue = 0;
          for (int i = v1; i <= v2; i++)
          {
              yvalue += yval(i);
          }
          yvalue /= (v2 - v1 + 1);

          float val = xvalue * yvalue;
          kernelImage(x + halfSize, y + halfSize) = val;
        }
      }

      //make it zero mean if it's a derivative of a gaussian
      if (deriv != 0)
      {
        makeZeroMean(kernelImage);
      }

      //make sure the filter has unit L1 norm 
      makeL1Normalized(kernelImage);
  }


  /////////////////////////////////////////////////////////////////////////////////////////////////////

  //
  // paste kernel into a bigger image
  //
  static void drawFilter(const Util::Image& filt, int cx, int cy, Util::Image& im)
  {
    int iwidth = im.size(0);
    int iheight = im.size(1);
    int w = filt.size(0);
    int h = filt.size(1);
    // (i,j) are coorindates local to the filter
    for (int i = 0; i < w; i++)
    {
      for (int j = 0; j < h; j++)
      {
        // (x,y) are image coordinates for (i,j)
        int x = cx + i - w / 2;
        int y = cy + j - h / 2;
        assert (x >= 0 && x < iwidth);
        assert (y >= 0 && y < iheight);
        im(x,y) = filt(i,j);
      }
    }
  }

  //
  // create a nice image of the filterbank for presentation or debugging
  //
  void FilterBank::getFilterBankImage (Util::Image& fbim) const
  {
      // figure out the max filter dimensions
      int maxd = 0;
      for (int i = 0; i < m_nfilters; i++)
      {
        maxd = Util::max (maxd,getFilter(i).size(0));
        maxd = Util::max (maxd,getFilter(i).size(1));
      }
      maxd++;

      // allocate an image large enough to hold images of all the filters
      int nscales = getNscales ();
      int norient = getNorient ();
      int iheight = maxd * nscales;
      int iwidth = maxd * (2 * norient + 1);
      fbim.resize(iwidth,iheight);
      fbim.init(0);

      // write the filters into the big image
      // (cx,cy) give the image coordinates for the center of
      // the this filter
      for (int scale = 0; scale < nscales; scale++)
      {
        // lay out filters for this scale like this: eoeoeoeoeoeoc
        for (int orient = 0; orient < norient; orient++)
        {
          for (int evenodd = 0; evenodd < 2; evenodd++)
          {
            int cx = (2 * orient + evenodd) * maxd + maxd / 2;
            int cy = scale * maxd + maxd / 2;
            const Util::Image filt = getFilter(getIndexElong (scale, orient, evenodd));
            drawFilter (filt, cx, cy, fbim);
          }
        }
        int cx = 2 * norient * maxd + maxd / 2;
        int cy = scale * maxd + maxd / 2;
        const Util::Image filt = getFilter(getIndexCenter (scale));
        drawFilter (filt, cx, cy, fbim);
      }
  }

} //namespace Group
