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
#include <time.h>
#include <math.h>
#include <float.h>
#include <vector>
#include <fstream>

#include "util.hh"
#include "exception.hh"
#include "string.hh"
#include "configure.hh"
#include "random.hh"
#include "array.hh"
#include "image.hh"
#include "filterbank.hh"
#include "kmeans.hh"

// Compute textons from a list of images.

// input images
static const char** imageFileNames;
static int numImageFiles;

// filter bank
static int numOrient, numScales;
static float scaleFactor, startScale;

// kmeans
static int k, iter;
static float tol;

// output
static const char* outFileName;

// misc
static float sampleRate;

// globals
static int numSamples, dim;
static Util::Array2D<float> samples;
static Util::Array2D<float> textons;

static void
registerConfig ()
{
    Util::registerConfig ();

    // filter bank
    Configure::registerInt (
        "numOrient", 6, 
        "Number of orientations for filterbank.");
    Configure::registerInt (
        "numScales", 1, 
        "Number of scales for filterbank.");
    Configure::registerFloat (
        "scaleFactor", sqrt(2), 
        "Scale factor for filterbank.");
    Configure::registerFloat (
        "startScale", 0.007, 
        "Beginning scale for filterbank, as fraction of image diagonal.");

    // kmeans
    Configure::registerInt (
        "k", 128, 
        "Number of textons.");
    Configure::registerFloat (
        "tol", 1e-5,
        "Kmeans stopping tolerance (fractional change in RMS error).");
    Configure::registerInt (
        "iter", 200,
        "Kmeans max iterations.");

    // output
    Configure::registerString (
        "out", NULL,
        "Output file for textons.");

    // misc
    Configure::registerFloat (
        "sampleRate", 1,
        "Fraction of pixels to use in the computation.");
}

static void
init (const int index, const int argc, const char** argv)
{
    Util::init ();
    
    // input images
    imageFileNames = &argv[index];
    numImageFiles = argc - index;

    // filter bank
    numOrient = Configure::getInt ("numOrient");
    numScales = Configure::getInt ("numScales");
    scaleFactor = Configure::getFloat ("scaleFactor");
    startScale = Configure::getFloat ("startScale");

    // kmeans
    k = Configure::getInt ("k");
    tol = Configure::getFloat ("tol");
    iter = Configure::getInt ("iter");

    // output
    outFileName = Configure::getString ("out");

    // misc
    sampleRate = Configure::getFloat ("sampleRate");

    // check arguments
    if (numImageFiles == 0) {
        throw Util::Exception ("No image files specified."); 
    }
    if (sampleRate <= 0 || sampleRate > 1) {
        throw Util::Exception (
            "Argument sampleRate out of range (0,1].");
    }
    if (k <= 0) {
        throw Util::Exception (
            "Argument k out of range (0,inf].");
    }
}

static void
readImages ()
{
    Util::Message::debug(Util::String("Number of images = %d\n",numImageFiles));
    int width = 0;
    int height = 0;

    // Read in the images to compute how many samples we'll need.
    Util::Array1D<int> sampPerImage (numImageFiles);
    Util::Message::startBlock(numImageFiles, "counting samples");
    for (int i = 0; i < numImageFiles; i++) 
    {
      Util::Message::stepBlock();

      Util::ImageStack image;
      if (!Util::readJpegFile(imageFileNames[i],image)) {
        throw Util::Exception (Util::String (
          "Could not read JPEG file '%s'.", imageFileNames[i]));
      }

      const int w = image.size(1);
      const int h = image.size(2);
      const int m = Util::minmax(0,(int)rint(sampleRate*(w*h-1)),w*h-1);
      sampPerImage(i) = m;
      numSamples += m;
      if (i == 0) 
      {
        width = w; 
        height = h;
      } 
      else 
      {
        /*
        if (w*h != width*height) 
        {
          throw Util::Exception ("Images are not all the same size");
        }
        */
      }
    }
    Util::Message::endBlock();
    Util::Message::debug(Util::String("Number of samples = %d",numSamples));
    if (numSamples <= 0) 
    {
      throw Util::Exception (Util::String ("Can't run with %d samples!", numSamples));
    }

    // Build filter bank.
    const float idiag = sqrt (width*width + height*height);
    const float filtStartScale = idiag * startScale;
    Util::Message::debug(Util::String("Filterbank %d %d %f %f",
      numOrient,numScales,filtStartScale,scaleFactor),2);
    Group::FilterBank fbank(numOrient, numScales, filtStartScale, scaleFactor);
    dim = fbank.getNumFilters();

    // Allocate space for the samples.
    std::cerr << "Allocating space for samples..." << std::endl;
    samples.resize(numSamples,dim);

    // Filter all the images and extract the samples.
    Util::Message::startBlock(numImageFiles, "filtering");
    int count = 0;
    for (int i = 0; i < numImageFiles; i++) 
    {
      Util::Message::stepBlock();

      // Read in the image.
      Util::ImageStack imageRGB;
      if (!Util::readJpegFile(imageFileNames[i],imageRGB)) {
        throw Util::Exception (Util::String (
          "Could not read JPEG file '%s'.", imageFileNames[i]));
      }

//      Util::ImageStack imageLAB;
//      rgb2lab(imageRGB,imageLAB); 
//      labNormalize(imageLAB);
//      Util::Image image = *imageLAB.slice(Util::LAB_L);
      Util::Image image = *imageRGB.slice(Util::RGB_R);

      const int w = image.size(0);
      const int h = image.size(1);
      const int m = sampPerImage(i);

      // Filter the image.
      //assert (w*h == width*height);
      Util::ImageStack filtered; 
      fbank.filter(image,filtered);

      // Sample the filter responses.
      Util::Array1D<uint> which(m);
      Util::kOfN (m, w*h, which.data());
      for (int i = 0; i < m; i++) 
      {
          const int j = which(i);
          const int x = j % w;
          const int y = j / w;
          
          assert (j >= 0 && j < w*h);
          assert (x >= 0 && x < w);
          assert (y >= 0 && y < h);

          assert (count < numSamples);
          for (int d = 0; d < dim; d++) {
            samples(count,d) = filtered(d,x,y);
          }
          count++;
      }
    }
    Util::Message::endBlock();
}

static void
computeTextons ()
{
    Util::Message::debug("computing textons");
    Util::Array1D<int> membership (numSamples);
    textons.resize (k,dim);
    const bool multilevel = true;
    kmeans(samples, k, tol, iter, multilevel, textons, membership);
}

static void
writeTextons ()
{
    if (outFileName == NULL) {
        Util::Message::debug("no texton output file given");
        return;
    }

    Util::Message::debug(Util::String("writing textons to: %s",outFileName));

    std::ofstream os(outFileName);
    if (!os.good()) 
    {
      throw Util::Exception (Util::String ("Error opening '%s' for writing.", outFileName));
    }
    os << textons;
    os.close();
}

static int
_main ()
{
    readImages();
    computeTextons();
    writeTextons();
    return 0;
}

                                                                                                                 
//////////////////////////////////////////////////////////////////////
//////////// Generic program wrapper below this point ////////////////
//////////////////////////////////////////////////////////////////////
                                                                                                                 
// Print some useful info to stderr.
static void
infoBegin (const int argc, const char** argv)
{
    static char hostname[256];
    if (gethostname (hostname, 256) == 0) {
        Util::Message::debug(Util::String("host: %s",hostname),1);
    }
                                                                                                                 
    static char cwd[1000];
    if (getcwd (cwd, 1000) != NULL) {
        Util::Message::debug(Util::String("pwd: %s",cwd),1);
    }
                                                                                                                 
    Util::String args ("command: %s", argv[0]);
    for (int i = 1; i < argc; i++) {
        args.append (" %s", argv[i]);
    }
    Util::Message::debug(args,1);
                                                                                                                 
    time_t t = time (NULL);
    Util::Message::debug(Util::String("start time: %s",ctime(&t)),1);
}
                                                                                                                 
int
main (const int argc, const char** argv)
{
    int status = 0;
    try
    {
        // Register configuration options.
        registerConfig();
                                                                                                                 
        // Initialize all modules.
        try
        {
            int index = Configure::init (argc, (const char**) argv);
            init (index, argc, argv);
            infoBegin (argc, argv);
            Configure::show();
        }
        catch (Util::Exception& e)
        {
            Util::Message::debug(Util::String("usage: %s  <image1.jpg> <image2.jpg> .... ",argv[0]));
            Configure::usage();
            Util::Message::error(e.msg());
            exit (1);
        }
        // Run the main program.
        status = _main();
    }
    catch (Util::Exception& e)
    {
        // Some sort of uncaught runtime error.
        Util::Message::error(e.msg());
        status = 1;
    }
    return status;
}


