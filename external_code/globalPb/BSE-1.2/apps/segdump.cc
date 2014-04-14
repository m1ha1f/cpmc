// Copyright (C) 2005 Charless C. Fowlkes <fowlkes@eecs.berkeley.edu>
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
#include <fstream>
#include <iostream>

#include "types.hh"
#include "array.hh"
#include "util.hh"
#include "exception.hh"
#include "string.hh"
#include "configure.hh"
#include "kmeans.hh"
#include "segmentation.hh"

#include "smatrix.hh"
#include "dmatrix.hh"

#include "affinity.hh"
#include "pb.hh"
#include "ic.hh"
#include "region.hh"


//read/write cached Pb edge map?
static bool readLattice;
static const char* latticeFileName;

//cues
static Group::cueModeType cueMode;
static const char* const cueModeStrings[] = {"ic","patch","both",NULL};
static bool useColor;
static int dthresh;
static float sigma;

//affinity matrix
static const char* smatrixFileName;
static const char* textonFileName;
static const char* brightnessFileName;
static const char* colorFileNameA;
static const char* colorFileNameB;

//the image
static const char* imageFileName;

static void
registerConfig ()
{
  Util::registerConfig();
  Group::Pb::registerConfig();
  Group::Region::registerConfig();

  Configure::registerString(
      "image", NULL,
      "Image file to proccess");

  Configure::registerBool(
      "readlattice", false,
      "Read in pb data?");

  Configure::registerString(
      "latticefile", NULL,
      "Pb lattice data file name");

  Configure::registerInt(
      "dthresh",3,"Connectivity distance threshold in pixels");

  Configure::registerFloat(
      "sigma",0.1,"Sigma non-linearity");

  Configure::registerBool(
      "usecolor", true,
      "Use color info.");

  Configure::registerString(
      "smatrixfile", NULL,
      "Write out affinity matrix");

  Configure::registerString(
      "textonmap", NULL,
      "Write out texton map");

  Configure::registerString(
      "brightnessmap", NULL,
      "Write out brightness histogram map");

  Configure::registerString(
      "colormapa", NULL,
      "Write out color histogram map (A channel)");

  Configure::registerString(
      "colormapb", NULL,
      "Write out color histogram map (B channel)");

}

static void
init (const int index, const int argc, const char** argv)
{
  Util::init();

  imageFileName = Configure::getString("image");

  latticeFileName = Configure::getString("latticefile"); 
  readLattice = Configure::getBool("readlattice");

  dthresh = Configure::getInt("dthresh");
  sigma = Configure::getFloat("sigma");
  useColor = Configure::getBool("usecolor");
  cueMode = Group::ic;
  smatrixFileName = Configure::getString("smatrixfile");
  textonFileName = Configure::getString("textonmap");
  brightnessFileName = Configure::getString("brightnessmap");
  colorFileNameA = Configure::getString("colormapa");
  colorFileNameB = Configure::getString("colormapb");

  if (imageFileName == NULL)
  {
    throw Util::Exception (Util::String ("image filename is required"));
  }
}

//////////////////////////////////////////////////////////////////////////////

static int
_main ()
{
    Util::Message::startBlock("segment_main");
    Util::Message::debug(Util::String("reading image %s",imageFileName));
    Util::ImageStack image;
    if (!Util::readJpegFile(imageFileName,image))
    {
      throw Util::Exception (Util::String ("Could not read JPEG file '%s'.", imageFileName));
    }

    const int width = image.size(1);
    const int height = image.size(2);
    Util::Message::debug(Util::String("image size [%d x %d]",width,height),1);

    if (smatrixFileName != NULL)
    {
      Group::DualLattice boundaries; 
      Util::Message::startBlock("acquiring contour features");
      if (readLattice)
      {
        if (latticeFileName != NULL)
        {
          Util::Message::debug(Util::String("reading latticefile %s",latticeFileName));
          std::ifstream in(latticeFileName);
          if (!in.good()) {
            throw Util::Exception (Util::String ("Error opening '%s' for reading.", latticeFileName));
          }
          in >> boundaries.H;
          in >> boundaries.V;
          boundaries.width = boundaries.H.size(0);
          boundaries.height = boundaries.V.size(1);
          in.close();
        }      
      }
      else
      {
        Group::Pb pbFactory;
        pbFactory.initialize(image,useColor);
        pbFactory.computePb(8,boundaries);
        if (latticeFileName != NULL)
        {
          Util::Message::debug("writing lattice file");
          std::ofstream out(latticeFileName);
          if (!out.good()) {
            throw Util::Exception(Util::String ("Error opening '%s' for writing.", latticeFileName));
          }
          out << boundaries.H;
          out << boundaries.V;
          out.close();
        }
      }

      Group::Region patchFactory;
      int smapRadius = dthresh;
      Util::Message::debug(Util::String("running intervening contour with map radius = %d", smapRadius));
      Group::SupportMap ic;
      Group::computeSupport(boundaries,smapRadius,1.0f,ic);
      Util::Message::endBlock();

      Util::Message::debug("building affinity matrix");
      SMatrix* affinities = NULL;
      Group::computeAffinities(ic,patchFactory,cueMode,sigma,dthresh,useColor,&affinities);
      assert(affinities != NULL);


      Util::Message::debug("writing affinity matrix");
      FILE* fp = fopen(smatrixFileName,"w");
      affinities->dump(fp);
      fclose(fp);
    }

    //////////////////////////////////////////////////////////////////////////////

    Group::Region patchFactory2;
    Util::Message::startBlock("acquiring region features");
    patchFactory2.initialize(image);
    Util::Message::endBlock();

    if (textonFileName != NULL)
    {
      Util::Message::debug(Util::String("writing texton map file %s",textonFileName));
      std::ofstream out(textonFileName);
      out << patchFactory2.getTextonHist(0);
      out.close();
    }

    if (brightnessFileName != NULL)
    {
      Util::Message::debug(Util::String("writing brightness map file %s",brightnessFileName));
      std::ofstream out(brightnessFileName);
      out << patchFactory2.getLHist(0);
      out.close();
    }

    if (colorFileNameA != NULL)
    {
      Util::Message::debug(Util::String("writing coloron A map file %s",colorFileNameA));
      std::ofstream out(colorFileNameA);
      out << patchFactory2.getAHist(0);
      out.close();
    }

    if (colorFileNameB != NULL)
    {
      Util::Message::debug(Util::String("writing coloron B map file %s",colorFileNameB));
      std::ofstream out(colorFileNameB);
      out << patchFactory2.getBHist(0);
      out.close();
    }

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
            Util::Message::debug(Util::String("usage: %s  -image <image.jpg> [options]",argv[0]));
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
        Util::Message::debug("EXITING WITH ERROR!!!!");
        status = 0; //exit with status 0 even though there was a problem so that milrun still goes
    }
    return status;
}


