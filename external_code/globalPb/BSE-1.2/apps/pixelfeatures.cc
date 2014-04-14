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
#include <fstream>
#include <iostream>

#include "configure.hh"
#include "string.hh"
#include "exception.hh"
#include "util.hh"
#include "image.hh"

#include "pb.hh"


static const char* imageFileName;
static const char* outputFileName;
static const char* pbFileName;
static const char* latticeFileName;
static const char* orientedPbFileName;
static const char* orientedNMSPbFileName;
// TODO: static const char* featureFileName;
static const char* textonMapFileName;

static bool useColor;
static int numOrient;

static void
registerConfig ()
{
    Util::registerConfig ();
    Group::Pb::registerConfig();
    Configure::registerString (
        "image", NULL,
        "Input image file.");
    Configure::registerString (
        "outfile", NULL,
        "Classifier features file name.");
    Configure::registerString (
        "latticefile", NULL,
        "Output duallattice file name.");
    Configure::registerString (
        "pbfile", NULL,
        "Output pb file name.");
    Configure::registerString (
        "opbfile", NULL,
        "Output oriented pb file name.");
    Configure::registerString (
        "onmspbfile", NULL,
        "Output oriented nonmax suppressed pb file name.");
    Configure::registerString (
        "textonfile", NULL,
        "Texton map file name.");
/*  feature names not yet implemented
    Configure::registerString (
        "features", NULL,
        "Classifier feature names list output file.");
*/
    Configure::registerBool (
        "usecolor", true, 
        "Use color cue in pb?");
    Configure::registerInt (
        "norient", 8, 
        "Number of orientations at which to sample.");
}

static void
init (const int index, const int argc, const char** argv)
{
    Util::init();

    imageFileName = Configure::getString ("image");
    outputFileName = Configure::getString ("outfile");
    orientedPbFileName = Configure::getString ("opbfile");
    orientedNMSPbFileName = Configure::getString ("onmspbfile");
    pbFileName = Configure::getString ("pbfile");
    latticeFileName = Configure::getString ("latticefile");
//    featureFileName = Configure::getString ("features");
    numOrient = Configure::getInt ("norient");
    textonMapFileName = Configure::getString ("textonfile");
    useColor = Configure::getBool("usecolor");

    // Check configuration.
    if (imageFileName == NULL) { 
        throw Util::Exception ("No image file specified."); 
    }
}

static int
_main ()
{
    Util::Message::startBlock("pixelfeatures_main",1);
    Util::Message::debug(Util::String("reading image %s",imageFileName),0);
    Util::ImageStack image;
    if (!Util::readJpegFile(imageFileName,image))
    {
      throw Util::Exception (Util::String ("Could not read JPEG file '%s'.", imageFileName));
    }

    const int width = image.size(1);
    const int height = image.size(2);
    Util::Message::debug(Util::String("image size [%d x %d]",width,height),1);
 
    Group::Pb pbFactory;
    pbFactory.initialize(image,useColor);


    if (textonMapFileName != NULL) 
    {
      Util::Message::debug(Util::String("writing texton map to %s",textonMapFileName),0);

      std::ofstream strm(textonMapFileName);
      if (strm.good())
      {
          strm << pbFactory.getTextons(0).map;
          strm.close();
      }
      else
      {
        throw Util::Exception(Util::String("Error opening %s for writing",textonMapFileName));    
      }
    }

    // write out pixel feature data
    if (outputFileName != NULL) 
    {
      Util::Message::debug(Util::String("writing features to %s",outputFileName),0);
      Util::Message::debug(Util::String("feature size [3 x %d x %d]",1,numOrient),1);

      std::ofstream strm(outputFileName);
      if (strm.good())
      {
        const int scale = 0;
        for (int orient = 0; orient < numOrient; orient++) 
        {
          strm << pbFactory.getBG(scale,orient);
          strm << pbFactory.getCG(scale,orient);
          strm << pbFactory.getTG(scale,orient);
        }
        strm.close();
      }
      else
      {
        throw Util::Exception(Util::String("Error opening %s for writing",outputFileName));    
      }
    }

    //write out oriented pb data if desired
    if (orientedPbFileName != NULL) 
    {
      //now extract features and estimate pb
      Util::ImageStack pbs;
      pbFactory.computeOrientedPb(numOrient,pbs);

      Util::Message::debug(Util::String("writing oriented pb to %s",orientedPbFileName),0);
      Util::Message::debug(Util::String("opb size [%d x %d x %d]",numOrient,width,height),1);

      std::ofstream strm(orientedPbFileName);
      if (strm.good())
      {
        strm << pbs;
        strm.close();
      }
      else
      {
        throw Util::Exception(Util::String("Error opening %s for writing",orientedPbFileName));    
      }
    }

    //write out oriented pb data if desired
    if (orientedNMSPbFileName != NULL) 
    {
      //now extract features and estimate pb
      Util::ImageStack pbs;
      pbFactory.computeOrientedNMSPb(numOrient,pbs);

      Util::Message::debug(Util::String("writing oriented pb to %s",orientedNMSPbFileName),0);
      Util::Message::debug(Util::String("opb size [%d x %d x %d]",numOrient,width,height),1);

      std::ofstream strm(orientedNMSPbFileName);
      if (strm.good())
      {
        strm << pbs;
        strm.close();
      }
      else
      {
        throw Util::Exception(Util::String("Error opening %s for writing",orientedNMSPbFileName));    
      }
    }

    // write out non-max suppressed pb data if desired 
    if (pbFileName != NULL) 
    {
      Util::Message::debug(Util::String("writing pb to %s",pbFileName),0);
      Util::Message::debug(Util::String("pb size [%d x %d]",width,height),1);
      
      Util::Image boundaries;
      pbFactory.computePb(numOrient,boundaries);

      std::ofstream strm(pbFileName);
      if (strm.good())
      {
        strm << boundaries;
        strm.close();
      }
      else
      {
        throw Util::Exception(Util::String("Error opening %s for writing",pbFileName));    
      }
    }

    // write out non-max suppressed pb data if desired 
    if (latticeFileName != NULL) 
    {
      Util::Message::debug(Util::String("writing lattice pb to %s",pbFileName),0);
      Util::Message::debug(Util::String("pb size [%d x %d]",width,height),1);
      
      Group::DualLattice boundaries;
      pbFactory.computePb(numOrient,boundaries);

      std::ofstream strm(latticeFileName);
      if (strm.good())
      {
        strm << boundaries.H;
        strm << boundaries.V;
        strm.close();
      }
      else
      {
        throw Util::Exception(Util::String("Error opening %s for writing",latticeFileName));    
      }
    }

    Util::Message::endBlock(1);

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
        status = 1;
    }
    return status;
}
