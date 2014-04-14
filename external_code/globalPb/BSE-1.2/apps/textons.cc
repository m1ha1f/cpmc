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
static const char* textonMapFileName;
static const char* universalTextonFile;

static int numOrient, numScales;
static float scaleFactor, startScale;
                                                                                                                                        

static void
registerConfig ()
{
    Util::registerConfig ();
    Configure::registerString (
        "image", NULL,
        "Input image file.");
    Configure::registerString (
        "textonfile", NULL,
        "Texton map file name.");
    Configure::registerString (
        "universaltextons", NULL,
        "Universal texton file.");

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
}

static void
init (const int index, const int argc, const char** argv)
{
    Util::init();

    imageFileName = Configure::getString ("image");
    textonMapFileName = Configure::getString ("textonfile");
    universalTextonFile = Configure::getString ("universaltextons");

    numOrient = Configure::getInt ("numOrient");
    numScales = Configure::getInt ("numScales");
    scaleFactor = Configure::getFloat ("scaleFactor");
    startScale = Configure::getFloat ("startScale");
                                                                                                                                        

    // Check configuration.
    if (imageFileName == NULL) { 
        throw Util::Exception ("No image file specified."); 
    }
}

static int
_main ()
{
    Util::Message::startBlock("textons_main",1);
    Util::Message::debug(Util::String("reading image %s",imageFileName),0);
    Util::ImageStack image;
    if (!Util::readJpegFile(imageFileName,image))
    {
      throw Util::Exception (Util::String ("Could not read JPEG file '%s'.", imageFileName));
    }

    const int width = image.size(1);
    const int height = image.size(2);
    Util::Message::debug(Util::String("image size [%d x %d]",width,height),1);

    Util::Image luminance;

    if (image.size(0) == 3)
    {
      Util::ImageStack lab;
      rgb2lab(image,lab);
      labNormalize(lab);
      luminance = *lab.slice(Util::LAB_L);
    }
    else
    {
      luminance = *image.slice(0);
    }

    const float idiag = sqrt (width*width + height*height);
    const float filtStartScale = idiag * startScale;
    Util::Message::debug(Util::String("Filterbank %d %d %f %f",
      numOrient,numScales,filtStartScale,scaleFactor),2);
    Group::FilterBank fbank(numOrient, numScales, filtStartScale, scaleFactor);

    Group::TextonMap textons; 
    computeTextons(luminance,fbank,universalTextonFile,textons); 

    if (textonMapFileName != NULL) 
    {
      Util::Message::debug(Util::String("writing texton map to %s",textonMapFileName),0);

      std::ofstream strm(textonMapFileName);
      if (strm.good())
      {
          strm << textons.map;
          strm.close();
      }
      else
      {
        throw Util::Exception(Util::String("Error opening %s for writing",textonMapFileName));    
      }
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
        status = 1;
    }
    return status;
}
