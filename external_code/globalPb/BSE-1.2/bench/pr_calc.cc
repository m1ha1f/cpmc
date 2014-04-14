#include <fstream>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <float.h>
#include "image.hh"
#include "util.hh"
#include "segmentation.hh"
#include "exception.hh"
#include "string.hh"
#include "configure.hh"
#include "array.hh"
#include "match.hh"

// Produce PR data given an image, a Pb image, human segmentations,
// and a matching thresholds.

static const char* pbFileName;
static const char* egradFileName;
static const char* prFileName;
static double distTolMatch, thetaTolMatch, outlierCost;
static int evecnum;
static int nthresh;
static bool halfSize;

static const char** segFiles;
static int numSegFiles;
static Util::Image pb;
static Util::ImageStack egrad;
static int width, height;
static double* threshvals;
static const char* const threshModeStrings[] = {"linear","adaptive","telescoping",NULL};
enum threshModes {linear,adaptive,telescoping};
static threshModes threshMode;
Util::Array1D< Util::Array2D<bool> > bmaps;

// for precision-recall
static Util::Array1D<int> humanHit;
static Util::Array1D<int> humanSum;
static Util::Array1D<int> machineHit;
static Util::Array1D<int> machineSum;

static void
registerConfig ()
{
    Util::registerConfig ();

    Configure::registerString (
        "pbfile", NULL,
        "Input pb file.");
    Configure::registerString (
        "egradfile", NULL,
        "Input egrad file.");
    Configure::registerString (
        "prfile", NULL,
        "Output precision-recall data file.");
    Configure::registerFloat (
        "distTolMatch", 0.01, 
        "Tolerance for boundary matching, as fraction of image diagonal.");
    Configure::registerFloat (
        "thetaTolMatch", 45, 
        "Tolerance for segmentation matching, in degrees.");
    Configure::registerInt (
        "nthresh", 30,
        "Num points for curve.");
    Configure::registerEnum(
        "threshmode", threshModeStrings,linear,
        "How to generate threshold levels?");
    Configure::registerInt (
        "evec", 8,
        "Which eigenvector to use from the egrad file.");
    Configure::registerBool (
        "halfsize", false,
        "Use half-size segmentations.");
}

static void
init (const int index, const int argc, const char** argv)
{
    Util::init ();

    pbFileName = Configure::getString ("pbfile");
    egradFileName = Configure::getString ("egradfile");
    prFileName = Configure::getString ("prfile");
    nthresh = Configure::getInt ("nthresh");
    evecnum = Configure::getInt ("evec");
    threshMode = (threshModes)Configure::getEnum("threshmode");
    halfSize= Configure::getBool("halfsize");
    segFiles = &argv[index];
    numSegFiles = argc - index;

    distTolMatch = Configure::getFloat ("distTolMatch");
    thetaTolMatch = M_PI/180 * Configure::getFloat ("thetaTolMatch");
    outlierCost = distTolMatch * 2;

    // Check configuration.
    if ((pbFileName == NULL) && (egradFileName == NULL)) { 
        throw Util::Exception ("No pb file specified."); 
    }
    if (prFileName == NULL) { 
        throw Util::Exception ("No pr file specified."); 
    }
    if (numSegFiles == 0) {
        throw Util::Exception ("No seg files specified."); 
    }
}



static void
readPb()
{
    if (pbFileName != NULL)
    {
        assert(pbFileName != NULL);
        std::ifstream is (pbFileName);
        if (!is.good()) {
            throw Util::Exception (Util::String (
                "Error opening '%s' for reading.", pbFileName));
        }
        is >> pb;
        is.close();
    }
    else
    {
        assert(egradFileName != NULL);
        std::ifstream is (egradFileName);
        if (!is.good()) 
        {
            throw Util::Exception (Util::String (
                "Error opening '%s' for reading.", egradFileName));
        }
        is >> egrad;
        int nvec = egrad.size(0);
        width = egrad.size(1);
        height = egrad.size(2);
        if (nvec < evecnum)
        {
            throw Util::Exception (Util::String (
                "Only %d available eigengradients.", nvec));
        }
        pb.resize(width,height);
        for (int y = 0; y < height; y++) 
        {
          for (int x = 0; x < width; x++) 
          {
            pb(x,y) = egrad(evecnum,x,y);
          }
        }
    }

    width = pb.size(0);
    height = pb.size(1);

    std::cerr << width << " " << height << std::endl;

    double* pblist = new double[width*height];
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            double p = pb(x,y);
            if (p < 0)
            {
                fprintf (stderr, "WARNING: pb(%d,%d)=%g is out of range [0,inf].\n", x, y, p);
            }
            pblist[y*width+x] = Util::max(p,0.0);
        }
    }

    //even increments over the usable part of the data
    Util::sort(pblist,width*height);
    threshvals = new double[nthresh+1];
    if (threshMode == telescoping)
    {
      double maxval = pblist[width*height-1];
      std::cerr << "max val = " << maxval << std::endl;
      for (int i = 0; i <= nthresh; i++)
      {
        threshvals[i] = maxval*((double)i / (double)nthresh);
      }
    }
    else if (threshMode == adaptive)
    {
      int zoffset = 0;
      while (pblist[zoffset] == 0.0)
      {
        zoffset++;
      }
      int nzarea = width*height - zoffset;
      int stepsize = ((int)((float)nzarea/(float)nthresh)) - 1;

      for (int i = 0; i <= nthresh; i++)
      {
        assert(i*stepsize + zoffset < width*height);
        threshvals[i] = pblist[i*stepsize + zoffset];      
      }
    }
    else
    {
      for (int i = 0; i <= nthresh; i++)
      {
        threshvals[i] = (double)i / (double)nthresh;
      }
    }

    delete pblist;
    std::cerr << std::endl;
}

static void
computePR()
{
    humanHit.resize(nthresh+1);
    humanSum.resize(nthresh+1);
    machineHit.resize(nthresh+1);
    machineSum.resize(nthresh+1);

    humanHit.init(0);
    humanSum.init(0);
    machineHit.init(0);
    machineSum.init(0);

    Util::Array2D<bool> machine (width,height);

    std::cerr << nthresh << " thresholds [";
    for (int t = 0; t <= nthresh; t++) {
        std::cerr << ".";
        // compute machine boundary map from pb image for this threshold
        //const double threshold = (double) t / nthresh;
        const double threshold = threshvals[t];
        std::cerr << " " << threshvals[t];
        for (int y = 0; y < height; y++) {
            for (int x = 0; x < width; x++) {
                const double p = pb(x,y);
                machine(x,y) = (p >= threshold);
            }
        }
        // correspond the machine boundary map with each valid
        // segmentation in turn
        Pixel outlier (-1,-1);
        Util::Array2D<Pixel> matchM (width,height); // machine matches
        Util::Array2D<Pixel> matchH (width,height); // human matches
        Util::Array2D<bool> matchMacc (width,height); // accumulated machine matches
        matchMacc.init(false);
        for (int i = 0; i < numSegFiles; i++) 
        {
            Util::Array2D<bool>& human = bmaps(i);
            const double idiag = sqrt (width*width + height*height);

            (void) matchEdgeMaps (
                width, height, machine, human,
                distTolMatch*idiag, outlierCost*idiag, matchM, matchH);

            for (int x = 0; x < width; x++) {
                for (int y = 0; y < height; y++) {
                    if (machine(x,y) && matchM(x,y) != outlier) {
                        matchMacc(x,y) = true;
                    }
                    if (human(x,y)) {
                        humanSum(t) += 1;
                        humanHit(t) += (matchH(x,y) != outlier);
                    }
                }
            }
        }
        for (int x = 0; x < width; x++) {
            for (int y = 0; y < height; y++) {
                if (machine(x,y)) {
                    machineSum(t) += 1;
                    machineHit(t) += matchMacc(x,y);
                }
            }
        }
    }    
    std::cerr << "]" << std::endl;
}

static void
writePR()
{
    assert (prFileName != NULL);
    std::cerr << "writing " << Util::String("%s",prFileName) << std::endl;
    FILE* strm = Util::openOutputStrm (Util::String("%s",prFileName));
    for (int t = 0; t <= nthresh; t++) {
        fprintf (strm, "%10.3g %10d %10d %10d %10d\n", 
                 threshvals[t],
                 humanHit(t), humanSum(t),
                 machineHit(t), machineSum(t));
    }    
    fclose(strm);
    std::cerr << "done " << std::endl;
}

static int
_main ()
{
    std::cerr << "Reading Pb image..." << std::endl;
    readPb();

    bmaps.resize(numSegFiles);
    std::cerr << "Reading segmentation files..." << std::endl;
    for (int i = 0; i < numSegFiles; i++) {
      std::cerr << "Reading seg file '" << segFiles[i] << "'..." << std::endl;
      Util::Segmentation seg;
      seg.readFile(segFiles[i]);
      if (halfSize)
      {
        seg.computeBoundaryMapHalf(bmaps(i));
      }
      else
      {
        seg.computeBoundaryMap(bmaps(i));
      }
    }

    std::cerr << "Computing precision/recall..." << std::endl;
    computePR();
    std::cerr << "Writing output file " << prFileName << "..." << std::endl;
    writePR();


    std::cerr << "Done." << std::endl;
    return 0;
}

//////////////////////////////////////////////////////////////////////
//////////// Generic program wrapper below this point ////////////////
//////////////////////////////////////////////////////////////////////

static Util::Timer programTimer;

// Print some useful info to stderr.
static void
infoBegin (const int argc, const char** argv)
{
    time_t t = time (NULL);
    std::cerr << "time: " << ctime(&t);
    
    static char hostname[256];
    if (gethostname (hostname, 256) == 0) {
        std::cerr << "host: " << hostname << std::endl;
    }
    
    char* pwd = getcwd (NULL, 0); 
    if (pwd != NULL) { 
        std::cerr << "pwd: " << pwd << std::endl;
    }
    
    Util::String args ("%s", argv[0]);
    for (int i = 1; i < argc; i++) {
        args.append (" %s", argv[i]);
    }
    std::cerr << "command: " << args << std::endl;

    programTimer.reset();
    programTimer.start();
}

// Print some useful info to stderr.
static void
infoEnd (const int status)
{
    programTimer.stop();

    double e = programTimer.elapsed();
    double c = programTimer.cpu();
    double u = programTimer.user();
    double s = programTimer.system();

    std::cerr << std::endl;
    fprintf (stderr, "%10s: %11s\n", 
             "elapsed", Util::Timer::formatTime(e,1));
    fprintf (stderr, "%10s: %11s  (%.0f%% of %s)\n", 
             "cpu", Util::Timer::formatTime(c,1), c/e*100, "elapsed");
    fprintf (stderr, "%10s: %11s  (%.0f%% of %s)\n", 
             "user", Util::Timer::formatTime(u,1), u/c*100, "cpu");
    fprintf (stderr, "%10s: %11s  (%.0f%% of %s)\n", 
             "system", Util::Timer::formatTime(s,1), s/c*100, "cpu");
    std::cerr << std::endl;

    time_t t = time(NULL);
    std::cerr << "end time: " << ctime(&t);

    std::cerr << "status: " << status << std::endl;
}

int
main (const int argc, const char** argv)
{
    int status = 0;
    infoBegin (argc, argv);
    try {
        // Register configuration options.
        registerConfig();
        // Initialize all modules.
        try {
            int index = Configure::init (argc, (const char**) argv);
            init (index, argc, argv);
            // Output configuration to stderr.
            std::cerr << std::endl << "CONFIGURATION BEGIN" << std::endl;
            Configure::show ();
            std::cerr << "CONFIGURATION END" << std::endl;
        } catch (Util::Exception& e) {
            // Initialization error.  Print usage message and exit.
            std::cerr << std::endl;
            std::cerr << "usage: " << argv[0] 
                 << " [options] <seg files...>" << std::endl;
            Configure::usage();
            std::cerr << std::endl;
            std::cerr << e << std::endl;
            std::cerr << "Program halted." << std::endl;
            exit (1);
        }
        // Run the main program.
        status = _main();
    } 
    catch (Util::Exception& e) {
        // Some sort of uncaught runtime error.
        std::cerr << e << std::endl;
        std::cerr << "Program halted." << std::endl;
        status = 1;
    }
    infoEnd(status);
    return status;
}

  /*
  Matrix*
  computeUnionDMap (
      const int numSegFiles, Matrix const* const* dmaps,
      Matrix const* const* homaps, Matrix* unionHOmap)
  {
      assert (numSegFiles > 0);
      // compute the union dmap: min dist to any segmentation boundary
      const int width = dmaps[0]->getWidth();
      const int height = dmaps[0]->getHeight();
      Matrix* unionDmap = new Matrix (width,height);
      for (int y = 0; y < height; y++) {
          for (int x = 0; x < width; x++) {
              double dmin = dmaps[0]->get(y,x);
              int minSeg = 0;
              for (int seg = 0; seg < numSegFiles; seg++) {
                  const double d = dmaps[seg]->get(y,x);
                  if (d < dmin) {
                      dmin = d;
                      minSeg = seg;
                  }
              }
              unionDmap->set(y,x,dmin);
              if (unionHOmap != NULL) {
                  unionHOmap->set(y,x,homaps[minSeg]->get(y,x));
              }
          }
      }
      return unionDmap;
  }
  */
   

