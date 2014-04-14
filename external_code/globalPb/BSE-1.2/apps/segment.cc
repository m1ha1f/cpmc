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

//write out intermediate results
static const char* evecFileName;
static const char* egradFileName;

//cues
static Group::cueModeType cueMode;
static const char* const cueModeStrings[] = {"ic","patch","both",NULL};
static bool useColor;
static int dthresh;
static float sigma;

//computing discrete ncut
static const char** segFiles;
static int numSegFiles;
static int dustThresh;

//affinity matrix
static const char* smatrixFileName;

//eigenvectors
static bool filterEvec;
static int numVectors;

//clustering
static int numSuperpixels;
static const char* segFileName;
static const char* graphFileName;

//static int numGroups;
//static const char* segFileName;

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
      "dthresh",5,"Connectivity distance threshold in pixels");

  Configure::registerFloat(
      "sigma",0.1,"Sigma non-linearity");

  Configure::registerBool(
      "usecolor", true,
      "Use color info.");

  Configure::registerEnum(
     "cues",cueModeStrings,Group::ic,
     "What cues to use in similarity computation");

  Configure::registerString(
      "smatrixfile", NULL,
      "Write out affinity matrix");

  Configure::registerInt(
      "numvectors", 16,
      "Number of eigenvectors.");

  Configure::registerBool(
      "filterevec", true,
      "Median filter eigenvectors.");

  Configure::registerString(
      "evecfile", NULL,
      "Write out eigenvector data");

  Configure::registerString(
      "egradfile", NULL,
      "Write out eigenvector gradient data");

  Configure::registerInt(
      "numsuperpixels", 50,
      "Number of superpixels.");

  Configure::registerString(
      "segfile", NULL,
      "Write out superpixel segmentation");

  Configure::registerString(
      "graphfile", NULL,
      "Write out collapsed graph data");

  Configure::registerInt(
      "dustthresh", 10,
      "How small before we merge a disconnected component.");
}

static void
init (const int index, const int argc, const char** argv)
{
  Util::init();

  segFiles = &argv[index];
  numSegFiles = argc - index;

  imageFileName = Configure::getString("image");

  latticeFileName = Configure::getString("latticefile"); 
  readLattice = Configure::getBool("readlattice");

  dthresh = Configure::getInt("dthresh");
  sigma = Configure::getFloat("sigma");
  useColor = Configure::getBool("usecolor");
  cueMode = (Group::cueModeType)Configure::getEnum("cues");
  smatrixFileName = Configure::getString("smatrixfile");

  numVectors = Configure::getInt("numvectors");
  filterEvec = Configure::getBool("filterevec");
  egradFileName = Configure::getString("egradfile");
  evecFileName = Configure::getString("evecfile");

  graphFileName = Configure::getString("graphfile");
  numSuperpixels = Configure::getInt("numsuperpixels");
  dustThresh = Configure::getInt("dustthresh");
  segFileName = Configure::getString("segfile");

  if (imageFileName == NULL)
  {
    throw Util::Exception (Util::String ("image filename is required"));
  }
}

//////////////////////////////////////////////////////////////////////////////
//
// compute gradients on telescoping sets of eigenvectors and run
// non-maxima supression on the results.
//
void computeEigengradients(const Util::ImageStack& evec, Util::ImageStack& pbs)
{
  const int numVectors = evec.size(0);
  const int width = evec.size(1);
  const int height = evec.size(2);

  Util::Image g1(width,height); g1.init(0);
  Util::Image g2(width,height); g2.init(0);
  Util::Image g12(width,height); g12.init(0);
  Util::Image mag(width,height); 
  Util::Image ori(width,height); 
  Util::Image dx(width,height); 
  Util::Image dy(width,height);
  Util::Image gradient(width,height);

  pbs.resize(numVectors,width,height);

  //create filters
  double smooth = 0.5;
  Util::Image hFilt;
  Util::Image vFilt;
  Group::FilterBank::createGaussianKernel(smooth, smooth, 10.0, 1.5708, 1, false,hFilt);
  Group::FilterBank::createGaussianKernel(smooth, smooth, 10.0, 0, 1, false,vFilt);

  //filter each set of vectors up to i
  for (int i = 0; i < numVectors; i++)
  {
    getFiltered(*evec.slice(i),hFilt,dx);
    getFiltered(*evec.slice(i),vFilt,dy);
    for (int x = 0; x < width; x++)
    {
      for (int y = 0; y < height; y++)
      {
        double f1 = g1(x,y) + dx(x,y)*dx(x,y);
        double f2 = g2(x,y) + dy(x,y)*dy(x,y);
        double f12 = g12(x,y) + dx(x,y)*dy(x,y);
        double orival = 0.5*atan(2*f12 / (f1 - f2));
        if ((f1-f2) < 0)
        {
          if (2*f12>0)
          {
            orival+=M_PI/2;
          }
          else
          {
            orival-=M_PI/2;
          }
        }
        mag(x,y) =  0.5*(f1 + f2 + sqrt((f1-f2)*(f1-f2) + 4*f12*f12));
        ori(x,y) =  fmod(-orival+(M_PI/2),M_PI);
        g1(x,y) = f1;
        g2(x,y) = f2;
        g12(x,y) = f12;
      }
    }
    Group::nonmaxSuppress(mag,ori,gradient);
    for (int x = 0; x < width; x++)
    {
      for (int y = 0; y < height; y++)
      {
        pbs(i,x,y) = gradient(x,y);
      }
    }
  }
}


//////////////////////////////////////////////////////////////////////////////

void 
densify(const Util::Array2D<float>& W, const Util::Segmentation& seg, const int desiredNumSP, 
        Util::Segmentation& coarseSeg, Util::Array2D<float>& Wcoarse)
{
    static int level = 0;
    const int width = seg.getWidth();
    const int height = seg.getHeight();
    const int numP = W.size(0);
    level++;

    Util::Message::debug("computing eigenvectors of G");
    DMatrix* DW = new DMatrix(W);
    double* rowSum = DW->getRowSum();
    DW->normalizedLaplacian(rowSum);
    double* devec;
    double* deval;
    DW->eig(numVectors+1,&devec,&deval);
    DW->undoNormalizedLaplacian(rowSum);

    Util::Message::startBlock("finding superpixels");
    Util::Message::debug("transposing eigenvectors");
    Util::Array2D<float> embed(numP,numVectors);
    for (int vnum = 0; vnum < numVectors; vnum++)
    {
      for (int pnum = 0; pnum < numP; pnum++)
      {
        embed(pnum,vnum) = devec[(vnum+1)*numP + pnum] / devec[0*numP + pnum];
      }
    }

    Util::Message::startBlock(50,Util::String("running kmeans k=%d",desiredNumSP));
    Util::Image means(desiredNumSP, numVectors);
    Util::Array1D<int> bestMembership(numP);
    Util::kmeans(embed, desiredNumSP, 1e-4, 30, true, means, bestMembership);
    float bestNcut = DW->computeNCut(rowSum,bestMembership);
    for (int iter = 0; iter < 50; iter++)
    {
      Util::Message::stepBlock();
      Util::Array1D<int> membership(width*height);
      Util::kmeans(embed, desiredNumSP, 1e-4, 30, true, means, membership);
      float ncut = DW->computeNCut(rowSum,membership);
      if (ncut < bestNcut)
      {
        ncut = bestNcut;
        bestMembership = membership;
      }
    }
    Util::Message::endBlock();

    for (int x = 0; x < width; x++)
    {
      for (int y = 0; y < height; y++)
      {
        coarseSeg(x,y) = bestMembership(seg(x,y));
      }
    }

    Util::Message::debug("post-processing segmentation");
    coarseSeg.fragment();
    coarseSeg.mergeDust(dustThresh);
    Util::Array1D<int> segmentSizes;
    coarseSeg.segSizes(segmentSizes);
    int numSP = segmentSizes.size(0);
    Util::Message::debug(Util::String("found %d superpixels",numSP));
    Util::Message::endBlock();

    Util::Message::debug("coarsening graph");
    Util::Array1D<int> p2sp(numP);
    for (int x = 0; x < width; x++)
    {
      for (int y = 0; y < height; y++)
      {
        p2sp(seg(x,y)) = coarseSeg(x,y);
      }
    }
    Util::Array2D<float> Wsum(numSP,numSP); 
    Wsum.init(0);
    Util::Array2D<float> D(numSP,numSP); 
    D.init(0);
    for (int i = 0; i < numP; i++)
    {
      for (int j = 0; j < numP; j++)
      {
        //collapse the edge from i to j into the superpixel edge
        Wsum(p2sp(i),p2sp(j)) += W(i,j);
        D(p2sp(i),p2sp(j))++; //TODO???? should we be adding a different degree here?
      } 
    }

    Util::Message::debug("densifying graph");
    Wcoarse.resize(numSP,numSP);
    for (int i = 0; i < numSP; i++)
    {
      for (int j = 0; j < numSP; j++)
      {
        double dense = 0;
        if (D(i,j) != 0)
        {
          dense = segmentSizes(i)*(Wsum(i,j) / D(i,j))*segmentSizes(j);
        }
        Wcoarse(i,j) = dense;
      } 
    }

    if (evecFileName != NULL)
    {
      Util::ImageStack evecImage(numVectors,width,height);
      for (int vnum = 0; vnum < numVectors; vnum++)
      {
        for (int x = 0; x < width; x++)
        {
          for (int y = 0; y < height; y++)
          {
            int pnum = seg(x,y);
            evecImage(vnum,x,y) = devec[vnum*numP + pnum];
          }
        }
      }
      char filename[1000];
      sprintf(filename,"%s.%d",evecFileName,level);
      Util::Message::debug("writing eigenvector file");
      std::ofstream out(filename);
      if (!out.good()) {
        throw Util::Exception (Util::String ("Error opening '%s' for writing.", evecFileName));
      }
      out << evecImage;
      out.close();
    }

    if (graphFileName != NULL)
    {
      char filename[1000];
      sprintf(filename,"%s.%d",graphFileName,level);
      Util::Message::debug("writing coarsened graph file");
      std::ofstream out(filename);
      out << seg << W << Wsum << D << Wcoarse << coarseSeg;
      out.close();
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
    const int idiag = (int)sqrt((double)(width*width+height*height));
    const int numPixels = width*height;
    Util::Message::debug(Util::String("image size [%d x %d]",width,height),1);

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
    int smapRadius = 0;
    if (cueMode == Group::ic)
    {
      smapRadius = dthresh;    
    }
    else
    {
      smapRadius = Util::max(patchFactory.getSupportMapRadius(idiag),dthresh);
    }
    Util::Message::debug(Util::String("running intervening contour with map radius = %d", smapRadius));
    Group::SupportMap ic;
    Group::computeSupport(boundaries,smapRadius,1.0f,ic);
    Util::Message::endBlock();

    Util::Message::startBlock("acquiring region features");
    if ((cueMode == Group::patch) || (cueMode == Group::both))
    {
      Util::Message::debug("region processing");
      patchFactory.initialize(image,ic);
    }
    Util::Message::endBlock();

    //////////////////////////////////////////////////////////////////////////////

    Util::Message::debug("building affinity matrix");
    SMatrix* affinities = NULL;
    Group::computeAffinities(ic,patchFactory,cueMode,sigma,dthresh,useColor,&affinities);
    assert(affinities != NULL);

    //////////////////////////////////////////////////////////////////////////////

    if (smatrixFileName != NULL)
    {
      Util::Message::debug("writing affinity matrix");
      FILE* fp = fopen(smatrixFileName,"w");
      affinities->dump(fp);
      fclose(fp);
    }

    if ((evecFileName != NULL) || (egradFileName != NULL) || (graphFileName != NULL) || (segFileName != NULL))
    {
      Util::Message::debug("computing normalized laplacian");
      double* rowSum = affinities->getRowSum();
      affinities->normalizedLaplacian(rowSum);
      
      Util::Message::debug("running eignesolver");
      double* eval;
      double* evec;
      affinities->eigs(numVectors+1,&eval,&evec);

      //output raw eigenvectors if desired
      if (evecFileName != NULL)
      {
        Util::ImageStack evecImage(numVectors,width,height);
        for (int vnum = 0; vnum < numVectors; vnum++)
        {
          for (int x = 0; x < width; x++)
          {
            for (int y = 0; y < height; y++)
            {
              int pnum = y*width+x;
              evecImage(vnum,x,y) = evec[vnum*numPixels + pnum];
            }
          }
        }

        Util::Message::debug("writing eigenvector file");
        std::ofstream out(evecFileName);
        if (!out.good()) {
          throw Util::Exception (Util::String ("Error opening '%s' for writing.", evecFileName));
        }
        out << evecImage;
        out.close();
      }

      //divide thru by leading eigenvector, then median filter as desired
      Util::ImageStack evecImage(numVectors,width,height);
      Util::Message::debug("scaling and filtering eigenvectors");
      for (int vnum = 1; vnum < (numVectors+1); vnum++)
      {
        Util::Image im(width,height);
        for (int x = 0; x < width; x++)
        {
          for (int y = 0; y < height; y++)
          {
            int pnum = y*width+x;
            im(x,y) = evec[vnum*numPixels + pnum] / evec[0*numPixels + pnum];
          }
        }
        float norm = 0.0f;
        if (filterEvec)
        { 
          Util::Image mim(width,height);
          getPctFiltered(im,3,0.5,mim);
          for (int x = 0;  x < width; x++)
          {
            for (int y = 0; y < height; y++)
            {
              evecImage(vnum-1,x,y) = mim(x,y);    
              norm = norm+im(x,y)*im(x,y);
            }
          }
        }
        else
        {
          for (int x = 0;  x < width; x++)
          {
            for (int y = 0; y < height; y++)
            {
              evecImage(vnum-1,x,y) = im(x,y);    
              norm = norm+im(x,y)*im(x,y);
            }
          }
        }
        norm = sqrt(norm);
        for (int x = 0;  x < width; x++)
        {
          for (int y = 0; y < height; y++)
          {
            evecImage(vnum-1,x,y) = evecImage(vnum-1,x,y)/norm;
          }
        }
      }

      if (egradFileName != NULL)
      {
        Util::ImageStack pbs;
        computeEigengradients(evecImage,pbs);
        Util::Message::debug("writing eigenvector-gradient file");
        std::ofstream out(egradFileName);
        if (!out.good()) {
          throw Util::Exception (Util::String ("Error opening '%s' for writing.", egradFileName));
        }
        out << pbs;
        out.close();
      }

      if ((graphFileName != NULL) || (segFileName != NULL))
      {
        Util::Message::debug("transposing eigenvectors");
        Util::Array2D<float> embed(width*height,numVectors);
        for (int x = 0; x < width; x++)
        {
          for (int y = 0; y < height; y++)
          {
            for (int vnum = 0; vnum < numVectors; vnum++)
            {
              embed(x*height+y,vnum) = evecImage(vnum,x,y);
            }
          }
        }

        affinities->undoNormalizedLaplacian(rowSum);

        Util::Message::startBlock(5,"finding superpixels");
        Util::Image means(numSuperpixels, numVectors);
        Util::Array1D<int> bestMembership(width*height);
        Util::kmeans(embed, numSuperpixels, 1e-4, 30, true, means, bestMembership);
        float bestNcut = affinities->computeNCut(rowSum,bestMembership,numSuperpixels);
        for (int iter = 0; iter < 10; iter++)
        {
          Util::Message::stepBlock();
          Util::Array1D<int> membership(width*height);
          Util::kmeans(embed, numSuperpixels, 1e-4, 30, true, means, membership);
          float ncut = affinities->computeNCut(rowSum,membership,numSuperpixels);
          if (ncut < bestNcut)
          {
            ncut = bestNcut;
            bestMembership = membership;
          }
        }
        Util::Message::endBlock();


        Util::Segmentation seg(width,height); 
        for (int x = 0; x < width; x++)
        {
          for (int y = 0; y < height; y++)
          {
            seg(x,y) = bestMembership(x*height+y);
          }
        }

        Util::Message::debug("post-processing segments");
        seg.fragment();
        seg.mergeDust(dustThresh);
        Util::Array1D<int> segmentSizes;
        seg.segSizes(segmentSizes);
        const int numSP = segmentSizes.size(0);
        Util::Message::debug(Util::String("found %d superpixels",numSP));

        if (segFileName != NULL)
        {
          Util::Message::debug(Util::String("writing segmentation file %s",segFileName));
          std::ofstream out(segFileName);
          out << seg;
          out.close();
        }

        Util::Message::debug("coarsening graph G_1 -> G_2");
        Util::Array2D<float> Wsum(numSP,numSP);
        Wsum.init(0);
        Util::Array2D<float> D(numSP,numSP);
        D.init(0);
        for (int r = 0; r < affinities->n; r++)
        {
          int rx = r % width;
          int ry = (r - rx) / width;
          int rseg = seg(rx,ry);
          for (int i = 0; i < affinities->nz[r]; i++)
          {
            int c = affinities->col[r][i];       
            int cx = c % width;
            int cy = (c - cx) / width;
            int cseg = seg(cx,cy);
            Wsum(rseg,cseg) += affinities->values[r][i];
            D(rseg,cseg)++;
          }
        }

        Util::Message::debug("densifying G_2");
        Util::Array2D<float> W2(numSP,numSP);
        for (int i = 0; i < numSP; i++)
        {
          for (int j = 0; j < numSP; j++)
          {
            float dense = 0;
            if (D(i,j) != 0)
            {
              dense = segmentSizes(i)*(Wsum(i,j) / D(i,j))*segmentSizes(j);
            }
            W2(i,j) = dense;
          } 
        }

        if (graphFileName != NULL)
        {
          Util::Message::debug(Util::String("writing coarsened graph file %s",graphFileName));
          std::ofstream out(graphFileName);
          out << seg << Wsum << D << W2;
          out.close();
        }
     
     /*
        Util::Array2D<float> W3;
        Util::Segmentation seg3(width,height);
        densify(W2, seg, 80, seg3, W3);

        Util::Array2D<float> W4;
        Util::Segmentation seg4(width,height);
        densify(W3, seg3, 60, seg4, W4);

        Util::Array2D<float> W5;
        Util::Segmentation seg5(width,height);
        densify(W4, seg4, 30, seg5, W5);

        Util::Array2D<float> W6;
        Util::Segmentation seg6(width,height);
        densify(W5, seg5, 20, seg6, W6);

        Util::Array2D<float> W7;
        Util::Segmentation seg7(width,height);
        densify(W6, seg6, 10, seg7, W7);
      */  

        std::cerr << "done." << std::endl;
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
        Util::Message::debug("EXITING WITH ERROR!!!!");
        status = 0; //exit with status 0 even though there was a problem so that milrun still goes
    }
    return status;
}


