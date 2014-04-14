#include <math.h>
#include "common.h"
#include "matrix.h"
#include "clique.h"
#include "Util.h"

//
// find cliques in the graph given by thresholding the intervening
// contour similarity.  only consider quadtree cells so we can exploit
// the geometric fact that a convex region containing no above threshold
// edges will be a clique in the final graph.
// (these are of course not maximal cliques) 
//
void 
findCliques(const Matrix* pb, const double threshold, 
             Matrix* cliqueMembership, int* numCliques)
{
  int width = cliqueMembership->getWidth();
  int height = cliqueMembership->getHeight();
  Matrix* assigned = new Matrix(width,height);
  Matrix* oldassigned = new Matrix(width,height);

  //mark each pixel as a singleton clique
  //if a pixel is above threshold, mark it as assigned
  int nextclique = 0;
  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x++)
    {
      cliqueMembership->set(y,x,nextclique);
      nextclique++;
      if (pb->get(y,x) < threshold)
      {
        oldassigned->set(y,x,0);
        assigned->set(y,x,0);
      }
      else
      {
        oldassigned->set(y,x,1);
        assigned->set(y,x,1);
      }
    }
  }  

  //now merge upwards searching for cliques
  //in quadtree cells.
  bool done = false;
  int csize = 2;
  while (csize < Util::min(width,height))
  {
    for (int y = 0; y < height-csize; y+=csize)
    {
      for (int x = 0; x < width-csize; x+=csize)
      {
        //compute the number of neighbors who have been assigned at
        //the previous level
        int ass = 0;
        for (int dy = Util::max(y-(csize/2),0); dy < Util::min(y+csize+(csize/2),height-1); dy++)
        {
          for (int dx = Util::max(x-(csize/2),0); dx < Util::min(x+csize+(csize/2),width-1); dx++)
          {
            ass = ass + oldassigned->get(dy,dx);
          }
        }

        //in order to assure smoothness
        // 1. if no neighbors have been assigned, merge
        // 2. otherwise, assign at this level
        if (ass > 0)
        {
          for (int dy = y; dy < y+csize; dy++)
          {
            for (int dx = x; dx < x+csize; dx++)
            {
              assigned->set(dy,dx,1);
            } 
          }
        }
        else
        {
          for (int dy = y; dy < y+csize; dy++)
          {
            for (int dx = x; dx < x+csize; dx++)
            {
              cliqueMembership->set(dy,dx,nextclique);
            } 
          }
          nextclique++;
        }
      }
    }
    csize = csize*2;
    for (int y = 0; y < height; y++)
    {
      for (int x = 0; x < width; x++)
      {
        oldassigned->set(y,x,assigned->get(y,x));
      }
    }  
  }

  //now relabel cliques to get a seqential list
  int* unique = new int[width*height];  //there can't be more labels than pixels
  memset(unique,-1,width*height*sizeof(int));
  int endlist = 0;
  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x++)
    {
      int label = cliqueMembership->get(y,x);
      bool found = false;
      int i = 0;
      while (!found && (i < endlist))
      {
        found = (unique[i] == label);
        i++;
      }
      if (found)
      {
        cliqueMembership->set(y,x,i-1);
      }
      else
      {
        unique[endlist] = label;
        cliqueMembership->set(y,x,endlist);
        endlist++;
      }
    }
  }
  *numCliques = endlist;
}


//
// compute the number of pixels in each clique
//
void
getCliqueSize(const Matrix* cliqueMembership, const int numCliques, int* cliqueSize)
{
  //zero out entries
  //memset(cliqueSize,0,numCliques*sizeof(int));

  //count the number of pixels in each clique 
  for (int y = 0; y < cliqueMembership->getHeight(); y++)
  {
    for (int x = 0; x < cliqueMembership->getWidth(); x++)
    {
      int cnum = (int)cliqueMembership->get(y,x);
      assert(cnum < numCliques);
      cliqueSize[cnum]++;
    }
  }
}


//add value to list if it's unique
void
addunique(int* nbrlist, const int len, const int newval, int exclude)
{
  if (newval != exclude)
  {
    bool found = false;
    bool end = false;
    int i = 0;
    while (!found && !end)
    {
      if (nbrlist[i] == newval)
      {
        found = true;
      }
      if (nbrlist[i] == -1)
      {
        end = true;
      }
      i++;
    }
    if (end)
    {
      nbrlist[i-1] = newval; 
    }
  }
}

//
// find neigboring cliques of each clique.  nbrhood does not include the clique itself
//
void
findNbrs(const Matrix* cliqueMembership, const int numCliques, int** nbrs, int* nnbrs)
{
  int width = cliqueMembership->getWidth();
  int height = cliqueMembership->getHeight();
  *nbrs = new int[numCliques*12];  //at most 12 neighbors
  *nnbrs = 12;
  memset(*nbrs,-1,numCliques*12*sizeof(int));
  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x++)
    {
      int thisclique = cliqueMembership->get(y,x); 
      addunique(&((*nbrs)[12*thisclique]),12,cliqueMembership->get(y,Util::min(x+1,width-1)),thisclique); 
      addunique(&((*nbrs)[12*thisclique]),12,cliqueMembership->get(Util::min(y+1,height-1),x),thisclique); 
      addunique(&((*nbrs)[12*thisclique]),12,cliqueMembership->get(y,Util::max(x-1,0)),thisclique); 
      addunique(&((*nbrs)[12*thisclique]),12,cliqueMembership->get(Util::max(y-1,0),x),thisclique); 
      addunique(&((*nbrs)[12*thisclique]),12,cliqueMembership->get(Util::min(y+1,height-1),Util::min(x+1,width-1)),thisclique); 
      addunique(&((*nbrs)[12*thisclique]),12,cliqueMembership->get(Util::min(y+1,height-1),Util::max(x-1,0)),thisclique); 
      addunique(&((*nbrs)[12*thisclique]),12,cliqueMembership->get(Util::max(y-1,0),Util::max(x-1,0)),thisclique); 
      addunique(&((*nbrs)[12*thisclique]),12,cliqueMembership->get(Util::max(y-1,0),Util::min(x+1,width-1)),thisclique); 
    }
  }
}


void
lift(const double* cliques, const Matrix* cliqueMembership, Matrix* pixels)
{
  int width = cliqueMembership->getWidth();
  int height = cliqueMembership->getHeight();
  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x++)
    {
      int cnum = (int)(cliqueMembership->get(y,x));
      pixels->set(y,x,cliques[cnum]);
    }
  }
}

void
project(const Matrix* pixels, const Matrix* cliqueMembership, const int numCliques, double* cliques)
{
  memset(cliques,0,numCliques*sizeof(double));
  int* count = new int[numCliques];
  memset(count,0,numCliques*sizeof(int));
  int width = cliqueMembership->getWidth();
  int height = cliqueMembership->getHeight();
  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x++)
    {
      int cnum = (int)(cliqueMembership->get(y,x));
      cliques[cnum] += pixels->get(y,x);
      count[cnum]++;
    }
  }

  for (int i = 0; i < numCliques; i++)
  {
    if (count[i] > 0)
    {
      cliques[i] = cliques[i] / count[i];  
    }
    else
    {
      cliques[i] = 0;
    }
  }
  delete[] count;
}


void
lift(const int* cliques, const Matrix* cliqueMembership, Matrix* pixels)
{
  int width = cliqueMembership->getWidth();
  int height = cliqueMembership->getHeight();
  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x++)
    {
      int cnum = (int)(cliqueMembership->get(y,x));
      pixels->set(y,x,cliques[cnum]);
    }
  }
}


void
project(const Matrix* pixels, const Matrix* cliqueMembership, const int numCliques, int* cliques)
{
  memset(cliques,0,numCliques*sizeof(int));
  int* count = new int[numCliques];
  memset(count,0,numCliques*sizeof(int));
  int width = cliqueMembership->getWidth();
  int height = cliqueMembership->getHeight();

  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x++)
    {
      int cnum = (int)(cliqueMembership->get(y,x));
      cliques[cnum] += pixels->get(y,x);
      count[cnum]++;
    }
  }

  for (int i = 0; i < numCliques; i++)
  {
    if (count[i] > 0)
    {
      cliques[i] = (int)floor(cliques[i] / count[i]);  
    }
    else
    {
      cliques[i] = 0;
    }
  }
  delete[] count;
}

