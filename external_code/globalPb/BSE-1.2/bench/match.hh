
#ifndef __match_hh__
#define __match_hh__

#include "point.hh"
#include "edgelattice.hh"

// Compute min-cost perfect match using Andrew Goldberg's CSA.
// There must be a perfect match possible, or CSA may not terminate.
void minCostMaxMatch (
    // number of nodes on each side
    const int n, 
    // number of edges
    const int m, 
    // m x 3 array of edges (i,j,w), where i,j in [0,n) and w is the
    // edge weight.
    const Util::Array2D<int>& edges, 
    // n x 3 array of edges in assignment, same structure as edges.
    Util::Array2D<double> match); 

// Correspond two boundary maps, and return the assignment cost.
double matchEdgeMaps (
    // Dimensions of the image.
    const int width, const int height,
    // Boundary maps, width x height.
    const Util::Array2D<bool>& map1,
    const Util::Array2D<bool>& map2,
    // Maximum distance at which pixels can be matched.
    const double maxDist,
    // The cost of not matching a boundary pixel.
    const double outlierCost,
    // Match edges, width x height.
    // Set to a valid pixel for matched items.
    // Set to (-1,-1) for unmatched items, for both non-boundary pixels
    // and unmatched boundary pixels.
    Util::Array2D<Pixel>& match1,
    Util::Array2D<Pixel>& match2);

// Correspond two segmentations, given as edge lattices along with
// edgel orientations.  Return the assignment cost.  Arguments similar
// to matchEdgeMaps.
double matchSegs (
    const EdgeLattice<bool>& map1,
    const EdgeLattice<bool>& map2,
    const EdgeLattice<double>& theta1,
    const EdgeLattice<double>& theta2,
    const double maxDist, 
    const double maxTheta,
    const double outlierCost,
    EdgeLattice< Point2D<double> >& match1,
    EdgeLattice< Point2D<double> >& match2);

#endif // __Match_hh__
