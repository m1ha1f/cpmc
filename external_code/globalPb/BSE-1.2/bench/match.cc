#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <fstream>
#include <sys/types.h>
#include <signal.h>
#include <sys/wait.h>
#include "util.hh"
#include "array.hh"
#include "match.hh"
#include "string.hh"
#include "exception.hh"

//////////////////////////////////////////////////////////////////////
// minCostMaxMatch
//////////////////////////////////////////////////////////////////////

// Uses Andrew Goldberg's CSA program, the prec_costs version with
// options QUICK_MIN and MIN_COST.

// CSA input ASN assignment problem format 
// p asn N M
// n i		i=[1,N]
// a i j w	i=[1,N] j=[N+1,2N]

// CSA output assignment solution format
// f i j -w	i=[1,N] j=[N+1,2N]

// Compute min-cost perfect match using CSA.
// There must be a perfect match possible, or CSA may not terminate.
void
minCostMaxMatch (
    const int n, const int m,
    const Util::Array2D<int>& edges,
    Util::Array2D<int>& match) {

    // Check the input graph.
    assert (edges.issize(m,3));
    for (int i = 0; i < m; i++) {
        const int a = edges(i,0);
        const int b = edges(i,1);
        const int w = edges(i,2);
        assert (a >= 0 && a < n);
        assert (b >= 0 && b < n);
        assert (w >= 0); // non-negative edge weights only
    }

    // Output match array should be n x 3.
    assert (match.issize(n,3));

    // No perfect assignment is possible if m < n.
    assert (m >= n);

    FILE* csaIn = NULL;
    FILE* csaOut = NULL;
    int pid = 0;
    try {
        // Start up CSA.
        char* const argv[] = { "csa", NULL };
        const char* program = "csa";
        pid = Util::systemPipe (program, argv, csaIn, csaOut);

        // Send the graph in ASN format to CSA.
        fprintf (csaIn, "p asn %d %d\n", 2*n, m); // problem line
        for (int i = 0; i < n; i++) { // node list
            fprintf (csaIn, "n %d\n", i+1);
        }
        for (int i = 0; i < m; i++) { // arc list
            fprintf (csaIn, "a %d %d %d\n",
                     1+edges(i,0), 1+n+edges(i,1), edges(i,2));
                 
        }
        fclose(csaIn);
        csaIn = NULL;

        // Read response from CSA.
        Util::String line;
        for (int i = 0; i < n; i++) {
            if (!line.nextLine(csaOut)) {
                throw Util::Exception ("minCostMaxMatch: unexpected EOF parsing CSA output.");
            }
            int a,b,c;
            int cnt = sscanf (line.text(), "f %d %d %d", &a, &b, &c);
            if (cnt != 3
                || a < 1 || a > n
                || b < n+1 || b > 2*n
                || c > 0) {
                throw Util::Exception (Util::String (
                    "minCostMaxMatch: invalid CSA output line '%s'.", line.text()));
            }
            assert (cnt == 3);
            assert (a >= 1 && a <= n);
            assert (b >= n+1 && b <= 2*n);
            assert (c <= 0);
            match(i,0) = a - 1;
            match(i,1) = b - n - 1;
            match(i,2) = -c;
        }
        if (line.nextLine(csaOut)) {
            throw Util::Exception (Util::String (
                "minCostMaxMatch: expected EOF, found '%s'.", line.text()));
        }
        fclose(csaOut);
        csaOut = NULL;

        // Done with CSA.  Make sure it's dead.
        kill(pid,SIGKILL);
        waitpid(pid,NULL,0);
        pid = 0;
    }
    catch (Util::Exception& e) {
        // Clean up and rethrow.
        if (csaIn != NULL) { fclose(csaIn); }
        if (csaOut != NULL) { fclose(csaOut); }
        if (pid != 0) { 
            kill(pid,SIGKILL); 
            waitpid(pid,NULL,0);
        }
        throw;
    }
}

//////////////////////////////////////////////////////////////////////
// matchEdgeMaps
//////////////////////////////////////////////////////////////////////

struct Edge {
    int i,j;	// node ids, 0-based
    double w;	// distance between pixels
};

// CSA code needs integer weights.  Use this multiplier to convert
// floating-point weights to integers.
static const int multiplier = 100;

// The degree of outlier connections.
static const int degree = 6;

double matchEdgeMaps (
    const int width, const int height,
    const Util::Array2D<bool>& map1,
    const Util::Array2D<bool>& map2,
    const double maxDist,
    const double outlierCost,
    Util::Array2D<Pixel>& match1,
    Util::Array2D<Pixel>& match2)
{
    // Check global constants.
    assert (degree > 0);
    assert (multiplier > 0);

    // Max distance must be non-negative.
    assert (maxDist >= 0);

    // Outlier cost must be larger than the max distance.
    assert (outlierCost > maxDist);

    // Check input array dimensions.
    assert (map1.issize(width,height));
    assert (map2.issize(width,height));

    // Check output array dimensions.
    assert (match1.issize(width,height));
    assert (match2.issize(width,height));

    // Initialize match[12] arrays to (-1,-1).
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            match1(x,y) = Pixel(-1,-1);
            match2(x,y) = Pixel(-1,-1);
        }
    }

    // Radius of search window.
    const int r = (int) ceil (maxDist);	

    // Figure out which nodes are matchable, i.e. within maxDist
    // of another node.
    Util::Array2D<bool> matchable1 (width,height);
    Util::Array2D<bool> matchable2 (width,height);
    matchable1.init(false);
    matchable2.init(false);
    for (int x1 = 0; x1 < width; x1++) {
        for (int y1 = 0; y1 < height; y1++) {
            if (!map1(x1,y1)) { continue; }
            for (int u = -r; u <= r; u++) {
                for (int v = -r; v <= r; v++) {
                    const double d2 = u*u + v*v;
                    if (d2 > maxDist*maxDist) { continue; }
                    const int x2 = x1 + u;
                    const int y2 = y1 + v;
                    if (x2 < 0 || x2 >= width) { continue; }
                    if (y2 < 0 || y2 >= height) { continue; }
                    if (!map2(x2,y2)) { continue; }
                    matchable1(x1,y1) = true;
                    matchable2(x2,y2) = true;
                }
            }
        }
    }

    // Count the number of nodes on each side of the match.
    // Construct nodeID->pixel and pixel->nodeID maps.
    // Node IDs range from [0,n1) and [0,n2).
    int n1=0, n2=0;
    std::vector<Pixel> nodeToPix1;
    std::vector<Pixel> nodeToPix2;
    Util::Array2D<int> pixToNode1 (width,height);
    Util::Array2D<int> pixToNode2 (width,height);
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            pixToNode1(x,y) = -1;
            pixToNode2(x,y) = -1;
            Pixel pix (x,y);
            if (matchable1(x,y)) {
                pixToNode1(x,y) = n1;
                nodeToPix1.push_back(pix);
                n1++;
            }
            if (matchable2(x,y)) {
                pixToNode2(x,y) = n2;
                nodeToPix2.push_back(pix);
                n2++;
            }
        }
    }

    // Construct the list of edges between pixels within maxDist.
    std::vector<Edge> edges;
    for (int x1 = 0; x1 < width; x1++) {
        for (int y1 = 0; y1 < height; y1++) {
            if (!matchable1(x1,y1)) { continue; }
            for (int u = -r; u <= r; u++) {
                for (int v = -r; v <= r; v++) {
                    const double d2 = u*u + v*v;
                    if (d2 > maxDist*maxDist) { continue; }
                    const int x2 = x1 + u;
                    const int y2 = y1 + v;
                    if (x2 < 0 || x2 >= width) { continue; }
                    if (y2 < 0 || y2 >= height) { continue; }
                    if (!matchable2(x2,y2)) { continue; }
                    Edge e; 
                    e.i = pixToNode1(x1,y1);
                    e.j = pixToNode2(x2,y2);
                    e.w = sqrt(d2);
                    assert (e.i != -1);
                    assert (e.j != -1);
                    assert (e.w < outlierCost);
                    edges.push_back(e);
                }
            }
        }
    }

    // The cardinality of the match is n.
    const int n = n1 + n2;
    const int nmin = Util::min(n1,n2);
    const int nmax = Util::max(n1,n2);

    // Compute the degree of various outlier connections.
    const int d1 = Util::minmax(0,degree,n1-1); // from map1
    const int d2 = Util::minmax(0,degree,n2-1); // from map2
    const int d3 = Util::min(degree,n1,n2); // between outliers
    const int dmax = Util::max(d1,d2,d3);

    assert (n1 == 0 || (d1 >= 0 && d1 < n1));
    assert (n2 == 0 || (d2 >= 0 && d2 < n2));
    assert (d3 >= 0 && d3 <= nmin);

    // Count the number of edges.
    int m = 0;
    m += edges.size(); 	// real connections
    m += d1 * n1;	// outlier connections
    m += d2 * n2;	// outlier connections
    m += d3 * nmax;	// outlier-outlier connections
    m += n; 		// high-cost perfect match overlay

    // Weight of outlier connections.
    const int ow = (int) ceil (outlierCost * multiplier);

    // Scratch array for outlier edges.
    Util::Array1D<uint> outliers (dmax);

    // Construct the input graph for the assignment problem.
    Util::Array2D<int> igraph (m,3);
    int count = 0;
    // real edges
    for (int a = 0; a < (int)edges.size(); a++) {
        igraph(count,0) = edges[a].i;
        igraph(count,1) = edges[a].j;
        igraph(count,2) = (int) rint (edges[a].w * multiplier);
        count++;
    }
    // outliers edges for map1, exclude diagonal
    for (int i = 0; i < n1; i++) {
        Util::kOfN(d1,n1-1,outliers.data());
        for (int a = 0; a < d1; a++) {
            int j = outliers(a);
            if (j >= i) { j++; }
            assert (i != j);
            assert (j >= 0 && j < n1);
            igraph(count,0) = i;
            igraph(count,1) = n2 + j;
            igraph(count,2) = ow;
            count++;
        }
    }
    // outliers edges for map2, exclude diagonal
    for (int j = 0; j < n2; j++) {
        Util::kOfN(d2,n2-1,outliers.data());
        for (int a = 0; a < d2; a++) {
            int i = outliers(a);
            if (i >= j) { i++; }
            assert (i != j);
            assert (i >= 0 && i < n2);
            igraph(count,0) = n1 + i;
            igraph(count,1) = j;
            igraph(count,2) = ow;
            count++;
        }
    }
    // outlier-to-outlier edges
    for (int i = 0; i < nmax; i++) {
        Util::kOfN(d3,nmin,outliers.data());
        for (int a = 0; a < d3; a++) {
            const int j = outliers(a);
            assert (j >= 0 && j < nmin);
            if (n1 < n2) {
                igraph(count,0) = n1 + i;
                igraph(count,1) = n2 + j;
            } else {
                igraph(count,0) = n1 + j;
                igraph(count,1) = n2 + i;
            }
            igraph(count,2) = ow;
            count++;
        }
    }
    // perfect match overlay (diagonal)
    for (int i = 0; i < n1; i++) {
        igraph(count,0) = i;
        igraph(count,1) = n2 + i;
        igraph(count,2) = ow * multiplier;
        count++;
    }
    for (int i = 0; i < n2; i++) {
        igraph(count,0) = n1 + i;
        igraph(count,1) = i;
        igraph(count,2) = ow * multiplier;
        count++;
    }
    assert (count == m);

    // Solve the assignment problem.
    Util::Array2D<int> ograph (n,3);
    (void) minCostMaxMatch (n, m, igraph, ograph);

    // Check the solution.
    // Count the number of high-cost edges from the perfect match
    // overlay that were used in the match.
    int overlayCount = 0;
    for (int a = 0; a < n; a++) {
        const int i = ograph(a,0);
        const int j = ograph(a,1);
        const int c = ograph(a,2);
        assert (i >= 0 && i < n);
        assert (j >= 0 && j < n);
        assert (c >= 0);
        // edge from high-cost perfect match overlay
        if (c == ow * multiplier) { overlayCount++; }
        // skip outlier edges
        if (i >= n1) { continue; }
        if (j >= n2) { continue; }
        // for edges between real nodes, check the edge weight
        const Pixel pix1 = nodeToPix1[i];
        const Pixel pix2 = nodeToPix2[j];
        const int dx = pix1.x - pix2.x;
        const int dy = pix1.y - pix2.y;
        const int w = (int) rint (sqrt(dx*dx+dy*dy)*multiplier);
        assert (w == c);
    }

    // Print a warning if any of the edges from the perfect match overlay
    // were used.  This should happen rarely.  If it happens frequently,
    // then the outlier connectivity should be increased.
    if (overlayCount > 0) {
        fprintf (stderr, "%s:%d: WARNING: The match includes %d "
                 "outliers from the perfect match overlay.\n",
                 __FILE__, __LINE__, overlayCount);
    }

    // Compute match arrays.
    for (int a = 0; a < n; a++) {
        // node ids
        const int i = ograph(a,0);
        const int j = ograph(a,1);
        // skip outlier edges
        if (i >= n1) { continue; }
        if (j >= n2) { continue; }
        // map node ids to pixels
        const Pixel pix1 = nodeToPix1[i];
        const Pixel pix2 = nodeToPix2[j];
        // record edges
        match1(pix1.x,pix1.y) = pix2;
        match2(pix2.x,pix2.y) = pix1;
    }

    // Compute the match cost.
    double cost = 0;
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            if (map1(x,y)) {
                if (match1(x,y) == Pixel(-1,-1)) {
                    cost += outlierCost;
                } else {
                    const int dx = x - match1(x,y).x;
                    const int dy = y - match1(x,y).y;
                    cost += 0.5 * sqrt (dx*dx + dy*dy);
                }
            }
            if (map2(x,y)) {
                if (match2(x,y) == Pixel(-1,-1)) {
                    cost += outlierCost;
                } else {
                    const int dx = x - match2(x,y).x;
                    const int dy = y - match2(x,y).y;
                    cost += 0.5 * sqrt (dx*dx + dy*dy);
                }
            }
        }
    }    

    // Return the match cost.
    return cost;
}

//////////////////////////////////////////////////////////////////////
// matchSegs
//////////////////////////////////////////////////////////////////////

double matchSegs (
    const EdgeLattice<bool>& map1,
    const EdgeLattice<bool>& map2,
    const EdgeLattice<double>& theta1,
    const EdgeLattice<double>& theta2,
    const double maxDist, 
    const double maxTheta,
    const double outlierCost,
    EdgeLattice< Point2D<double> >& match1,
    EdgeLattice< Point2D<double> >& match2)
{
    const int width = map1.getWidth();
    const int height = map1.getHeight();

    // Check global constants.
    assert (degree > 0);
    assert (multiplier > 0);

    // Max distance must be non-negative.
    assert (maxDist >= 0);

    // Max theta must be non-negative.
    assert (maxTheta >= 0);

    // Outlier cost must be larger than the max distance.
    assert (outlierCost > maxDist);

    // Check input array dimensions.
    assert (map1.issize(width,height));
    assert (map2.issize(width,height));
    assert (theta1.issize(width,height));
    assert (theta2.issize(width,height));

    // Check output array dimensions.
    assert (match1.issize(width,height));
    assert (match2.issize(width,height));

    // Initialize match[12] arrays to (-1,-1).
    for (double x = 0; x <= width; x += 0.5) {
        for (double y = 0; y <= height; y += 0.5) {
            if (!map1.valid(x,y)) { continue; }
            match1(x,y) = Point2D<double>(-1,-1);
            match2(x,y) = Point2D<double>(-1,-1);
        }
    }

    // Radius of search window.
    const int r = (int) ceil (maxDist);	

    // Figure out which edgels are matchable, i.e. within maxDist
    // and maxTheta of another edgel.
    EdgeLattice<bool> matchable1 (width,height);
    EdgeLattice<bool> matchable2 (width,height);
    // initialize arrays to all false
    for (double x = 0; x <= width; x += 0.5) {
        for (double y = 0; y <= height; y += 0.5) {
            // skip non-lattice points
            if (!map1.valid(x,y)) { continue; }
            // initialize lattice points to false
            matchable1(x,y) = false;
        }
    }
    // mark matchable edgels
    for (double x1 = 0; x1 <= width; x1 += 0.5) {
        for (double y1 = 0; y1 <= height; y1 += 0.5) {
            if (!map1.valid(x1,y1)) { continue; }
            if (!map1(x1,y1)) { continue; }
            for (double u = -r; u <= r; u += 0.5) {
                for (double v = -r; v <= r; v += 0.5) {
                    const double x2 = x1 + u;
                    const double y2 = y1 + v;
                    if (!map2.valid(x2,y2)) { continue; }
                    if (!map2(x2,y2)) { continue; }
                    const double d2 = u*u + v*v;
                    if (d2 > maxDist*maxDist) { continue; }
                    const double t1 = theta1(x1,y1);
                    const double t2 = theta2(x2,y2);
                    const double dt = Util::thetaDiff(t1,t2);
                    if (dt > maxTheta) { continue; }
                    matchable1(x1,y1) = true;
                    matchable2(x2,y2) = true;
                }
            }
        }
    }

    // Count the number of nodes on each side of the match.
    // Construct nodeID->edgel and edgel->nodeID maps.
    // Node IDs range from [0,n1) and [0,n2).
    int n1=0, n2=0;
    std::vector< Point2D<double> > nodeToEdgel1;
    std::vector< Point2D<double> > nodeToEdgel2;
    EdgeLattice<int> edgelToNode1 (width,height);
    EdgeLattice<int> edgelToNode2 (width,height);
    for (double x = 0; x <= width; x += 0.5) {
        for (double y = 0; y <= height; y += 0.5) {
            if (!map1.valid(x,y)) { continue; }
            edgelToNode1(x,y) = -1;
            edgelToNode2(x,y) = -1;
            Point2D<double> edgel (x,y);
            if (matchable1(x,y)) {
                edgelToNode1(x,y) = n1;
                nodeToEdgel1.push_back(edgel);
                n1++;
            }
            if (matchable2(x,y)) {
                edgelToNode2(x,y) = n2;
                nodeToEdgel2.push_back(edgel);
                n2++;
            }
        }
    }

    // Construct the list of edges between edgels within maxDist and
    // maxTheta.
    std::vector<Edge> edges;
    for (double x1 = 0; x1 <= width; x1 += 0.5) {
        for (double y1 = 0; y1 <= height; y1 += 0.5) {
            if (!map1.valid(x1,y1)) { continue; }
            if (!matchable1(x1,y1)) { continue; }
            for (double u = -r; u <= r; u += 0.5) {
                for (double v = -r; v <= r; v += 0.5) {
                    const double x2 = x1 + u;
                    const double y2 = y1 + v;
                    if (!map2.valid(x2,y2)) { continue; }
                    if (!matchable2(x2,y2)) { continue; }
                    const double d2 = u*u + v*v;
                    if (d2 > maxDist*maxDist) { continue; }
                    const double t1 = theta1(x1,y1);
                    const double t2 = theta2(x2,y2);
                    const double dt = Util::thetaDiff(t1,t2);
                    if (dt > maxTheta) { continue; }
                    Edge e; 
                    e.i = edgelToNode1(x1,y1);
                    e.j = edgelToNode2(x2,y2);
                    e.w = sqrt(d2) / maxDist + dt / maxTheta;
                    assert (e.i != -1);
                    assert (e.j != -1);
                    assert (e.w < outlierCost);
                    edges.push_back(e);
                }
            }
        }
    }

    // The cardinality of the match is n.
    const int n = n1 + n2;
    const int nmin = Util::min(n1,n2);
    const int nmax = Util::max(n1,n2);

    // Compute the degree of various outlier connections.
    const int d1 = Util::minmax(0,degree,n1-1); // from map1
    const int d2 = Util::minmax(0,degree,n2-1); // from map2
    const int d3 = Util::min(degree,n1,n2); // between outliers
    const int dmax = Util::max(d1,d2,d3);

    assert (n1 == 0 || (d1 >= 0 && d1 < n1));
    assert (n2 == 0 || (d2 >= 0 && d2 < n2));
    assert (d3 >= 0 && d3 <= nmin);

    // Count the number of edges.
    int m = 0;
    m += edges.size(); 	// real connections
    m += d1 * n1;	// outlier connections
    m += d2 * n2;	// outlier connections
    m += d3 * nmax;	// outlier-outlier connections
    m += n; 		// high-cost perfect match overlay

    // Weight of outlier connections.
    const int ow = (int) ceil (outlierCost * multiplier);

    // Scratch array for outlier edges.
    Util::Array1D<uint> outliers (dmax);

    // Construct the input graph for the assignment problem.
    Util::Array2D<int> igraph (m,3);
    int count = 0;
    // real edges
    for (int a = 0; a < (int)edges.size(); a++) {
        igraph(count,0) = edges[a].i;
        igraph(count,1) = edges[a].j;
        igraph(count,2) = (int) rint (edges[a].w * multiplier);
        count++;
    }
    // outliers edges for map1, exclude diagonal
    for (int i = 0; i < n1; i++) {
        Util::kOfN(d1,n1-1,outliers.data());
        for (int a = 0; a < d1; a++) {
            int j = outliers(a);
            if (j >= i) { j++; }
            assert (i != j);
            assert (j >= 0 && j < n1);
            igraph(count,0) = i;
            igraph(count,1) = n2 + j;
            igraph(count,2) = ow;
            count++;
        }
    }
    // outliers edges for map2, exclude diagonal
    for (int j = 0; j < n2; j++) {
        Util::kOfN(d2,n2-1,outliers.data());
        for (int a = 0; a < d2; a++) {
            int i = outliers(a);
            if (i >= j) { i++; }
            assert (i != j);
            assert (i >= 0 && i < n2);
            igraph(count,0) = n1 + i;
            igraph(count,1) = j;
            igraph(count,2) = ow;
            count++;
        }
    }
    // outlier-to-outlier edges
    for (int i = 0; i < nmax; i++) {
        Util::kOfN(d3,nmin,outliers.data());
        for (int a = 0; a < d3; a++) {
            const int j = outliers(a);
            assert (j >= 0 && j < nmin);
            if (n1 < n2) {
                igraph(count,0) = n1 + i;
                igraph(count,1) = n2 + j;
            } else {
                igraph(count,0) = n1 + j;
                igraph(count,1) = n2 + i;
            }
            igraph(count,2) = ow;
            count++;
        }
    }
    // perfect match overlay (diagonal)
    for (int i = 0; i < n1; i++) {
        igraph(count,0) = i;
        igraph(count,1) = n2 + i;
        igraph(count,2) = ow * multiplier;
        count++;
    }
    for (int i = 0; i < n2; i++) {
        igraph(count,0) = n1 + i;
        igraph(count,1) = i;
        igraph(count,2) = ow * multiplier;
        count++;
    }
    assert (count == m);

    // Solve the assignment problem.
    Util::Array2D<int> ograph (n,3);
    (void) minCostMaxMatch (n, m, igraph, ograph);

    // Check the solution.
    // Count the number of high-cost edges from the perfect match
    // overlay that were used in the match.
    int overlayCount = 0;
    for (int a = 0; a < n; a++) {
        const int i = ograph(a,0);
        const int j = ograph(a,1);
        const int c = ograph(a,2);
        assert (i >= 0 && i < n);
        assert (j >= 0 && j < n);
        assert (c >= 0);
        // edge from high-cost perfect match overlay
        if (c == ow * multiplier) { overlayCount++; }
        // skip outlier edges
        if (i >= n1) { continue; }
        if (j >= n2) { continue; }
        // for edges between real nodes, check the edge weight
        const Point2D<double> edgel1 = nodeToEdgel1[i];
        const Point2D<double> edgel2 = nodeToEdgel2[j];
        const double dx = edgel1.x - edgel2.x;
        const double dy = edgel1.y - edgel2.y;
        const double t1 = theta1(edgel1.x,edgel1.y);
        const double t2 = theta2(edgel2.x,edgel2.y);
        const double d2 = dx*dx + dy*dy;
        const double dt = Util::thetaDiff(t1,t2);
        const int w = (int) rint ((sqrt(d2)/maxDist + dt/maxTheta)*multiplier);
        assert (w == c);
    }

    // Print a warning if any of the edges from the perfect match overlay
    // were used.  This should happen rarely.  If it happens frequently,
    // then the outlier connectivity should be increased.
    if (overlayCount > 0) {
        fprintf (stderr, "%s:%d: WARNING: The match includes %d "
                 "outliers from the perfect match overlay.\n",
                 __FILE__, __LINE__, overlayCount);
    }

    // Compute match arrays.
    for (int a = 0; a < n; a++) {
        // node ids
        const int i = ograph(a,0);
        const int j = ograph(a,1);
        // skip outlier edges
        if (i >= n1) { continue; }
        if (j >= n2) { continue; }
        // map node ids to pixels
        const Point2D<double> edgel1 = nodeToEdgel1[i];
        const Point2D<double> edgel2 = nodeToEdgel2[j];
        // record edges
        match1(edgel1.x,edgel1.y) = edgel2;
        match2(edgel2.x,edgel2.y) = edgel1;
    }

    // Compute the match cost.
    double cost = 0;
    for (double x = 0; x <= width; x += 0.5) {
        for (double y = 0; y <= height; y += 0.5) {
            if (!map1.valid(x,y)) { continue; }
            if (map1(x,y)) {
                if (match1(x,y) == Point2D<double>(-1,-1)) {
                    cost += outlierCost;
                } else {
                    const double dx = x - match1(x,y).x;
                    const double dy = y - match1(x,y).y;
                    cost += 0.5 * sqrt (dx*dx + dy*dy);
                }
            }
            if (map2(x,y)) {
                if (match2(x,y) == Point2D<double>(-1,-1)) {
                    cost += outlierCost;
                } else {
                    const double dx = x - match2(x,y).x;
                    const double dy = y - match2(x,y).y;
                    cost += 0.5 * sqrt (dx*dx + dy*dy);
                }
            }
        }
    }    

    // Return the match cost.
    return cost;
}
