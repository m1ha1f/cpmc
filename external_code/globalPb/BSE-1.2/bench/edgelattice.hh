
#ifndef __EdgeLattice_hh__
#define __EdgeLattice_hh__

#include <math.h>
#include <assert.h>

// Pixel edge lattice template class.
// Pixel centers are at (0,0)..(w-1,h-1).
// Horizontal edges are at (0,.5)..(w-1,h-1-.5).
// Vertical edges are at (.5,0)..(w-1-.5,h-1).

template<class T>
class EdgeLattice
{
public:
    EdgeLattice () 
        : h(0,0), v(0,0), width(0), height(0)
        {}
    EdgeLattice (int _width, int _height) 
        : h(_width,_height-1),
          v(_width-1,_height),
          width(_width),
          height(_height)
        {}
    int getWidth () const { return width; }
    int getHeight () const { return height; }
    bool issize (const int _width, const int _height) const {
        return (width == _width) && (height == _height);
    }
    void init (T x) {
        h.init(x); v.init(x);
    }
    void resize (const int _width, const int _height) {
        h.resize(_width,_height-1);
        v.resize(_width-1,_height);
        width = _width;
        height = _height;
    }
    // Is (x,y) a valid pixel edge coordinate?
    bool valid (const double x, const double y) const { 
        return x >= 0 && x <= width-1
            && y >= 0 && y <= height-1
            && (isH(x,y) || isV(x,y));
    }
    // Does (x,y) denote a horizontal pixel edge?
    static bool isH (const double x, const double y) {
        return (x - floor(x) == 0) && (y - floor(y) == 0.5);
    }
    // Does (x,y) denote a vertical pixel edge?
    static bool isV (const double x, const double y) {
        return (x - floor(x) == 0.5) && (y - floor(y) == 0);
    }
    // Read-only accessor.
    const T& operator() (const double x, const double y) const {
        assert(valid(x,y));
        const int i = (int) floor (x);
        const int j = (int) floor (y);
        return isH(x,y) ? h(i,j) : v(i,j);
    }
    // Read-write accessor.
    T& operator() (const double x, const double y) {
        assert(valid(x,y));
        const int i = (int) floor (x);
        const int j = (int) floor (y);
        return isH(x,y) ? h(i,j) : v(i,j);
    }
private:
    Util::Array2D<T> h,v;
    int width, height;
};

#endif // __EdgeLattice_hh__
