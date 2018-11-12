#ifndef LOCAL_GRID_HH
#define LOCAL_GRID_HH

#include "Tuple.hh"
#include <cassert>

class LocalGrid
{
 public:

   LocalGrid(int nx, int ny, int nz, int xOff, int yOff, int zOff);

   int nx() const; 
   int ny() const; 
   int nz() const; 
   
   Tuple localTuple(const Tuple& globalTuple) const;
   Tuple globalTuple(const Tuple& localTuple) const;

 private:
   int nx_; // The size of the grid in each direction
   int ny_;
   int nz_;

   int xOffset_; // These are the coordinates of the minimum
   int yOffset_; // corner of the local grid on the global 
   int zOffset_; // grid.  Yes, they can be negative.
   
};

inline
LocalGrid::LocalGrid(int nx, int ny, int nz, int xOff, int yOff, int zOff)
: nx_(nx), ny_(ny), nz_(nz), xOffset_(xOff), yOffset_(yOff), zOffset_(zOff)
{}

inline int LocalGrid::nx() const {return nx_;}
inline int LocalGrid::ny() const {return ny_;}
inline int LocalGrid::nz() const {return nz_;}



inline Tuple
LocalGrid::localTuple(const Tuple& globalTuple) const
{
   return globalTuple - Tuple(xOffset_, yOffset_, zOffset_);
}

inline Tuple
LocalGrid::globalTuple(const Tuple& localTuple) const
{
   return localTuple + Tuple(xOffset_, yOffset_, zOffset_);
}

#endif
