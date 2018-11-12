#ifndef BOUNDING_BOX_HH
#define BOUNDING_BOX_HH

#include <vector>
#include "Tuple.hh"

class BoundingBox
{
 public:
   BoundingBox(const Tuple& minCorner, const Tuple& maxCorner);
   BoundingBox(const std::vector<Tuple>& tt);
   
   bool overlap(const BoundingBox& b2) const;
   bool contains(const Tuple& tt) const;

 private:
   Tuple minCorner_;
   Tuple maxCorner_;

};

#endif
