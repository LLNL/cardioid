#ifndef DOMAIN_INFO_HH
#define DOMAIN_INFO_HH

#include <vector>
#include "Long64.hh"
#include "Vector.hh"
#include "Tuple.hh"
#include "BoundingBox.hh"

class DomainInfo
{
 public:
   DomainInfo(const std::vector<Long64>& gid, int nx, int ny, int nz);
   
   const Vector& center() const {return center_;}
   double radius() const {return radius_;}
   int nCells() const {return nCells_;}
   BoundingBox boundingBox() const {return BoundingBox(minCorner_, maxCorner_);}

 private:

   int    nCells_;
   double radius_;
   Vector center_;
   Tuple  minCorner_;
   Tuple  maxCorner_;
   
};

#endif
