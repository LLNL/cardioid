#include "BoundingBox.hh"
#include <cmath>

using namespace std;

BoundingBox::BoundingBox(const Tuple& minCorner, const Tuple& maxCorner)
: minCorner_(minCorner),
  maxCorner_(maxCorner){}

BoundingBox::BoundingBox(const vector<Tuple>& tt)
: minCorner_(Tuple(0, 0, 0)),
  maxCorner_(Tuple(0, 0, 0))
{
   if (tt.size() == 0)
      return;
   minCorner_ = tt[0];
   maxCorner_ = tt[0];
   for (unsigned ii=1; ii<tt.size(); ++ii)
   {
      minCorner_.x() = min(tt[ii].x(), minCorner_.x());
      minCorner_.y() = min(tt[ii].y(), minCorner_.y());
      minCorner_.z() = min(tt[ii].z(), minCorner_.z());
      maxCorner_.x() = max(tt[ii].x(), maxCorner_.x());
      maxCorner_.y() = max(tt[ii].y(), maxCorner_.y());
      maxCorner_.z() = max(tt[ii].z(), maxCorner_.z());
   }
}
   
bool BoundingBox::overlap(const BoundingBox& b2) const 
{
   if (b2.minCorner_.x() > this->maxCorner_.x()) return false;
   if (b2.maxCorner_.x() < this->minCorner_.x()) return false;
   if (b2.minCorner_.y() > this->maxCorner_.y()) return false;
   if (b2.maxCorner_.y() < this->minCorner_.y()) return false;
   if (b2.minCorner_.z() > this->maxCorner_.z()) return false;
   if (b2.maxCorner_.z() < this->minCorner_.z()) return false;
   return true;
}

bool BoundingBox::contains(const Tuple& tt) const 
{
   if ( tt.x() < minCorner_.x() ) return false;
   if ( tt.x() > maxCorner_.x() ) return false;
   if ( tt.y() < minCorner_.y() ) return false;
   if ( tt.y() > maxCorner_.y() ) return false;
   if ( tt.z() < minCorner_.z() ) return false;
   if ( tt.z() > maxCorner_.z() ) return false;

   return true;
}
