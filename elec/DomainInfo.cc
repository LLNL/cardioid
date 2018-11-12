#include "DomainInfo.hh"
#include <cmath>
#include <algorithm>
#include "IndexToTuple.hh"
#include "IndexToVector.hh"
#include "three_algebra.h"

using namespace std;

DomainInfo::DomainInfo(const vector<Long64>& gid, int nx, int ny, int nz)
: nCells_(gid.size()),
  radius_(0.0),
  center_(0., 0., 0.),
  minCorner_(Tuple(0, 0, 0)),
  maxCorner_(Tuple(0, 0, 0))
{

   if (gid.size() == 0)
      return;
   
   IndexToVector indexToVector(nx, ny, nz);

   for (unsigned ii=0; ii<gid.size(); ++ii)
      center_ += indexToVector(gid[ii]);

   center_ /= double(gid.size());
   
   for (unsigned ii=0; ii<gid.size(); ++ii)
   {
      Vector vv = indexToVector(gid[ii]);
      double rSq = dot( (vv-center_), (vv-center_) );
      radius_ = max(radius_, rSq);
   }
   radius_ = sqrt(radius_);

   IndexToTuple indexToTuple(nx, ny, nz);

   minCorner_ = indexToTuple(gid[0]);
   maxCorner_ = indexToTuple(gid[0]);
   for (unsigned ii=1; ii<nCells_; ++ii)
   {
      Tuple tt = indexToTuple(gid[ii]);
      minCorner_.x() = min(tt.x(), minCorner_.x());
      minCorner_.y() = min(tt.y(), minCorner_.y());
      minCorner_.z() = min(tt.z(), minCorner_.z());
      maxCorner_.x() = max(tt.x(), maxCorner_.x());
      maxCorner_.y() = max(tt.y(), maxCorner_.y());
      maxCorner_.z() = max(tt.z(), maxCorner_.z());
   }
}


