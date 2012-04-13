#include "DomainInfo.hh"
#include <cmath>
#include <algorithm>
#include "IndexToTuple.hh"
#include "IndexToVector.hh"
#include "three_algebra.h"

using namespace std;

DomainInfo::DomainInfo(const vector<Long64>& gid, int nx, int ny, int nz)
: radius_(0.0), center_(0., 0., 0.)
{
   IndexToVector indexToVector(nx, ny, nz);

   for (unsigned ii=0; ii<gid.size(); ++ii)
      center_ += indexToVector(gid[ii]);

   if (gid.size() > 0)
      center_ /= double(gid.size());
   
   for (unsigned ii=0; ii<gid.size(); ++ii)
   {
      Vector vv = indexToVector(gid[ii]);
      double rSq = dot( (vv-center_), (vv-center_) );
      radius_ = max(radius_, rSq);
   }
   radius_ = sqrt(radius_);
}

