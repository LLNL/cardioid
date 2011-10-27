#include "DomainInfo.hh"
#include <cmath>
#include "IndexToTuple.hh"
#include "IndexToVector.hh"
#include "three_algebra.h"

using namespace std;

DomainInfo::DomainInfo(const vector<Long64>& gid, int nx, int ny, int nz)
: radius_(0.0), center_(3, 0.0)
{
   IndexToTuple indexToTuple(nx, ny, nz);

   {//limit scope
      Tuple tt(0, 0, 0);
      for (unsigned ii=0; ii<gid.size(); ++ii)
	 tt += indexToTuple(gid[ii]);
      
      center_[0] = tt.x()/double(gid.size());
      center_[1] = tt.y()/double(gid.size());
      center_[2] = tt.z()/double(gid.size());
   }
   

   THREE_VECTOR center;
   center.x = center_[0];
   center.y = center_[1];
   center.z = center_[2];
   
   IndexToVector indexToVector(nx, ny, nz);
   
   for (unsigned ii=0; ii<gid.size(); ++ii)
   {
      THREE_VECTOR vv = indexToVector(gid[ii]);
      double rSq = DIFFSQ(vv, center);
      radius_ = max(radius_, rSq);
   }
   radius_ = sqrt(radius_);
}

