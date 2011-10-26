#ifndef DOMAIN_INFO_HH
#define DOMAIN_INFO_HH

#include <vector>
#include "Long64.hh"

class DomainInfo
{
 public:
   DomainInfo(const std::vector<Long64>& gid, int nx, int ny, int nz);
   
   const std::vector<double>& center() const {return center_;}
   double radius() const {return radius_;}
 private:
   double radius_;
   std::vector<double> center_;
};

#endif
