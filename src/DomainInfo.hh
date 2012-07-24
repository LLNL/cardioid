#ifndef DOMAIN_INFO_HH
#define DOMAIN_INFO_HH

#include <vector>
#include "Long64.hh"
#include "Vector.hh"

class DomainInfo
{
 public:
    DomainInfo(const std::vector<Long64>& gid, int nx, int ny, int nz);
   
    const Vector& center() const {return center_;}
    double radius() const {return radius_;}
    int ncells() const {return ncells_;}
private:

    int ncells_;
    double radius_;
    Vector center_;
};

#endif
