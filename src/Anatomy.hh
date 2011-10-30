#ifndef ANATOMY_HH
#define ANATOMY_HH

#include <vector>
#include <cassert>
#include "AnatomyCell.hh"
#include "Tuple.hh"
#include "IndexToTuple.hh"
class Anatomy
{
 public:
   unsigned size() const;
   unsigned nLocal() const;
   unsigned& nLocal();
   unsigned nRemote() const;
   unsigned& nRemote();

   unsigned& nx();
   unsigned& ny();
   unsigned& nz();

   unsigned nx() const;
   unsigned ny() const;
   unsigned nz() const;

   double& dx();
   double& dy();
   double& dz();

   double dx() const;
   double dy() const;
   double dz() const;
   
   int gid(unsigned ii) const;
   int theta(unsigned ii) const;
   int phi(unsigned ii) const;
   int cellType(unsigned ii) const;
   Tuple globalTuple(unsigned ii) const;

   std::vector<AnatomyCell>& cellArray();
   const std::vector<AnatomyCell>& cellArray() const;
   
   
 private:
   unsigned nx_, ny_, nz_;
   double dx_, dy_, dz_;
   unsigned nLocal_;
   unsigned nRemote_;
   std::vector<AnatomyCell> cell_; 

};

inline unsigned  Anatomy::size() const { return cell_.size();}
inline unsigned  Anatomy::nLocal() const { return nLocal_;}
inline unsigned& Anatomy::nLocal()       { return nLocal_;}
inline unsigned  Anatomy::nRemote() const { return nRemote_;}
inline unsigned& Anatomy::nRemote()       { return nRemote_;}

inline unsigned&  Anatomy::nx() { return nx_;}
inline unsigned&  Anatomy::ny() { return ny_;}
inline unsigned&  Anatomy::nz() { return nz_;}

inline unsigned  Anatomy::nx() const { return nx_;}
inline unsigned  Anatomy::ny() const { return ny_;}
inline unsigned  Anatomy::nz() const { return nz_;}

inline double&  Anatomy::dx()  { return dx_;}
inline double&  Anatomy::dy()  { return dy_;}
inline double&  Anatomy::dz()  { return dz_;}

inline double  Anatomy::dx() const { return dx_;}
inline double  Anatomy::dy() const { return dy_;}
inline double  Anatomy::dz() const { return dz_;}

inline int  Anatomy::gid(unsigned ii) const { return cell_[ii].gid_;}
inline int  Anatomy::theta(unsigned ii) const { return cell_[ii].theta_;}
inline int  Anatomy::phi(unsigned ii) const { return cell_[ii].phi_;}
inline int  Anatomy::cellType(unsigned ii) const { return cell_[ii].cellType_;}

inline Tuple Anatomy::globalTuple(unsigned ii) const
{
   IndexToTuple i2t(nx_, ny_, nz_);
   return i2t(cell_[ii].gid_);
}

inline       std::vector<AnatomyCell>& Anatomy::cellArray()       {return cell_;}
inline const std::vector<AnatomyCell>& Anatomy::cellArray() const {return cell_;}

#endif
