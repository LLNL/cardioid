#ifndef ANATOMY_HH
#define ANATOMY_HH

#include <vector>
#include <cassert>
#include "AnatomyCell.hh"
#include "Tuple.hh"
#include "SymmetricTensor.hh"
#include "IndexToTuple.hh"
#include "three_algebra.h"

class Conductivity;

class Anatomy
{
 public:

   Anatomy()
   : i2t_(0, 0, 0), nRemote_(0){};

   unsigned size() const;
   unsigned nLocal() const;
   unsigned nRemote() const;
   unsigned& nRemote();
   unsigned nGlobal() const;
   unsigned& nGlobal();

//    unsigned nx();
//    unsigned ny();
//    unsigned nz();

   const unsigned& nx() const;
   const unsigned& ny() const;
   const unsigned& nz() const;


   void setGridSize(int nx, int ny, int nz);

   double& dx();
   double& dy();
   double& dz();

   double dx() const;
   double dy() const;
   double dz() const;

   Long64 gid(unsigned ii) const;
   SymmetricTensor conductivity(unsigned ii) const;
   int cellType(unsigned ii) const;
   Tuple globalTuple(unsigned ii) const;

   std::vector<AnatomyCell>& cellArray();
   const std::vector<AnatomyCell>& cellArray() const;

   double& offset_x();
   double& offset_y();
   double& offset_z();

   THREE_VECTOR pointFromGid(unsigned ii) const;

 private:
   unsigned nx_, ny_, nz_;
   double dx_, dy_, dz_;
   double offset_x_, offset_y_, offset_z_;
   unsigned nRemote_;
   unsigned nGlobal_;
   IndexToTuple i2t_;

   std::vector<AnatomyCell> cell_;

};

inline unsigned  Anatomy::size() const { return cell_.size();}
inline unsigned  Anatomy::nLocal() const { return cell_.size()-nRemote_;}
inline unsigned  Anatomy::nRemote() const { return nRemote_;}
inline unsigned& Anatomy::nRemote()       { return nRemote_;}
inline unsigned  Anatomy::nGlobal() const { return nGlobal_;}
inline unsigned& Anatomy::nGlobal()       { return nGlobal_;}

// inline unsigned  Anatomy::nx() { return nx_;}
// inline unsigned  Anatomy::ny() { return ny_;}
// inline unsigned  Anatomy::nz() { return nz_;}

inline const unsigned&  Anatomy::nx() const { return nx_;}
inline const unsigned&  Anatomy::ny() const { return ny_;}
inline const unsigned&  Anatomy::nz() const { return nz_;}

inline void Anatomy::setGridSize(int nx, int ny, int nz)
{
   nx_ = nx; ny_ = ny; nz_ = nz; i2t_ = IndexToTuple(nx, ny, nz);
}

inline double&  Anatomy::dx()  { return dx_;}
inline double&  Anatomy::dy()  { return dy_;}
inline double&  Anatomy::dz()  { return dz_;}

inline double  Anatomy::dx() const { return dx_;}
inline double  Anatomy::dy() const { return dy_;}
inline double  Anatomy::dz() const { return dz_;}

inline Long64  Anatomy::gid(unsigned ii) const { return cell_[ii].gid_;}
inline SymmetricTensor Anatomy::conductivity(unsigned ii) const {return cell_[ii].sigma_;}
inline int  Anatomy::cellType(unsigned ii) const { return cell_[ii].cellType_;}

inline Tuple Anatomy::globalTuple(unsigned ii) const
{
   return i2t_(cell_[ii].gid_);
}

inline       std::vector<AnatomyCell>& Anatomy::cellArray()       {return cell_;}
inline const std::vector<AnatomyCell>& Anatomy::cellArray() const {return cell_;}

inline double& Anatomy::offset_x() {return offset_x_;}
inline double& Anatomy::offset_y() {return offset_y_;}
inline double& Anatomy::offset_z() {return offset_z_;}

inline THREE_VECTOR Anatomy::pointFromGid(unsigned ii) const
{
   int x=ii%nx_;
   int y=(ii/nx_) %ny_;
   int z=ii/nx_/ny_;

   // We retrive the center of the cell as a 3D point
   double xcoor=(x + 0.5)*dx_ + offset_x_;
   double ycoor=(y + 0.5)*dy_ + offset_y_;
   double zcoor=(z + 0.5)*dz_ + offset_z_;

   THREE_VECTOR point = { xcoor, ycoor, zcoor };

   return point;
}

#endif
