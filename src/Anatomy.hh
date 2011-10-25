#ifndef ANATOMY_HH
#define ANATOMY_HH

#include <vector>
#include <cassert>
#include "AnatomyCell.hh"
#include "Tuple.hh"

class Anatomy
{
 public:
   unsigned size() const;
   unsigned nLocal() const;
   

   double dx() const;
   double dy() const;
   double dz() const;
   

   
   int theta(unsigned ii) const;
   int phi(unsigned ii) const;
   int cellType(unsigned ii) const;
   Tuple globalTuple(unsigned ii) const;

   std::vector<AnatomyCell>& cellArray();
   const std::vector<AnatomyCell>& cellArray() const;
   
   
 private:
   unsigned nCellLocal_;
   unsigned nCellRemote_;
   std::vector<AnatomyCell> cell_; 

};

inline unsigned Anatomy::size() const { assert (1==0); return 0;}
inline unsigned Anatomy::nLocal() const { assert (1==0); return 0;}

inline double  Anatomy::dx() const { assert (1==0); return 0;}
inline double  Anatomy::dy() const { assert (1==0); return 0;}
inline double  Anatomy::dz() const { assert (1==0); return 0;}

inline int  Anatomy::theta(unsigned ii) const { assert (1==0); return 0;}
inline int  Anatomy::phi(unsigned ii) const { assert (1==0); return 0;}
inline int  Anatomy::cellType(unsigned ii) const { assert (1==0); return 0;}

inline Tuple Anatomy::globalTuple(unsigned ii) const { assert (1==0); return Tuple(0, 0, 0);}

inline       std::vector<AnatomyCell>& Anatomy::cellArray()       {return cell_;}
inline const std::vector<AnatomyCell>& Anatomy::cellArray() const {return cell_;}

#endif
