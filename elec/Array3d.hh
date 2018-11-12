#ifndef ARRAY3D_HH
#define ARRAY3D_HH

#include <vector>
#include <cassert>
#include <cstdlib>
#include <inttypes.h>
#include <iostream>

#include "AlignedAllocator.hh"
template <class T>
class Array3d
{
 public:

   Array3d();
   Array3d(unsigned nx, unsigned ny, unsigned nz);
   Array3d(unsigned nx, unsigned ny, unsigned nz, const T& initValue);
   
   T& operator()(unsigned index);
   const T& operator()(unsigned index) const;
   T& operator()(unsigned ix, unsigned iy, unsigned iz);
   const T& operator()(unsigned ix, unsigned iy, unsigned iz) const;

   T*** cArray();
   const T*** cArray() const;
   T* cBlock() {return &(data_[0]);};
   const T* cBlock() const {return &(data_[0]);};

   
   unsigned size() const;
   unsigned nx() const;
   unsigned ny() const;
   unsigned nz() const;
   
   void resize(unsigned nx, unsigned ny, unsigned nz);
   void resize(unsigned nx, unsigned ny, unsigned nz, const T& initValue);
   
   unsigned const tupleToIndex(unsigned ix, unsigned iy, unsigned iz) const;
   
   void dump(unsigned nx, unsigned ny, unsigned nz);

 private:
   void pointerSetup();
   
   unsigned nx_;
   unsigned ny_;
   unsigned nz_;
   std::vector<T, AlignedAllocator<T> > data_;
   T*** cArray_;
};

template <class T> 
Array3d<T>::Array3d()
: nx_(0), ny_(0), nz_(0), data_(0), cArray_(0)
{
}

template <class T> 
Array3d<T>::Array3d(unsigned nx, unsigned ny, unsigned nz)
: nx_(nx), ny_(ny), nz_(nz), data_(nx*ny*nz)
{
   pointerSetup();
}

template <class T> 
Array3d<T>::Array3d(unsigned nx, unsigned ny, unsigned nz, const T& initValue)
: nx_(nx), ny_(ny), nz_(nz), data_(nx*ny*nz, initValue)
{
   pointerSetup();
}


template <class T> inline
T& Array3d<T>::operator()(unsigned index)
{
   return data_[index];
}

template <class T> inline
const T& Array3d<T>::operator()(unsigned index) const
{
   return data_[index];
}

template <class T> inline
T& Array3d<T>::operator()(unsigned ix, unsigned iy, unsigned iz)
{
   return data_[tupleToIndex(ix, iy, iz)];
}

template <class T> inline
const T& Array3d<T>::operator()(unsigned ix, unsigned iy, unsigned iz) const
{
   return data_[tupleToIndex(ix, iy, iz)];
}

template <class T> inline
T*** Array3d<T>::cArray()
{
   return cArray_;
}

template <class T> inline
const T*** Array3d<T>::cArray() const
{
   return cArray_;
}

template <class T> inline
unsigned Array3d<T>::size() const
{
//   return nx_*ny_*nz_;
   return data_.size();
}

template <class T> inline
unsigned Array3d<T>::nx() const
{
   return nx_;
}

template <class T> inline
unsigned Array3d<T>::ny() const
{
   return ny_;
}

template <class T> inline
unsigned Array3d<T>::nz() const
{
   return nz_;
}


template <class T>  
void Array3d<T>::dump(unsigned nx, unsigned ny, unsigned nz)
{
   for(int ii=0;ii<nx;ii++)
   {
   for(int jj=0;jj<ny;jj++)
   {
   for(int kk=0;kk<nz;kk++)
     std::cout << (*this)(ii,jj,kk) << " ";
   std::cout << std::endl;
   }
   }
}

template <class T>  
void Array3d<T>::resize(unsigned nx, unsigned ny, unsigned nz)
{
   assert(nx_ == 0 && ny_ == 0 && nz_ == 0);
   data_.resize(nx*ny*nz);
   nx_ = nx;
   ny_ = ny;
   nz_ = nz;
   pointerSetup();
}
                        
template <class T>  
void Array3d<T>::resize(unsigned nx, unsigned ny, unsigned nz, const T& initValue)
{
   assert(nx_ == 0 && ny_ == 0 && nz_ == 0);
   data_.resize(nx*ny*nz, initValue);
   nx_ = nx;
   ny_ = ny;
   nz_ = nz;
   pointerSetup();
}

template <class T> inline 
unsigned const Array3d<T>::tupleToIndex(unsigned ix, unsigned iy, unsigned iz) const
{
   return iz + nz_*(iy + ny_*(ix));
//   return ix + nx_*(iy + ny_*(iz));
}

template <class T> 
void Array3d<T>::pointerSetup()
{
   cArray_ = (T***) malloc(sizeof(T**)*nx_);
   for (unsigned ii=0; ii<nx_; ++ii)
   {
      cArray_[ii] = (T**) malloc(sizeof(T*)*ny_);
      for (unsigned jj=0; jj<ny_; ++jj)
         cArray_[ii][jj] = &data_[0] + nz_*(jj + ny_*(ii)); 
   }
}


#endif

