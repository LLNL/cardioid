#pragma once

#include <vector>
#include <cassert>
#include "Ledger.hh"


#ifdef NVCC
#define HOST_DEVICE __host__ __device__
#else
#define HOST_DEVICE
#endif

template <typename TTT>
class ArrayView
{
 public:
   HOST_DEVICE ArrayView() : data_(nullptr), size_(0) {}   
   HOST_DEVICE ArrayView(TTT * data, const std::size_t size) : data_(data), size_(size) {}

   template <typename AAA>
   ArrayView(std::vector<TTT,AAA>& vvv) : ArrayView(&vvv[0],vvv.size()) {}

   HOST_DEVICE std::size_t size() const { return size_; }
   HOST_DEVICE TTT& operator[](const int ii) { return data_[ii]; }
   HOST_DEVICE const TTT& operator[](const int ii) const { return data_[ii]; }

   //operator TTT* () { return data_; }
   //operator const TTT*() const { return data_; }

   HOST_DEVICE TTT* begin() { return data_; }
   HOST_DEVICE TTT* end() { return data_+size(); }

   HOST_DEVICE const TTT* cbegin() const { return data_; }
   HOST_DEVICE const TTT* cend() const { return data_+size(); }

   HOST_DEVICE const TTT* begin() const { return cbegin(); }
   HOST_DEVICE const TTT* end() const { return cend(); }

   HOST_DEVICE operator TTT*() { return data_; }
   HOST_DEVICE operator const TTT*() const { return data_; }
   
   HOST_DEVICE void fill(TTT value)
   {
      for (std::size_t ii=0; ii<size(); ii++)
      {
         data_[ii] = value;
      }
   }
   HOST_DEVICE void assign(int numelem, TTT value)
   {
      assert(numelem == size());
      fill(value);
   }

 protected:
   TTT * data_;
   std::size_t size_;
};


/*
  ArrayView(const SimpleArray<TTT>& vvv) : data_(vvv), size_(vvv.size()) {}
  ArrayView(const SimpleArray<TTT>* pv) : ArrayView(*pv) {}
*/


template <typename TTT>
class ConstArrayView : public ArrayView<const TTT>
{
 public:
   HOST_DEVICE ConstArrayView(const TTT * const data, const std::size_t size) : ArrayView<const TTT>(data, size) {};
   HOST_DEVICE ConstArrayView(ArrayView<TTT> other) : ConstArrayView(&other[0], other.size()) {}

   template <typename AAA>
   ConstArrayView(const std::vector<TTT,AAA>& vvv) : ConstArrayView(&vvv[0], vvv.size()) {}
};
