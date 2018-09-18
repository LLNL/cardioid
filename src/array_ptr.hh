#pragma once

#include <cstddef>
#include <iostream>
#include <cassert>
#include <vector>

#if defined(USE_CUDA) and defined(__NVCC__) 
#define HOST_DEVICE __device__ __host__
#define HOST __host__
#else
#define HOST_DEVICE
#define HOST
#endif

template <typename TTT>
class OnlyAssignable
{
 public:
   HOST_DEVICE inline OnlyAssignable(TTT& dataPoint) : dataPoint_(dataPoint) {}
   HOST_DEVICE inline TTT& operator=(const TTT value)
   {
      dataPoint_ = value;
      return dataPoint_;
   }
 private:
   TTT& dataPoint_;
};

template <typename TTT>
class AssignableIterator
{
 public:
   typedef std::output_iterator_tag iterator_category;
   typedef TTT value_type;
   typedef std::ptrdiff_t difference_type;
   typedef TTT* pointer;
   typedef TTT& reference;

   HOST_DEVICE inline AssignableIterator() : ptr_(NULL) {}
   HOST_DEVICE inline AssignableIterator(TTT* ptr) : ptr_(ptr) {}
   
   HOST_DEVICE inline OnlyAssignable<TTT> operator*() { return OnlyAssignable<TTT>(*ptr_); }
   HOST_DEVICE inline AssignableIterator<TTT>& operator++() { ++ptr_; return *this; }
   HOST_DEVICE inline AssignableIterator<TTT> operator++(int) { AssignableIterator<TTT> temp(*this); ptr_++; return temp; }
   HOST_DEVICE inline AssignableIterator<TTT>& operator--() { --ptr_; return *this; }
   HOST_DEVICE inline AssignableIterator<TTT> operator--(int) { AssignableIterator<TTT> temp(*this); ptr_--; return temp; }

   HOST_DEVICE inline bool operator==(const AssignableIterator<TTT> other) const { return ptr_==other.ptr_; }
   HOST_DEVICE inline bool operator!=(const AssignableIterator<TTT> other) const { return ptr_!=other.ptr_; }

 private:
   pointer ptr_;
};

template <typename External,typename _value_type>
class basic_impl
{
 public:
   typedef _value_type value_type;
   HOST_DEVICE inline External slice(const int begin, const int end)
   {
      value_type* newData = data_;
      if (newData != NULL) { newData += begin; }
      return External(newData, end-begin);
   }
   HOST_DEVICE inline std::size_t size() const { return size_; }
 protected:
   HOST_DEVICE inline value_type* deref() { return data_; }
   HOST_DEVICE void create(value_type* newData, const std::size_t newSize)
   {
      data_ = newData;
      size_ = newSize;
   }
   value_type* data_;
   std::size_t size_;
};

template <typename Accessor, typename assign_type, typename iterator>
class ro_array_ptr_interface : public Accessor
{
 public:
   typedef typename Accessor::value_type value_type;
   HOST_DEVICE inline assign_type operator[](const int ii) { return assign_type(this->deref()[ii]); }
   HOST_DEVICE inline iterator begin() { return iterator(this->deref()); }
   HOST_DEVICE inline iterator end() { return iterator(this->deref()+this->size()); }
   HOST_DEVICE inline value_type* raw() { return this->deref(); }   
};

template <typename Accessor, typename assign_type, typename iterator>
class rw_array_ptr_interface : public ro_array_ptr_interface<Accessor,assign_type,iterator>
{
 public:
   typedef typename Accessor::value_type value_type;
   HOST_DEVICE inline void fill(const value_type value)
   {
      value_type* data=this->deref();
      for (std::size_t ii=0; ii<this->size(); ii++)
      {
         data[ii] = value;
      }
   }
   HOST_DEVICE inline void assign(const std::size_t numelem, const value_type value)
   {
      assert(numelem == this->size());
      fill(value);
   }
};

template <typename TTT>
class ro_array_ptr : public ro_array_ptr_interface<basic_impl<ro_array_ptr<TTT>,const TTT>,const TTT&,const TTT*>
{
 public:
   HOST_DEVICE ro_array_ptr() : ro_array_ptr(NULL, 0) {}
   template <typename Allocator>
   HOST_DEVICE ro_array_ptr(const std::vector<TTT,Allocator>& newData) : ro_array_ptr(&newData[0], newData.size()) {}
   HOST_DEVICE ro_array_ptr(const TTT* newData, const std::size_t newSize) { this->create(newData,newSize); }
};

template <typename TTT>
class wo_array_ptr : public rw_array_ptr_interface<basic_impl<wo_array_ptr<TTT>,TTT>,OnlyAssignable<TTT>,AssignableIterator<TTT> >
{
 public:
   HOST_DEVICE wo_array_ptr() : wo_array_ptr(NULL, 0) {}
   template <typename Allocator>
   HOST_DEVICE wo_array_ptr(std::vector<TTT,Allocator>& newData) : wo_array_ptr(&newData[0], newData.size()) {}
   HOST_DEVICE wo_array_ptr(TTT* newData, const std::size_t newSize) { this->create(newData,newSize); }
};

template <typename TTT>
class rw_array_ptr : public rw_array_ptr_interface<basic_impl<rw_array_ptr<TTT>,TTT>,TTT&,TTT*>
{
 public:
   HOST_DEVICE rw_array_ptr() : rw_array_ptr(NULL, 0) {}
   template <typename Allocator>
   HOST_DEVICE rw_array_ptr(std::vector<TTT,Allocator>& newData) : rw_array_ptr(&newData[0], newData.size()) {}
   HOST_DEVICE rw_array_ptr(TTT* newData, const std::size_t newSize) { this->create(newData,newSize); }
   HOST_DEVICE inline operator ro_array_ptr<TTT>()
   {
      return ro_array_ptr<TTT>(this->data_, this->size());
   }
   HOST_DEVICE inline operator wo_array_ptr<TTT>()
   {
      return wo_array_ptr<TTT>(this->data_, this->size());
   }
};
