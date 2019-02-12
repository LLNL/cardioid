#pragma once

#include <cstddef>
#include <iostream>
#include <cassert>
#include <vector>

#if defined(USE_CUDA) and defined(__NVCC__) 
#define LAZY_HOST_DEVICE __device__ __host__
#define LAZY_HOST_ONLY __host__
#else
#define LAZY_HOST_DEVICE
#define LAZY_HOST_ONLY
#endif


//! Expression template object used for assignments.
template <typename TTT>
class OnlyAssignable
{
 public:
   LAZY_HOST_DEVICE inline OnlyAssignable(TTT& dataPoint) : dataPoint_(dataPoint) {}
   LAZY_HOST_DEVICE inline TTT& operator=(const TTT value)
   {
      dataPoint_ = value;
      return dataPoint_;
   }
 private:
   TTT& dataPoint_;
};

//! Iterator for write-only arrays.  Prevents entries from being read, but allows them to be written.
template <typename TTT>
class AssignableIterator
{
 public:
   typedef std::output_iterator_tag iterator_category;
   typedef TTT value_type;
   typedef std::ptrdiff_t difference_type;
   typedef TTT* pointer;
   typedef TTT& reference;

   LAZY_HOST_DEVICE inline AssignableIterator() : ptr_(NULL) {}
   LAZY_HOST_DEVICE inline AssignableIterator(TTT* ptr) : ptr_(ptr) {}
   
   LAZY_HOST_DEVICE inline OnlyAssignable<TTT> operator*() { return OnlyAssignable<TTT>(*ptr_); }
   LAZY_HOST_DEVICE inline AssignableIterator<TTT>& operator++() { ++ptr_; return *this; }
   LAZY_HOST_DEVICE inline AssignableIterator<TTT> operator++(int) { AssignableIterator<TTT> temp(*this); ptr_++; return temp; }
   LAZY_HOST_DEVICE inline AssignableIterator<TTT>& operator--() { --ptr_; return *this; }
   LAZY_HOST_DEVICE inline AssignableIterator<TTT> operator--(int) { AssignableIterator<TTT> temp(*this); ptr_--; return temp; }

   LAZY_HOST_DEVICE inline bool operator==(const AssignableIterator<TTT> other) const { return ptr_==other.ptr_; }
   LAZY_HOST_DEVICE inline bool operator!=(const AssignableIterator<TTT> other) const { return ptr_!=other.ptr_; }

 private:
   pointer ptr_;
};

/**
 * OK, this is evil, but I can't find another way to do it.
 * I'm using some unholy variation of CRTP here to ease the implementation.
 * I've tried multiple versions, and this is the easiest way to implement things
 * that avoids unnecessary code duplication.
 *
 * Basically this is a bridge design pattern using static
 * polymorphism.  *_interface is the abstraction layer, and *_impl is
 * the implementation layer.
 *
 * *_array_ptr_interface presents an interface that looks like vectors
 * except you can't change their size.  It implements things like
 * operator[], iterators, so on and so forth.  To do this, they use
 * size(), deref(), and create() on an _impl class.
 *
 * The basic_impl class does access to bare vectors, while the
 * lazy_array_impl class offers something closer to CHAI access.
 *
 * These get combined together to make their external types.
 */

template <typename External,typename _value_type>
class basic_impl
{
 public:
   typedef _value_type value_type;
   LAZY_HOST_DEVICE inline External slice(const int begin, const int end)
   {
      value_type* newData = data_;
      if (newData != NULL) { newData += begin; }
      return External(newData, end-begin);
   }
   LAZY_HOST_DEVICE inline std::size_t size() const { return size_; }
 protected:
   LAZY_HOST_DEVICE inline value_type* deref() const { return data_; }
   LAZY_HOST_DEVICE void create(value_type* newData, const std::size_t newSize)
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
   LAZY_HOST_DEVICE inline assign_type operator[](const int ii) const { return assign_type(this->deref()[ii]); }
   LAZY_HOST_DEVICE inline iterator begin() const { return iterator(this->deref()); }
   LAZY_HOST_DEVICE inline iterator end() const { return iterator(this->deref()+this->size()); }
   LAZY_HOST_DEVICE inline value_type* raw() const { return this->deref(); }   
};

template <typename Accessor, typename assign_type, typename iterator>
class rw_array_ptr_interface : public ro_array_ptr_interface<Accessor,assign_type,iterator>
{
 public:
   typedef typename Accessor::value_type value_type;
   LAZY_HOST_DEVICE inline void fill(const value_type value)
   {
      value_type* data=this->deref();
      for (std::size_t ii=0; ii<this->size(); ii++)
      {
         data[ii] = value;
      }
   }
   LAZY_HOST_DEVICE inline void assign(const std::size_t numelem, const value_type value)
   {
      assert(numelem == this->size());
      fill(value);
   }
};

template <typename TTT>
class ro_array_ptr : public ro_array_ptr_interface<basic_impl<ro_array_ptr<TTT>,const TTT>,const TTT&,const TTT*>
{
 public:
   LAZY_HOST_DEVICE ro_array_ptr() : ro_array_ptr(NULL, 0) {}
   template <typename Allocator>
   LAZY_HOST_ONLY ro_array_ptr(const std::vector<TTT,Allocator>& newData) : ro_array_ptr(&newData[0], newData.size()) {}
   LAZY_HOST_DEVICE ro_array_ptr(const TTT* newData, const std::size_t newSize) { this->create(newData,newSize); }
};

template <typename TTT>
class wo_array_ptr : public rw_array_ptr_interface<basic_impl<wo_array_ptr<TTT>,TTT>,OnlyAssignable<TTT>,AssignableIterator<TTT> >
{
 public:
   LAZY_HOST_DEVICE wo_array_ptr() : wo_array_ptr(NULL, 0) {}
   template <typename Allocator>
   LAZY_HOST_ONLY wo_array_ptr(std::vector<TTT,Allocator>& newData) : wo_array_ptr(&newData[0], newData.size()) {}
   LAZY_HOST_DEVICE wo_array_ptr(TTT* newData, const std::size_t newSize) { this->create(newData,newSize); }
};

template <typename TTT>
class rw_array_ptr : public rw_array_ptr_interface<basic_impl<rw_array_ptr<TTT>,TTT>,TTT&,TTT*>
{
 public:
   LAZY_HOST_DEVICE rw_array_ptr() : rw_array_ptr(NULL, 0) {}
   template <typename Allocator>
   LAZY_HOST_ONLY rw_array_ptr(std::vector<TTT,Allocator>& newData) : rw_array_ptr(&newData[0], newData.size()) {}
   LAZY_HOST_DEVICE rw_array_ptr(TTT* newData, const std::size_t newSize) { this->create(newData,newSize); }
   LAZY_HOST_DEVICE inline operator ro_array_ptr<TTT>()
   {
      return ro_array_ptr<TTT>(this->data_, this->size());
   }
   LAZY_HOST_DEVICE inline operator wo_array_ptr<TTT>()
   {
      return wo_array_ptr<TTT>(this->data_, this->size());
   }
};
