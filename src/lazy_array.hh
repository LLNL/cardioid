#pragma once

#include <cstddef>
#include <cassert>
#include <iostream>
#include "SpaceAllocator.hh"
#include "ContextManager.hh"

#if defined(USE_CUDA) and defined(__NVCC__) 
#define HOST_DEVICE __device__ __host__
#define HOST __host__
#else
#define HOST_DEVICE
#define HOST
#endif

template <typename TTT>
class lazy_array;

template <typename TTT>
class ro_larray_ptr
{
  public:
   typedef const TTT* iterator;
   
   HOST_DEVICE inline std::size_t size() const { return size_; }
   HOST_DEVICE inline const TTT& operator[](const int ii) const { return activePointer_[ii]; }
   HOST_DEVICE inline const TTT* begin() const { return activePointer_; }
   HOST_DEVICE inline const TTT* end() const { return activePointer_+size(); }
   HOST_DEVICE inline const TTT* raw() const { return activePointer_; }

   HOST_DEVICE inline ro_larray_ptr(const ro_larray_ptr<TTT>& other)
   : ro_larray_ptr(other.dataManager_, other.activePointer_, other.offset_, other.size_) {}
   HOST_DEVICE inline ro_larray_ptr() : ro_larray_ptr(NULL, NULL, 0, 0) {}
   HOST_DEVICE inline ro_larray_ptr(lazy_array<TTT>* dataManager, const TTT* activePointer, std::ptrdiff_t offset, const std::size_t newSize)
   : dataManager_(dataManager), activePointer_(activePointer), offset_(offset), size_(newSize) {
#if !defined(__CUDA_ARCH__)
      use();
#endif
   }
   HOST_DEVICE inline ro_larray_ptr(const TTT* rawData, std::size_t newSize)
   : ro_larray_ptr(NULL, rawData, 0, newSize) {}
   template <typename Allocator>
   HOST_DEVICE inline ro_larray_ptr(const std::vector<TTT,Allocator>& rawData)
   : ro_larray_ptr(NULL, &rawData[0], 0, rawData.size()) {}
      
   HOST_DEVICE inline ro_larray_ptr<TTT> slice(const int begin, const int end)
   {
      const TTT* newActivePointer = (activePointer_ == NULL) ? NULL : activePointer_+begin;
      return ro_larray_ptr<TTT>(dataManager_, newActivePointer, (activePointer_+begin) - activePointer_, end-begin);
   }

   HOST void use()
   {
      if (dataManager_ != NULL)
      {
         ExecutionSpace space = dataManager_->getContext();
         if (space >= 0) { useOn(space); }
      }
   }
   
   HOST void useNowhere()
   {
      activePointer_ = NULL;
   }
   
  private:
   HOST void useOn(const ExecutionSpace space)
   {
      activePointer_ = dataManager_->readonlyRaw(space,offset_);
   }

   const TTT* activePointer_;
   std::ptrdiff_t offset_;
   std::size_t size_;
   lazy_array<TTT>* dataManager_;
};

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

template <typename TTT>
class wo_larray_ptr
{
  public:
   typedef AssignableIterator<TTT> iterator;
   
   HOST_DEVICE inline std::size_t size() const { return size_; }
   HOST_DEVICE inline OnlyAssignable<TTT> operator[](const int ii) { return OnlyAssignable<TTT>(activePointer_[ii]); }
   HOST_DEVICE inline iterator begin() { return iterator(activePointer_); }
   HOST_DEVICE inline iterator end() { return iterator(activePointer_+size()); }
   HOST_DEVICE inline TTT* raw() { return activePointer_; }

   HOST_DEVICE inline void fill(const TTT value)
   {
      TTT* data=activePointer_;
      for (std::size_t ii=0; ii<size(); ii++)
      {
         data[ii] = value;
      }
   }
   HOST_DEVICE inline void assign(const std::size_t numelem, const TTT value)
   {
      assert(numelem == size());
      fill(value);
   }

   HOST_DEVICE inline wo_larray_ptr(const wo_larray_ptr<TTT>& other)
   : wo_larray_ptr(other.dataManager_, other.activePointer_, other.offset_, other.size_) {}
   HOST_DEVICE inline wo_larray_ptr() : wo_larray_ptr(NULL, NULL, 0, 0) {}
   HOST_DEVICE inline wo_larray_ptr(lazy_array<TTT>* dataManager, TTT* activePointer, std::ptrdiff_t offset, const std::size_t newSize)
   : dataManager_(dataManager), activePointer_(activePointer), offset_(offset), size_(newSize)
   {
#if !defined(__CUDA_ARCH__)
      use();
#endif
   }
   HOST_DEVICE inline wo_larray_ptr(TTT* rawData, const std::size_t newSize)
   : wo_larray_ptr(NULL, rawData, 0, newSize) {}
   HOST_DEVICE inline wo_larray_ptr(std::vector<TTT>& rawData)
   : wo_larray_ptr(NULL, &rawData[0], 0, rawData.size()) {}


   HOST_DEVICE inline wo_larray_ptr<TTT> slice(const int begin, const int end)
   {
      TTT* newActivePointer = (activePointer_ == NULL) ? NULL : activePointer_+begin;
      return wo_larray_ptr<TTT>(dataManager_, newActivePointer, (activePointer_+begin) - activePointer_, end-begin);
   }

   HOST void use()
   {
      if (dataManager_ != NULL)
      {
         ExecutionSpace space = dataManager_->getContext();
         if (space >= 0) { useOn(space); }
      }
   }

   HOST void useNowhere()
   {
      activePointer_ = NULL;
   }
   
  private:
   HOST void useOn(const ExecutionSpace space)
   {
      activePointer_ = dataManager_->writeonlyRaw(space,offset_);
   }

   TTT* activePointer_;
   std::ptrdiff_t offset_;
   std::size_t size_;
   lazy_array<TTT>* dataManager_;
};

template <typename TTT>
class rw_larray_ptr
{
  public:
   typedef TTT* iterator;
   
   HOST_DEVICE inline std::size_t size() const { return size_; }
   HOST_DEVICE inline TTT& operator[](const int ii) { return activePointer_[ii]; }
   HOST_DEVICE inline TTT* begin() { return activePointer_; }
   HOST_DEVICE inline TTT* end() { return activePointer_+size(); }
   HOST_DEVICE inline TTT* raw() { return activePointer_; }

   HOST_DEVICE inline void fill(const TTT value)
   {
      TTT* data=activePointer_;
      for (std::size_t ii=0; ii<size(); ii++)
      {
         data[ii] = value;
      }
   }
   HOST_DEVICE inline void assign(const std::size_t numelem, const TTT value)
   {
      assert(numelem == size());
      fill(value);
   }

   HOST_DEVICE inline operator ro_larray_ptr<TTT>()
   {
      return ro_larray_ptr<TTT>(dataManager_, activePointer_, offset_, size_);
   }
   HOST_DEVICE inline operator wo_larray_ptr<TTT>()
   {
      return wo_larray_ptr<TTT>(dataManager_, activePointer_, offset_, size_);
   }


   HOST_DEVICE inline rw_larray_ptr(const rw_larray_ptr<TTT>& other)
   : rw_larray_ptr(other.dataManager_, other.activePointer_, other.offset_, other.size_) {}
   HOST_DEVICE inline rw_larray_ptr() : rw_larray_ptr(NULL, NULL, 0, 0) {}
   HOST_DEVICE inline rw_larray_ptr(lazy_array<TTT>* dataManager, TTT* activePointer, std::ptrdiff_t offset, const std::size_t newSize)
   : dataManager_(dataManager), activePointer_(activePointer), offset_(offset), size_(newSize)
   {
#if !defined(__CUDA_ARCH__)
      use();
#endif
   }
   HOST_DEVICE inline rw_larray_ptr(TTT* rawData, const std::size_t newSize)
   : rw_larray_ptr(NULL, rawData, 0, newSize) {}
   HOST_DEVICE inline rw_larray_ptr(std::vector<TTT>& rawData)
   : rw_larray_ptr(NULL, &rawData[0], 0, rawData.size()) {}

   HOST_DEVICE inline rw_larray_ptr<TTT> slice(const int begin, const int end)
   {
      TTT* newActivePointer = (activePointer_ == NULL) ? NULL : activePointer_+begin;
      return rw_larray_ptr<TTT>(dataManager_, newActivePointer, (activePointer_+begin) - activePointer_, end-begin);
   }

   HOST void use()
   {
      if (dataManager_ != NULL)
      {
         ExecutionSpace space = dataManager_->getContext();
         if (space >= 0) { useOn(space); }
      }
   }

   HOST void useNowhere()
   {
      activePointer_ = NULL;
   }
   
 private:
   HOST void useOn(const ExecutionSpace space)
   {
      activePointer_ = dataManager_->readwriteRaw(space,offset_);
   }

   TTT* activePointer_;
   std::ptrdiff_t offset_;
   std::size_t size_;
   lazy_array<TTT>* dataManager_;
};


template <typename TTT>
class lazy_array
{
 public:

   lazy_array() {
      for (int ii=0; ii<NUMSPACES; ii++)
      {
         pointerRecord_[ii] = NULL;
      }
      clear();
      contextManager_ = getContextManager();
   }
   ~lazy_array()
   {
      clear();
   }
   template <typename Allocator>
   lazy_array(const std::vector<TTT, Allocator>& vvv) : lazy_array()
   {
      resize(vvv.size());
      ContextRegion region(CPU);
      wo_larray_ptr<TTT> myData = writeonly();
      for (unsigned int ii=0; ii<vvv.size(); ii++)
      {
         myData[ii] = vvv[ii];
      }
   }

   
   lazy_array(const lazy_array<TTT>& other)
   : size_(other.size_)
   {
      for (int ispace=0; ispace<NUMSPACES; ispace++)
      {
         isValid_[ispace] = other.isValid_[ispace];
         if (other.pointerRecord_[ispace] == NULL) { pointerRecord_[ispace] = NULL; }
         else
         {
            spaceMalloc(ispace,&pointerRecord_[ispace],size());
            spaceMemcpy(ispace,pointerRecord_[ispace],
                        ispace,other.pointerRecord_[ispace],
                        size());
         }
      }
      contextManager_ = other.contextManager_;
   }   
   friend inline void swap(lazy_array<TTT>& first, lazy_array<TTT>& second)
   {
      using std::swap;

      swap(first.size_,second.size_);
      for (int ispace=0; ispace<NUMSPACES; ispace++)
      {
         swap(first.pointerRecord_[ispace],second.pointerRecords[ispace]);
         swap(first.isValid_[ispace],second.isValid_[ispace]);
      }
      swap(first.contextManager_, second.contextManager_);
   }
   lazy_array<TTT>& operator=(lazy_array<TTT> other)
   {
      swap(*this,other);
      return *this;
   }
   
   HOST inline bool empty() const { return size_==0; }
   HOST inline std::size_t size() const { return size_; }

   HOST void clear() {
      for (int ispace=0; ispace<NUMSPACES; ispace++)
      {
         if (pointerRecord_[ispace] != NULL) { spaceFree(ispace,pointerRecord_[ispace]); }
         pointerRecord_[ispace] = NULL;
         isValid_[ispace] = false;
      }
      size_ = 0;
   }
   
   HOST void resize(const std::size_t elems)
   {
      if (elems != size())
      {
         std::size_t smallest = (elems < size()) ? elems : size();
         for (int ispace=0; ispace<NUMSPACES; ispace++)
         {
            if (pointerRecord_[ispace] != NULL)
            {
               TTT* newPointer;
               spaceMalloc(ispace, &newPointer, elems);
               spaceMemcpy(ispace, newPointer, ispace, pointerRecord_[ispace], smallest);
               spaceFree(ispace, pointerRecord_[ispace]);
               pointerRecord_[ispace] = newPointer;
            }
         }
      }
      size_ = elems;
   }

   inline ExecutionSpace getContext() const { return contextManager_->current(); }
   
   inline operator rw_larray_ptr<TTT>() { return readwrite(); }
   inline operator wo_larray_ptr<TTT>() { return writeonly(); } 
   inline operator ro_larray_ptr<TTT>() const { return readonly(); }

   inline rw_larray_ptr<TTT> readwrite() { return rw_larray_ptr<TTT>(this, NULL, 0, size()); }
   inline wo_larray_ptr<TTT> writeonly() { return wo_larray_ptr<TTT>(this, NULL, 0, size()); }
   inline ro_larray_ptr<TTT> readonly() const { return ro_larray_ptr<TTT>(const_cast<lazy_array<TTT>*>(this), NULL, 0, size()); }
   
   TTT* readwriteRaw(ExecutionSpace space,std::ptrdiff_t offset=0)
   {
      allocateFor(space);
      moveTo(space);
      writeOn(space);
      
      return lookup(space,offset);
   }
   TTT* writeonlyRaw(ExecutionSpace space, std::ptrdiff_t offset=0)
   {
      allocateFor(space);
      writeOn(space);

      return lookup(space,offset);
   }
   TTT* readonlyRaw(ExecutionSpace space, std::ptrdiff_t offset=0)
   {
      allocateFor(space);
      moveTo(space);

      return lookup(space,offset);
   }   
   
 private:
   TTT* pointerRecord_[NUMSPACES];
   bool isValid_[NUMSPACES];
   std::size_t size_;
   ContextManager* contextManager_;

   inline void allocateFor(ExecutionSpace space)
   {
      if (pointerRecord_[space] == NULL)
      {
         spaceMalloc(space,&pointerRecord_[space],size());
         isValid_[space] = false;         
      }
   }
   inline void writeOn(ExecutionSpace space)
   {
      for (int ispace=0; ispace<NUMSPACES; ispace++)
      {
         isValid_[ispace] = false;
      }
      isValid_[space] = true;
   }
   inline void moveTo(ExecutionSpace space)
   {
      if (! isValid_[space])
      {
          
         if (NUMSPACES == 1) {}
         else if (NUMSPACES == 2) {
            ExecutionSpace otherSpace;
            if (space == GPU) { otherSpace = CPU; }
            else { otherSpace = GPU; }
            if (isValid_[otherSpace])
            {
               spaceMemcpy(space,pointerRecord_[space],
                           otherSpace, pointerRecord_[otherSpace],
                           size());
            }
         } else {
            assert(false && "Need to figure out an algorithm for more than 2 spaces");
         }
         isValid_[space] = true;
      }
   }
   inline TTT* lookup(ExecutionSpace space, std::ptrdiff_t offset)
   {
      return pointerRecord_[space]+offset;
   }
   
};

