#pragma once

#include "array_ptr.hh"
#include "SpaceAllocator.hh"

template <typename TTT>
class lazy_array;

template< class T > struct my_remove_const          { typedef T type; };
template< class T > struct my_remove_const<const T> { typedef T type; };


template <typename External, typename Basic, typename TTT>
class managed_array_impl
{
 public:
   typedef TTT value_type;
   typedef typename my_remove_const<TTT>::type noconst_value_type;
   LAZY_HOST_DEVICE inline External slice(const int begin, const int end)
   {
      value_type* nullPointer = NULL;
      return External(dataManager_, (nullPointer+begin) - nullPointer, end-begin);
   }
   LAZY_HOST_DEVICE inline std::size_t size() const { return size_; }
   LAZY_HOST_ONLY inline Basic useOn(const ExecutionSpace space)
   {
      return Basic(this->dataManager_->raw(this, space, offset_), size_);
   }
 protected:
   LAZY_HOST_DEVICE inline void create(lazy_array<noconst_value_type>* dataManager, std::ptrdiff_t offset, const std::size_t newSize)
   {
      dataManager_ = dataManager;
      offset_ = offset;
      size_ = newSize;
   }
   std::ptrdiff_t offset_;
   std::size_t size_;
   lazy_array<noconst_value_type>* dataManager_;
};

template <typename TTT>
class ro_mgarray_ptr : public managed_array_impl<ro_mgarray_ptr<TTT>,ro_array_ptr<TTT>,const TTT>
{
 public:
   LAZY_HOST_DEVICE inline ro_mgarray_ptr() : ro_mgarray_ptr(NULL, 0, 0) {}
   LAZY_HOST_DEVICE inline ro_mgarray_ptr(const ro_mgarray_ptr<TTT>& other)
   : ro_mgarray_ptr(other.dataManager_, other.offset_, other.size_) {}
   LAZY_HOST_DEVICE inline ro_mgarray_ptr(lazy_array<TTT>* dataManager, const std::ptrdiff_t offset, const std::size_t newSize)
   { this->create(dataManager,offset,newSize); }
};

template <typename TTT>
class wo_mgarray_ptr : public managed_array_impl<wo_mgarray_ptr<TTT>,wo_array_ptr<TTT>,TTT>
{
 public:
   LAZY_HOST_DEVICE inline wo_mgarray_ptr() : wo_mgarray_ptr(NULL, 0, 0) {}
   LAZY_HOST_DEVICE inline wo_mgarray_ptr(const wo_mgarray_ptr<TTT>& other)
   : wo_mgarray_ptr(other.dataManager_, other.offset_, other.size_) {}
   LAZY_HOST_DEVICE inline wo_mgarray_ptr(lazy_array<TTT>* dataManager, const std::ptrdiff_t offset, const std::size_t newSize)
   { this->create(dataManager,offset,newSize); }
};

template <typename TTT>
class rw_mgarray_ptr : public managed_array_impl<rw_mgarray_ptr<TTT>,rw_array_ptr<TTT>,TTT>
{
 public:
   LAZY_HOST_DEVICE inline rw_mgarray_ptr() : rw_mgarray_ptr(NULL, 0, 0) {}
   LAZY_HOST_DEVICE inline rw_mgarray_ptr(const rw_mgarray_ptr<TTT>& other)
   : rw_mgarray_ptr(other.dataManager_, other.offset_, other.size_) {}
   LAZY_HOST_DEVICE inline rw_mgarray_ptr(lazy_array<TTT>* dataManager, const std::ptrdiff_t offset, const std::size_t newSize)
   { this->create(dataManager,offset,newSize); }
   LAZY_HOST_DEVICE inline operator ro_mgarray_ptr<TTT>()
   {
      return ro_mgarray_ptr<TTT>(this->dataManager_, this->offset_, this->size());
   }
   LAZY_HOST_DEVICE inline operator wo_mgarray_ptr<TTT>()
   {
      return wo_mgarray_ptr<TTT>(this->dataManager_, this->offset_, this->size());
   }
};

#if 0
//Version of lazy arrays that work more like CHAI, can be indexed directly

#include "ContextManager.hh"

template <typename External,typename TTT>
class lazy_array_impl
{
 public:
   typedef TTT value_type;
   typedef typename my_remove_const<TTT>::type noconst_value_type;
   LAZY_HOST_DEVICE inline External slice(const int begin, const int end)
   {
      value_type* newData = activePointer_;
      if (newData != NULL) { newData += begin; }
      return External(dataManager_, newData, (activePointer_+begin) - activePointer_, end-begin);
   }
   LAZY_HOST_ONLY inline void use()
   {
      if (dataManager_ != NULL)
      {
         ExecutionSpace space = dataManager_->getContext();
         if (space >= 0) { useOn(space); }
      }
   }
   LAZY_HOST_DEVICE inline std::size_t size() const { return size_; }
   
 protected:
   LAZY_HOST_DEVICE inline value_type* deref() const {
#if !defined(__CUDA_ARCH__) && defined(LAZY_ARRAY_CHECK_USAGE)
      assert(dataManager_->getContext() != CPU);
#endif
      return activePointer_;
   }
   LAZY_HOST_DEVICE inline void create(lazy_array<noconst_value_type>* dataManager, TTT* activePointer, std::ptrdiff_t offset, const std::size_t newSize)
   {
      dataManager_ = dataManager;
      activePointer_ = activePointer;
      offset_ = offset;
      size_ = newSize;
#if !defined(__CUDA_ARCH__)
      use();
#endif      
   }
   LAZY_HOST_ONLY inline void useOn(const ExecutionSpace space)
   {
      this->activePointer_ = this->dataManager_->raw(this, space, offset_);
   }
   TTT* activePointer_;
   std::ptrdiff_t offset_;
   std::size_t size_;
   lazy_array<noconst_value_type>* dataManager_;
};


template <typename TTT>
class ro_lzarray_ptr : public ro_array_ptr_interface<lazy_array_impl<ro_lzarray_ptr<TTT>,const TTT>,const TTT&,const TTT*>
{
 public:
   LAZY_HOST_DEVICE inline ro_lzarray_ptr() : ro_lzarray_ptr(NULL, NULL, 0, 0) {}
   LAZY_HOST_DEVICE inline ro_lzarray_ptr(const ro_lzarray_ptr<TTT>& other)
   : ro_lzarray_ptr(other.dataManager_, other.activePointer_, other.offset_, other.size_) {}
   LAZY_HOST_DEVICE inline ro_lzarray_ptr(lazy_array<TTT>* dataManager, const TTT* activePointer, const std::ptrdiff_t offset, const std::size_t newSize)
   { this->create(dataManager,activePointer,offset,newSize); }
   LAZY_HOST_DEVICE inline operator ro_array_ptr<TTT>() { return ro_array_ptr<TTT>(this->raw(), this->size()); }
};

template <typename TTT>
class wo_lzarray_ptr : public rw_array_ptr_interface<lazy_array_impl<wo_lzarray_ptr<TTT>,TTT>,OnlyAssignable<TTT>,AssignableIterator<TTT> >
{
 public:
   LAZY_HOST_DEVICE inline wo_lzarray_ptr() : wo_lzarray_ptr(NULL, NULL, 0, 0) {}
   LAZY_HOST_DEVICE inline wo_lzarray_ptr(const wo_lzarray_ptr<TTT>& other)
   : wo_lzarray_ptr(other.dataManager_, other.activePointer_, other.offset_, other.size_) {}
   LAZY_HOST_DEVICE inline wo_lzarray_ptr(lazy_array<TTT>* dataManager, TTT* activePointer, const std::ptrdiff_t offset, const std::size_t newSize)
   { this->create(dataManager,activePointer,offset,newSize); }
   LAZY_HOST_DEVICE inline operator wo_array_ptr<TTT>() { return wo_array_ptr<TTT>(this->raw(), this->size()); }
};

template <typename TTT>
class rw_lzarray_ptr : public rw_array_ptr_interface<lazy_array_impl<rw_lzarray_ptr<TTT>,TTT>,TTT&,TTT*>
{
 public:
   LAZY_HOST_DEVICE inline rw_lzarray_ptr() : rw_lzarray_ptr(NULL, NULL, 0, 0) {}
   LAZY_HOST_DEVICE inline rw_lzarray_ptr(const rw_lzarray_ptr<TTT>& other)
   : rw_lzarray_ptr(other.dataManager_, other.activePointer_, other.offset_, other.size_) {}
   LAZY_HOST_DEVICE inline rw_lzarray_ptr(lazy_array<TTT>* dataManager, TTT* activePointer, const std::ptrdiff_t offset, const std::size_t newSize)
   { this->create(dataManager,activePointer,offset,newSize); }
   LAZY_HOST_DEVICE inline operator ro_lzarray_ptr<TTT>()
   {
      return ro_lzarray_ptr<TTT>(this->dataManager_, this->activePointer_, this->offset_, this->size());
   }
   LAZY_HOST_DEVICE inline operator wo_lzarray_ptr<TTT>()
   {
      return wo_lzarray_ptr<TTT>(this->dataManager_, this->activePointer_, this->offset_, this->size());
   }
   LAZY_HOST_DEVICE inline operator rw_array_ptr<TTT>() { return rw_array_ptr<TTT>(this->raw(), this->size()); }

};

#endif

template <typename TTT>
class DelayedUse
{
 public:
   DelayedUse(lazy_array<TTT>* dataManager, const ExecutionSpace space)
   : dataManager_(dataManager), space_(space) {}

   inline rw_array_ptr<TTT> readwrite() { return dataManager_->readwrite(space_); }
   inline wo_array_ptr<TTT> writeonly() { return dataManager_->writeonly(space_); }
   inline ro_array_ptr<TTT> readonly() const { return dataManager_->readonly(space_); }

   inline operator rw_array_ptr<TTT>() { return readwrite(); }
   inline operator wo_array_ptr<TTT>() { return writeonly(); }
   inline operator ro_array_ptr<TTT>() { return readonly(); }
 private:
   ExecutionSpace space_;
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
   }
   ~lazy_array()
   {
      clear();
   }
   template <typename Allocator>
   lazy_array(const std::vector<TTT, Allocator>& vvv) : lazy_array()
   {
      resize(vvv.size());
      wo_array_ptr<TTT> myData = writeonly().useOn(CPU);
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
   }
   lazy_array<TTT>& operator=(lazy_array<TTT> other)
   {
      swap(*this,other);
      return *this;
   }
   
   LAZY_HOST_ONLY inline bool empty() const { return size_==0; }
   LAZY_HOST_ONLY inline std::size_t size() const { return size_; }

   LAZY_HOST_ONLY void clear() {
      for (int ispace=0; ispace<NUMSPACES; ispace++)
      {
         if (pointerRecord_[ispace] != NULL) { spaceFree(ispace,pointerRecord_[ispace]); }
         pointerRecord_[ispace] = NULL;
         isValid_[ispace] = false;
      }
      size_ = 0;
   }
   
   LAZY_HOST_ONLY void resize(const std::size_t elems)
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

   inline operator rw_mgarray_ptr<TTT>() { return readwrite(); }
   inline operator wo_mgarray_ptr<TTT>() { return writeonly(); } 
   inline operator ro_mgarray_ptr<TTT>() const { return readonly(); }

   inline rw_mgarray_ptr<TTT> readwrite() { return rw_mgarray_ptr<TTT>(this, 0, size()); }
   inline wo_mgarray_ptr<TTT> writeonly() { return wo_mgarray_ptr<TTT>(this, 0, size()); }
   inline ro_mgarray_ptr<TTT> readonly() const { return ro_mgarray_ptr<TTT>(const_cast<lazy_array<TTT>*>(this), 0, size()); }

   inline rw_array_ptr<TTT> readwrite(const ExecutionSpace space) { return readwrite().useOn(space); }
   inline wo_array_ptr<TTT> writeonly(const ExecutionSpace space) { return writeonly().useOn(space); }
   inline ro_array_ptr<TTT> readonly(const ExecutionSpace space) const { return readonly().useOn(space); }

   inline DelayedUse<TTT> useOn(const ExecutionSpace space) { return DelayedUse<TTT>(this, space); }
   ro_array_ptr<TTT> useOn(const ExecutionSpace space) const { return readonly(space); }
   
   /*
   LAZY_HOST_ONLY inline const TTT* raw(lazy_array_impl<ro_lzarray_ptr<TTT>,const TTT>*,
                              const ExecutionSpace space, const std::ptrdiff_t offset=0)
   {
      allocateFor(space);
      moveTo(space);
      
      return lookup(space,offset);      
   }
   LAZY_HOST_ONLY inline TTT* raw(lazy_array_impl<wo_lzarray_ptr<TTT>,TTT>*,
                        const ExecutionSpace space, const std::ptrdiff_t offset=0)
   {
      allocateFor(space);
      writeOn(space);
      
      return lookup(space,offset);      
   }
   LAZY_HOST_ONLY inline TTT* raw(lazy_array_impl<rw_lzarray_ptr<TTT>,TTT>*,
                        const ExecutionSpace space, const std::ptrdiff_t offset=0)
   {
      allocateFor(space);
      moveTo(space);
      writeOn(space);
      
      return lookup(space,offset);
   }
   */

   LAZY_HOST_ONLY inline const TTT* raw(managed_array_impl<ro_mgarray_ptr<TTT>,ro_array_ptr<TTT>,const TTT>*,
                                        const ExecutionSpace space, const std::ptrdiff_t offset=0)
   {
      allocateFor(space);
      moveTo(space);
      
      return lookup(space,offset);      
   }
   LAZY_HOST_ONLY inline TTT* raw(managed_array_impl<wo_mgarray_ptr<TTT>,wo_array_ptr<TTT>,TTT>*,
                                  const ExecutionSpace space, const std::ptrdiff_t offset=0)
   {
      allocateFor(space);
      writeOn(space);
      
      return lookup(space,offset);      
   }
   LAZY_HOST_ONLY inline TTT* raw(managed_array_impl<rw_mgarray_ptr<TTT>,rw_array_ptr<TTT>,TTT>*,
                                  const ExecutionSpace space, const std::ptrdiff_t offset=0)
   {
      allocateFor(space);
      moveTo(space);
      writeOn(space);
      
      return lookup(space,offset);
   }

 private:
   TTT* pointerRecord_[NUMSPACES];
   bool isValid_[NUMSPACES];
   std::size_t size_;

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
#ifdef USE_CUDA
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
         }
#endif
         else
         {
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
