#ifndef PINNED_ALLOCATOR_HH
#define PINNED_ALLOCATOR_HH

#include <new>
#include <limits>
#include <cstdlib>
#include <cstddef>
#include <vector>
#include <cuda.h>
#include <cuda_runtime_api.h>

template <class T> class PinnedAllocator;
 
// specialize for void:
template <> class PinnedAllocator<void>
{
 public:
   typedef void*       pointer;
   typedef const void* const_pointer;
   // reference to void members are impossible.
   typedef void value_type;
   template <class U> struct rebind { typedef PinnedAllocator<U>
                                      other; };
};

template <class T> class PinnedAllocator
{
 public:
   typedef size_t    size_type;
   typedef ptrdiff_t difference_type;
   typedef T*        pointer;
   typedef const T*  const_pointer;
   typedef T&        reference;
   typedef const T&  const_reference;
   typedef T         value_type;
   template <class U> struct rebind { typedef PinnedAllocator<U>
                                      other; };
   
   PinnedAllocator() throw() {}
   PinnedAllocator(const PinnedAllocator&) throw() {}
   template <class U> PinnedAllocator(const PinnedAllocator<U>&) throw() {}
   ~PinnedAllocator() throw() {}
   
   pointer address(reference x) const {return &x;}
   const_pointer address(const_reference x) const {return &x;}
   
   pointer allocate(size_type size,
                    PinnedAllocator<void>::const_pointer hint = 0)
   {
      pointer ptr;
      cudaError_t status = cudaMallocHost((void**)&ptr, size*sizeof(T));
      if (status == cudaSuccess)
      {
         return ptr;
      }
      throw std::bad_alloc();
   }
   
   void deallocate(pointer p, size_type n) {cudaFreeHost(p);}
   size_type max_size() const throw()
   {
      return std::numeric_limits<size_type>::max() / sizeof(T);
   }
   
   void construct(pointer p, const T& val){ new((void*)p) T(val);}
   void destroy(pointer p) {p->~T();}
};


template <class T1, class T2>
bool operator==(const PinnedAllocator<T1>&, const PinnedAllocator<T2>&)
   throw() {return true;}

template <class T1, class T2>
bool operator!=(const PinnedAllocator<T1>&, const PinnedAllocator<T2>&)
   throw() {return false;}

template <typename TTT>
using PinnedVector = std::vector<TTT, PinnedAllocator<TTT>>;

#endif
