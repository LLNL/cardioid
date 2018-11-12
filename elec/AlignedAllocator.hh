#ifndef ALIGNED_ALLOCATOR_HH
#define ALIGNED_ALLOCATOR_HH

#include <new>
#include <limits>
#include <cstdlib>
#include <cstddef>

template <class T> class AlignedAllocator;
 
// specialize for void:
template <> class AlignedAllocator<void>
{
 public:
   typedef void*       pointer;
   typedef const void* const_pointer;
   // reference to void members are impossible.
   typedef void value_type;
   template <class U> struct rebind { typedef AlignedAllocator<U>
                                      other; };
};

template <class T> class AlignedAllocator
{
 public:
   typedef size_t    size_type;
   typedef ptrdiff_t difference_type;
   typedef T*        pointer;
   typedef const T*  const_pointer;
   typedef T&        reference;
   typedef const T&  const_reference;
   typedef T         value_type;
   template <class U> struct rebind { typedef AlignedAllocator<U>
                                      other; };
   
   AlignedAllocator() throw() {}
   AlignedAllocator(const AlignedAllocator&) throw() {}
   template <class U> AlignedAllocator(const AlignedAllocator<U>&) throw() {}
   ~AlignedAllocator() throw() {}
   
   pointer address(reference x) const {return &x;}
   const_pointer address(const_reference x) const {return &x;}
   
   pointer allocate(size_type size,
                    AlignedAllocator<void>::const_pointer hint = 0)
   {
      pointer ptr;
      if (posix_memalign((void**)&ptr, 32, size*sizeof(T)) == 0)
         return ptr;
      throw std::bad_alloc();
   }
   
   void deallocate(pointer p, size_type n) {free(p);}
   size_type max_size() const throw()
   {
      return std::numeric_limits<size_type>::max() / sizeof(T);
   }
   
   void construct(pointer p, const T& val){ new((void*)p) T(val);}
   void destroy(pointer p) {p->~T();}
};


template <class T1, class T2>
bool operator==(const AlignedAllocator<T1>&, const AlignedAllocator<T2>&)
   throw() {return true;}

template <class T1, class T2>
bool operator!=(const AlignedAllocator<T1>&, const AlignedAllocator<T2>&)
   throw() {return false;}


#endif
