
#include <algorithm>
#include <cstddef>
#include <limits>

template <class T, int alignment> class AlignedAllocator;

//specialize for void.
template <int alignment> class AlignedAllocator<void,alignment>
{
 public:
   typedef void* pointer;
   typedef const void* const_pointer;
   typedef void value_type;
   template <class U> struct rebind { typedef AlignedAllocator<U,alignment> other; };
};

template <class T, int alignment> class AlignedAllocator
{
 public:
   typedef std::size_t    size_type;
   typedef std::ptrdiff_t difference_type;
   typedef T*        pointer;
   typedef const T*  const_pointer;
   typedef T&        reference;
   typedef const T&  const_reference;
   typedef T         value_type;
   template <class U> struct rebind { typedef AlignedAllocator<U,alignment>
                                      other; };
   
   AlignedAllocator() throw() {}
   AlignedAllocator(const AlignedAllocator&) throw() {}
   template <class U> AlignedAllocator(const AlignedAllocator<U,alignment>&) throw() {}
   ~AlignedAllocator() throw() {}
   
   pointer address(reference x) const {return &x;}
   const_pointer address(const_reference x) const {return &x;}
   
   pointer allocate(size_type size,
                    typename AlignedAllocator<void,alignment>::const_pointer hint = 0)
   {
      pointer ptr;
      if (posix_memalign((void**)&ptr, alignment, size*sizeof(T)) == 0)
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
