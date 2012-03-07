#ifndef IBM_INTRINSICS_HH
#define IBM_INTRINSICS_HH

/** This is just a wrapper around portable versions of the IBM intrinsic
 *  functions for BGQ.  It handles the logic of deciding how to best
 *  implement the intrisics for non-Q platforms.  That means you can
 *  just #include his file and be confident that the intrinsics will
 *  work the best way possible.
 */

// use the portable replacement funtions by default.
#define USE_PORTABLE_SIMD


// BGQ has the real SIMD functions.  Don't use replacements
#ifdef BGQ
#undef USE_PORTABLE_SIMD
#endif

#ifdef USE_PORTABLE_SIMD
#include "portableSIMD.hh"
#endif


#endif // IBM_INTRINSICS_HH
