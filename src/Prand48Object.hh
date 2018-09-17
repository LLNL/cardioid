#ifndef PRAND48_OBJECT_HH
#define PRAND48_OBJECT_HH

#include <cstdlib>
#include <stdint.h> // uint64_t

#include <iostream>

/** Random number generator that combines a particle gid, a seed, and a
 *  call site id to generate a random number.  Since the underlying
 *  random engine is erand48, only the 48 low order bits in the three
 *  parameters are used.
 *
 *  The rationale for generating random numbers this way is that for a
 *  fixed user specified seed, we want reproducible results, regardless
 *  of the number of tasks the code is running on.  Since we generate a
 *  new seed for erand based on the user seed and the particle gid, we
 *  get the same number (for a given seed) regardless of the calling
 *  history of erand.  The callSite parameter allows for different call
 *  sites in the code to get different pseudo-random streams even for
 *  the same particle and user seed.
 *
 *  Issues: Even when we discard a fairly large number of initial random
 *  numbers there tends to be correlations between the random numbers
 *  for nearly equal gids.  If the gids are somehow correlated and good
 *  random numbers are important for your application, you will likely
 *  want to chose a different method, or increase the discard size.
 */
class Prand48Object
{
 public:

   Prand48Object(uint64_t gid, uint64_t seed, uint64_t callSite)
   {
      unsigned nDiscard = gid%1000;
      seed ^= callSite;
      gid ^= seed;
      
      seedVec_[0] = seedVec_[1] = seedVec_[2] = 0;
      unsigned short mask[16] = {0x0001, 0x0002, 0x0004, 0x0008,
                                 0x0010, 0x0020, 0x0040, 0x0080,
                                 0x0100, 0x0200, 0x0400, 0x0800,
                                 0x1000, 0x2000, 0x4000, 0x8000};
      for (unsigned ii=0; ii<16; ++ii)
      {
         seedVec_[1] += (gid & mask[ii]); gid = gid >> 1;
         seedVec_[0] += (gid & mask[ii]); gid = gid >> 1;
         seedVec_[2] += (gid & mask[ii]); 
      }

      for (unsigned ii=0; ii<nDiscard; ++ii)
          erand48(seedVec_);
   }

   double operator()(void)
   {
      return erand48(seedVec_);
   }

   int operator()(int maxVal)
   {
      return int(erand48(seedVec_) * maxVal);
   }
   
 private:
   unsigned short seedVec_[3];
};
   

#endif
