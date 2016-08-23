#include "EndianSwap.hh"

#include <cstring>
using std::memcpy;

/** The value of the local endian key is:
 *  big endian:    825373492
 *  little endian: 875770417
 *
 * Typically intel hardware is little endian.
 * IBM hardware is usually big endian.
*/

EndianSwap::EndianSwap(unsigned endianKey)
{
   unsigned keyLocal;
   memcpy(&keyLocal, "1234", 4);
   swapNeeded_ = false;
   if (keyLocal != endianKey)
      swapNeeded_ = true;
}

void
EndianSwap::operator()(int& i)
{
   if (swapNeeded_)
      endianSwap(&i, 4);
}


void
EndianSwap::endianSwap(void* data, int size)
{ 
   char save; 
   char* b = (char*) (data) ; 
   int j = size; 
   for (int i=0; i<size/2; i++) 
   {
      --j; 
      save = b[i] ;  
      b[i] = b[j]; 
      b[j] = save; 
   }
} 

