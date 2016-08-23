#ifndef ENDIAN_SWAP_HH
#define ENDIAN_SWAP_HH

class EndianSwap
{
 public:

   EndianSwap(unsigned endianKey);

   void operator()(int& i);
      
   
 private:

   void endianSwap(void* data, int size);
   
   bool swapNeeded_;

};

#endif
