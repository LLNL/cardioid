#ifndef HALO_EXCHANGE_HH
#define HALO_EXCHANGE_HH

#include <vector>
#include "CommTable.hh"
#include "PerformanceTimers.hh"

template <class T>
class HaloExchange
{
 public:
   HaloExchange(const std::vector<int>& sendMap,
                const CommTable& comm);

   void execute(std::vector<T>& data, unsigned nLocal);
   
 private:
   unsigned width_;
   CommTable commTable_;
   std::vector<int> sendMap_;
};

template <class T>
HaloExchange<T>::HaloExchange(const std::vector<int>& sendMap,
                              const CommTable& comm)
: width_(sizeof(T)),
  commTable_(comm),
  sendMap_(sendMap)
{
}

template <class T>
void HaloExchange<T>::execute(std::vector<T>& data, unsigned nLocal)
{
   // Create and fill send buffer
   static TimerHandle bufFillHandle = profileGetHandle("Buffer Fill");
   profileStart(bufFillHandle);
   std::vector<T> sendBuf;
   sendBuf.reserve(sendMap_.size());
   for (unsigned ii=0; ii<sendMap_.size(); ++ii)
   {
      assert(sendMap_[ii] < data.size());
      sendBuf.push_back(data[sendMap_[ii]]);
   }
   // make room for received data
   data.resize(nLocal+commTable_.nRemote());
   profileStop(bufFillHandle);

   void* sendPtr = (void*) &sendBuf[0];
   void* recvPtr = (void*) &data[nLocal];

   commTable_.execute(sendPtr, recvPtr, width_);
}
#endif
