#include "HaloExchange.hh"

using namespace std;

template <class T>
HaloExchange<T>::HaloExchange(const vector<int>& sendMap, const CommTable& comm)
: width_(sizeof(T)),
  commTable_(comm),
  sendMap_(sendMap_)
{
}

template <class T>
void HaloExchange<T>::execute(vector<T>& data, unsigned nLocal)
{
   // Create and fill send buffer
   vector<T> sendBuf;
   sendBuf.reserve(sendMap_.size());
   for (unsigned ii=0; ii<sendMap_.size(); ++ii)
   {
      assert(sendMap_[ii] < data.size());
      sendBuf.push_back(data[sendMap_[ii]]);
   }
   
   // make room for received data
   data.resize(nLocal+commTable_.nRemote());

   void* sendPtr = (void*) &sendBuf[0];
   void* recvPtr = (void*) &data[nLocal];

   commTable_.execute(sendPtr, recvPtr, width_);
}
