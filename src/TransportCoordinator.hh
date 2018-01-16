#ifndef TRANSPORTCOORDINATOR_HH
#define TRANSPORTCOORDINATOR_HH

#include <vector>
#include <cassert>
#include "Ledger.hh"

template <typename TTT>
inline void copyToHost(const TTT& data) {};
template <typename TTT>
inline void copyToDevice(const TTT& data) {};
template <typename TTT>
inline void allocOnDevice(const TTT& data) {};
template <typename TTT>
inline void freeOnDevice(const TTT& data) {};

template <typename TTT,typename AAA>
inline void copyToDevice(const std::vector<TTT, AAA>& data)
{
   ledger_copyToDevice(&data[0]);
   for (const TTT& item : data) {
      copyToDevice(item);
   }
}
template <typename TTT,typename AAA>
inline void copyToHost(std::vector<TTT, AAA>& data)
{
   ledger_copyToHost(&data[0]);
   for (const TTT& item : data) {
      copyToHost(item);
   }
}
template <typename TTT, typename AAA>
inline void allocOnDevice(const std::vector<TTT, AAA>& data)
{
   ledger_alloc(&data[0], data.size()*sizeof(TTT));
   for (const TTT& item : data) {
      allocOnDevice(item);
   }
}
template <typename TTT, typename AAA>
inline void freeOnDevice(const std::vector<TTT,AAA>& data)
{
   for (const TTT& item : data) {
      freeOnDevice(item);
   }
   ledger_free(&data[0]);
}

template <typename TTT>
class TransportCoordinator {
 public:
   TransportCoordinator()
   {
      allocated_ = false;
      isHostValid_ = false;
      isDeviceValid_ = false;
   }
   TransportCoordinator(TTT&& initData)
   {
      setup(std::move(initData));
   }
   ~TransportCoordinator()
   {
      teardown();
   }
   inline void setup(TTT&& initData)
   {
      data_ = initData;
      allocOnDevice(data_);
      allocated_ = true;
      isHostValid_ = true;
      isDeviceValid_ = false;
   }
   inline void teardown()
   {
      if (allocated_)
      {
         freeOnDevice(data_);
         allocated_ = false;
         isHostValid_ = false;
         isDeviceValid_ = false;
      }
   }      
   inline TTT& modifyOnHost()
   {
      readOnHost();
      isDeviceValid_ = false;
      return data_;
   }
   inline TTT& modifyOnDevice()
   {
      readOnDevice();
      isHostValid_ = false;
      return data_;
   }
   inline const TTT& readOnHost() const
   {
      if (isHostValid_==false)
      {
         assert(isDeviceValid_ == true);
         copyToHost(const_cast<TTT&>(data_));
         const_cast<bool&>(isHostValid_) = true;
      }
      return data_;
   }
   inline const TTT& readOnDevice() const
   {
      if (isDeviceValid_==false)
      {
         assert(isHostValid_ == true);
         copyToDevice(data_);
         const_cast<bool&>(isDeviceValid_) = true;
      }
      return data_;
   }
 private:
   TTT data_;
   bool allocated_;
   bool isHostValid_;
   bool isDeviceValid_;
};

#endif //TRANSPORTCOORDINATOR_HH
