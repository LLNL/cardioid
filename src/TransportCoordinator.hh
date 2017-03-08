#ifndef TRANSPORTCOORDINATOR_HH
#define TRANSPORTCOORDINATOR_HH

#include <vector>
#include <cassert>

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
   const TTT* dataAlias = &data[0];
#pragma omp target update to(dataAlias[0:data.size()])
   for (const TTT& item : data) {
      copyToDevice(item);
   }
}
template <typename TTT,typename AAA>
inline void copyToHost(const std::vector<TTT, AAA>& data)
{
   const TTT* dataAlias = &data[0];
#pragma omp target update from(dataAlias[0:data.size()])
   for (const TTT& item : data) {
      copyToHost(item);
   }
}
template <typename TTT, typename AAA>
inline void allocOnDevice(const std::vector<TTT, AAA>& data)
{
   const TTT* dataAlias = &data[0];
#if _OPENMP >= 201511
#pragma omp target enter data map(alloc: dataAlias[0:data.size()])
#endif
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
   const TTT* dataAlias = &data[0];
#if _OPENMP >= 201511
#pragma omp target exit data map(release: dataAlias[0:data.size()])
#endif   
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
         copyToHost(data_);
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
