#ifndef TRANSPORTCOORDINATOR_HH
#define TRANSPORTCOORDINATOR_HH

#include <vector>
#include <cassert>
#include "Ledger.hh"
#include "ArrayView.hh"
#include "PinnedAllocator.hh"


template <typename TTT>
inline void copyToHost(const TTT* data) {};
template <typename TTT>
inline void copyToDevice(const TTT* data) {};
template <typename TTT>
inline void allocOnDevice(const TTT* data) {};
template <typename TTT>
inline void freeOnDevice(const TTT* data) {};

template <typename TTT>
inline void copyToDevice(const ConstArrayView<TTT> data)
{
   ledger_copyToDevice(&data[0]);
   for (const TTT& item : data) {
      copyToDevice(&item);
   }
}
template <typename TTT>
inline void copyToHost(ArrayView<TTT> data)
{
   ledger_copyToHost(&data[0]);
   for (const TTT& item : data) {
      copyToHost(&item);
   }
}
template <typename TTT>
inline void allocOnDevice(ConstArrayView<TTT> data)
{
   ledger_alloc(&data[0], data.size()*sizeof(TTT));
   for (const TTT& item : data) {
      allocOnDevice(&item);
   }
}
template <typename TTT>
inline void freeOnDevice(ConstArrayView<TTT> data)
{
   for (const TTT& item : data) {
      freeOnDevice(&item);
   }
   ledger_free(&data[0]);
}


template <typename View>
class OnDevice;

template <typename TTT>
class OnDevice<TTT*>
{
 public:
   TTT* weak_;
   explicit OnDevice(TTT* view) : weak_(view) {};
   
   operator TTT*() {
      return weak_;
   }
   operator const TTT*() const {
      return weak_;
   }
   TTT& operator*() { return *(weak_); }
   const TTT& operator*() const { return *(weak_); }
};

template <typename TTT>
class OnDevice<ArrayView<TTT>> : public ArrayView<TTT>
{
  public:
   explicit OnDevice(ArrayView<TTT> view) : ArrayView<TTT>(ledger_lookup(&view[0]),view.size()) {}
};

template <typename TTT>
class OnDevice<ConstArrayView<TTT>> : public ConstArrayView<TTT>
{
  public:
   explicit OnDevice(ConstArrayView<TTT> view) : ConstArrayView<TTT>(ledger_lookup(&view[0]),view.size()) {}
};


//////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename TTT>
class ViewForType {
 public:
   typedef TTT* result;
   static inline result create(TTT& ttt) { return &ttt; }
};

template <typename TTT, typename AAA>
class ViewForType<std::vector<TTT,AAA>>
{
  public:
   typedef ArrayView<TTT> result;
   static inline result create(std::vector<TTT,AAA>& ttt) { return result(ttt); }
};

template <typename TTT>
class ConstViewFromView;

template <typename TTT>
class ConstViewFromView<TTT*> {
 public:
   typedef const TTT* result;
   static inline TTT* remove_const(const TTT* arg) { return const_cast<TTT*>(arg); }
};

template <typename TTT>
class ConstViewFromView<ArrayView<TTT>> {
  public:
   typedef ConstArrayView<TTT> result;
   static inline ArrayView<TTT> remove_const(ConstArrayView<TTT> arg) { return ArrayView<TTT>(const_cast<TTT*>(&arg[0]),arg.size()); }
};

////////////////////////////////////////////////////////////////////////////////////////////////////


template <typename View>
class ViewManager
{
 public:
   typedef typename ConstViewFromView<View>::result ConstView;
   
   ViewManager()
   {
      allocated_ = false;
      isHostValid_ = false;
      isDeviceValid_ = false;
   }
   ViewManager(View initData)
   {
      setup(initData);
   }
   ~ViewManager()
   {
      teardown();
   }
   inline void setup(View initData)
   {
      data_ = initData;
      allocOnDevice(ConstView(data_));
      allocated_ = true;
      isHostValid_ = true;
      isDeviceValid_ = false;
   }
   inline void teardown()
   {
      if (allocated_)
      {
         freeOnDevice(ConstView(data_));
         allocated_ = false;
         isHostValid_ = false;
         isDeviceValid_ = false;
      }
   }
   inline View modifyOnHost()
   {
      readOnHost();
      isDeviceValid_ = false;
      return View(data_);
   }
   inline OnDevice<View> modifyOnDevice()
   {
      readOnDevice();
      isHostValid_ = false;
      return OnDevice<View>(data_);
   }
   inline ConstView readOnHost() const
   {
      if (isHostValid_==false)
      {
         assert(isDeviceValid_ == true);
         copyToHost(ConstViewFromView<View>::remove_const(data_));
         const_cast<bool&>(isHostValid_) = true;
      }
      return ConstView(data_);
   }
   inline OnDevice<ConstView> readOnDevice() const
   {
      if (isDeviceValid_==false)
      {
         assert(isHostValid_ == true);
         copyToDevice(ConstView(data_));
         const_cast<bool&>(isDeviceValid_) = true;
      }
      return OnDevice<ConstView>(data_);
   }
   inline View raw() {
      return data_;
   }
   inline ConstView raw() const {
      return ConstView(data_);
   }  
 private:
   View data_;
   bool allocated_;
   bool isHostValid_;
   bool isDeviceValid_;   
};

template <typename View>
class Managed
{
 public:
   typedef typename ConstViewFromView<View>::result ConstView;

   explicit Managed(ViewManager<View>* p_manager) : p_manager_(p_manager) {}
   explicit Managed(const ViewManager<View>* p_manager) : p_manager_(const_cast<ViewManager<View>*>(p_manager)) {}
   
   inline View modifyOnHost() { return p_manager_->modifyOnHost(); }
   inline OnDevice<View> modifyOnDevice() { return p_manager_->modifyOnDevice(); }
   inline ConstView readOnHost() const { return p_manager_->readOnHost(); }
   inline OnDevice<ConstView> readOnDevice() const { return p_manager_->readOnDevice(); }

   inline operator View() { return modifyOnHost(); }
   inline operator OnDevice<View>() { return modifyOnDevice(); }
   inline operator ConstView() const { return readOnHost(); }
   inline operator OnDevice<ConstView>() const { return readOnDevice(); }

   //inline View host() { return modifyOnHost(); }
   //inline ConstView host() const { return readOnHost(); }

   //inline View device() { return modifyOnDevice(); }
   //inline ConstView device() const { return readOnDevice(); }

   inline View raw() { return p_manager_->raw(); }
   inline ConstView raw() const { return p_manager_->raw(); }
   
 private:
   ViewManager<View> * p_manager_;
};

template <typename TTT>
class TransportCoordinator {
 public:
   typedef typename ViewForType<TTT>::result View;
   typedef typename ConstViewFromView<View>::result ConstView;

   TransportCoordinator() {}
   TransportCoordinator(TTT&& initData)
   {
      setup(std::move(initData));
   }
   inline void setup(TTT&& initData)
   {
      data_ = initData;
      manager_.setup(ViewForType<TTT>::create(data_));
      
   }
   inline void teardown()
   {
      manager_.teardown();
   }
   
   inline Managed<View> view() { return Managed<View>(&manager_); }
   inline const Managed<View> view() const { return Managed<View>(&manager_); }

   inline operator Managed<View>() { return view(); }
   inline operator const Managed<View>() const { return view(); }

   inline operator View() { return view(); }
   inline operator OnDevice<View>() { return view(); }
   inline operator ConstView() const { return view(); }
   inline operator OnDevice<ConstView>() const { return view(); }

   
   inline View modifyOnHost() { return view().modifyOnHost(); }
   inline OnDevice<View> modifyOnDevice() { return view().modifyOnDevice(); }
   inline ConstView readOnHost() const { return view().readOnHost(); }
   inline OnDevice<ConstView> readOnDevice() const { return view().readOnDevice(); }

   //inline View host() { return view().host(); }
   //inline ConstView host() const { return view().host(); }

   //inline View device() { return view().device(); }
   //inline ConstView device() const { return view().device(); }

   inline View raw() { return manager_.raw(); }
   inline ConstView raw() const { return manager_.raw(); }
 private:
   ViewManager<View> manager_;
   TTT data_;
};




#endif //TRANSPORTCOORDINATOR_HH
