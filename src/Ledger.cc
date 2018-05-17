
#include "Ledger.hh"
#include <map>
#include <cassert>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <iostream>

#define CUDA_VERIFY(x) do { cudaError_t error = x; if (error != cudaSuccess) { cout << error << endl; assert(error == cudaSuccess && #x ); } } while(0)

using namespace std;

static map<const void*, pair<size_t, void*> > ledgerMap_;
static int nullCount = 0;

void ledger_init()
{
   ledgerMap_.clear();
}

void ledger_alloc(const void* host, std::size_t size)
{
   if (host ==NULL) {
      assert(size==0);
      nullCount++;
      return;
   }
   assert(ledgerMap_.find(host) == ledgerMap_.end());
   void* devicePtr=NULL;
   CUDA_VERIFY(cudaMalloc(&devicePtr, size));
   ledgerMap_[host] = make_pair(size, devicePtr);
}
void ledger_free(const void* host)
{
   if (host==NULL) {
      assert(nullCount > 0);
      nullCount--;
      return;
   }
   map<const void*, pair<size_t, void* > >::iterator iter = ledgerMap_.find(host);
   assert(iter != ledgerMap_.end());
   void* devicePtr = iter->second.second;
   ledgerMap_.erase(iter);
   CUDA_VERIFY(cudaFree(devicePtr));
}
void ledger_copyToDevice(const void* host)
{
   if (nullCount && host==NULL) { return; }
   map<const void*, pair<size_t, void* > >::iterator iter = ledgerMap_.find(host);
   assert(iter != ledgerMap_.end());
   CUDA_VERIFY(cudaMemcpy(iter->second.second, host, iter->second.first, cudaMemcpyHostToDevice));
}
void ledger_copyToHost(void* host)
{
   if (nullCount && host==NULL) { return; }
   map<const void*, pair<size_t, void* > >::iterator iter = ledgerMap_.find(host);
   assert(iter != ledgerMap_.end());
   CUDA_VERIFY(cudaMemcpy(host, iter->second.second, iter->second.first, cudaMemcpyDeviceToHost));   
}
template<>
void* ledger_lookup<void>(const void* host)
{
   if (nullCount && host==NULL) { return NULL; }
   assert(ledgerMap_.size() >= 1);
   map<const void*, pair<size_t, void* > >::iterator iter = ledgerMap_.upper_bound(host);
   assert(iter != ledgerMap_.begin());
   --iter;
   size_t diff = (char*)(host) - (char*)(iter->first);
   assert(diff <= iter->second.first);
   return (void*)((char*)iter->second.second + diff);
}
void ledger_deviceZero(const void* host)
{
   if (nullCount && host==NULL) { return; }
   map<const void*, pair<size_t, void* > >::iterator iter = ledgerMap_.find(host);
   assert(iter != ledgerMap_.end());
   CUDA_VERIFY(cudaMemset(iter->second.second, 0, iter->second.first));
}
