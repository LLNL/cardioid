
#include "Ledger.hh"
#include <map>
#include <cassert>
#include <cuda.h>
#include <cuda_runtime_api.h>

using namespace std;

static map<const void*, pair<size_t, void*> > ledgerMap_;

void ledger_init()
{
   ledgerMap_.clear();
}

void ledger_alloc(const void* host, std::size_t size)
{
   assert(ledgerMap_.find(host) == ledgerMap_.end());
   void* devicePtr;
   cudaMalloc(&devicePtr, size);
   ledgerMap_[host] = make_pair(size, devicePtr);
}
void ledger_free(const void* host)
{
   ledgerMap_.erase(host);
}
void ledger_copyToDevice(const void* host)
{
   map<const void*, pair<size_t, void* > >::iterator iter = ledgerMap_.find(host);
   assert(iter != ledgerMap_.end());
   cudaMemcpy(iter->second.second, host, iter->second.first, cudaMemcpyHostToDevice);
}
void ledger_copyToHost(void* host)
{
   map<const void*, pair<size_t, void* > >::iterator iter = ledgerMap_.find(host);
   assert(iter != ledgerMap_.end());
   cudaMemcpy(host, iter->second.second, iter->second.first, cudaMemcpyDeviceToHost);   
}
void* ledger_lookup(const void* host)
{
   map<const void*, pair<size_t, void* > >::iterator iter = ledgerMap_.find(host);
   assert(iter != ledgerMap_.end());
   return iter->second.second;
}
void ledger_deviceZero(const void* host)
{
   map<const void*, pair<size_t, void* > >::iterator iter = ledgerMap_.find(host);
   assert(iter != ledgerMap_.end());
   cudaMemset(iter->second.second, 0, iter->second.first);
}
