#pragma once

#include "ExecutionSpace.hh"

template <typename TTT>
bool spaceMalloc(const ExecutionSpace space, TTT** dst, const std::size_t size);
template <>
bool spaceMalloc<void>(const ExecutionSpace space,void** dst, const std::size_t size);
template <typename TTT>
bool spaceMalloc(const ExecutionSpace space, TTT** dst, const std::size_t size)
{
   return spaceMalloc<void>(space,reinterpret_cast<void**>(dst),size*sizeof(TTT));
}


template <typename TTT>
void spaceMemcpy(const ExecutionSpace dstSpace, TTT* dst,
                 const ExecutionSpace srcSpace, const TTT* src,
                 const std::size_t size);
template <>
void spaceMemcpy<void>(const ExecutionSpace dstSpace, void* dst,
                       const ExecutionSpace srcSpace, const void* src,
                       const std::size_t size);
template <typename TTT>
void spaceMemcpy(const ExecutionSpace dstSpace, TTT* dst,
                 const ExecutionSpace srcSpace, const TTT* src,
                 const std::size_t size)
{
   spaceMemcpy<void>(dstSpace,reinterpret_cast<void*>(dst),
                     srcSpace,reinterpret_cast<const void*>(src),
                     size*sizeof(TTT));
}


template <typename TTT>
void spaceFree(const ExecutionSpace space,TTT* dst);
template <>
void spaceFree<void>(const ExecutionSpace space,void* dst);
template <typename TTT>
void spaceFree(const ExecutionSpace space,TTT* dst)
{
   spaceFree<void>(space,reinterpret_cast<void*>(dst));
}


