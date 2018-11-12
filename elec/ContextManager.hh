#pragma once

#include <vector>
#include "ExecutionSpace.hh"

class ContextManager
{
 public:
   inline ExecutionSpace current() const { return current_; }
   void push(const ExecutionSpace context);
   void pop();

   ContextManager();
   
 private:
   std::vector<ExecutionSpace> stack_;
   ExecutionSpace current_;
};

ContextManager* getContextManager();

class ContextRegion
{
 public:
   ContextRegion(const ExecutionSpace space);
   ~ContextRegion();
 private:
   ContextRegion(const ContextRegion& region);
   ContextRegion& operator=(const ContextRegion& region);
   ContextManager* manager_;
};
