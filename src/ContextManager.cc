
#include "ContextManager.hh"

void ContextManager::push(const ExecutionSpace context)
{
   stack_.push_back(context);
   current_ = context;
}
void ContextManager::pop()
{
   stack_.pop_back();
   current_ = stack_.back();
}

ContextManager::ContextManager()
{
   push(NONE);
}

static ContextManager g_contextManager;

ContextManager* getContextManager() { return &g_contextManager; }

ContextRegion::ContextRegion(const ExecutionSpace space)
{
   manager_ = getContextManager();
   manager_->push(space);
}
ContextRegion::~ContextRegion()
{
   manager_->pop();
}
