#include "ThreadUtils.hh"
#include <cassert>
#include "ThreadServer.hh"

using namespace std;

void mkOffsets(vector<int>& offset, int nItems, const ThreadTeam& threadInfo)
{
   int chunkSize = nItems / threadInfo.nThreads();
   int leftOver  = nItems % threadInfo.nThreads();
   offset.resize(threadInfo.nThreads()+1);
   offset[0] = 0;
   for (int ii=0; ii<threadInfo.nThreads(); ++ii)
   {
      offset[ii+1] = offset[ii] + chunkSize;
      if (ii < leftOver)
         ++offset[ii+1];
   }
   assert(offset[threadInfo.nThreads()] == nItems);
}
