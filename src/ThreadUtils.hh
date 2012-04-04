#ifndef THREAD_UTILS_HH
#define THREAD_UTILS_HH

#include <vector>
class ThreadTeam;

void mkOffsets(std::vector<int>& offset, int nItems, const ThreadTeam& threadInfo);

#endif
