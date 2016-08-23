#ifndef HEART_IO_HH
#define HEART_IO_HH

#include <set>
#include <string>

struct HeaderInfo
{
   std::string dataType;
   int nFiles;
   unsigned long long nRec;
   int lRec;
   int nFields;
   std::string fieldNames;
   std::string fieldTypes;
};

void writeHeader(const std::string& filename,
		 const HeaderInfo& hdr,
		 int nx, int ny, int nz,
		 const std::set<int>& cellSet);


std::string headSpace(int size);


#endif
