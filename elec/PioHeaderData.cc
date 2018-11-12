#include "PioHeaderData.hh"
#include <cassert>
#include <sstream>
#include <cstring>

#include "Version.hh"
#include "utilities.h"
#include "pio.h"

using namespace std;

void PioHeaderData::addItem(const string& keyword, const string& value)
{
   otherItems_[keyword] = value;
}

void PioHeaderData::addItem(const string& keyword, int value)
{
   stringstream buf;
   buf << value;
   addItem(keyword, buf.str());
}
      
void PioHeaderData::addItem(const string& keyword, unsigned value)
{
   stringstream buf;
   buf << value;
   addItem(keyword, buf.str());
}

void PioHeaderData::addItem(const string& keyword, double value)
{
   stringstream buf;
   buf << value;
   addItem(keyword, buf.str());
}


void PioHeaderData::writeHeader(PFILE* file, int loop, double time)
{
   string dataType;
   switch (dataType_)
   {
     case PioHeaderData::ASCII:
      dataType = "FIXRECORDASCII";
      break;
     case PioHeaderData::BINARY:
      dataType = "FIXRECORDBINARY";
      break;
     default:
      assert(false);
   }
   int endianKey;
   memcpy(&endianKey, "1234", 4);
   int nFiles = file->ngroup;
   
   const Version& vv = Version::getInstance();
   
   Pprintf(file, "%s %s {\n", objectName_.c_str(), className_.c_str());
   Pprintf(file, "   create_time = %s;\n",  timestamp_string());
   Pprintf(file, "   user = %s;",           vv.user().c_str());
   Pprintf(file, "   host = %s;\n",         vv.host().c_str());
   Pprintf(file, "   exe_version = %s;",    vv.version().c_str());
   Pprintf(file, "   srcpath = %s;\n",      vv.srcPath().c_str());
   Pprintf(file, "   compile_time = %s;\n", vv.compileTime().c_str());
   Pprintf(file, "   datatype = %s;\n", dataType.c_str());
   Pprintf(file, "   nfiles = %d;\n", nFiles);
   Pprintf(file, "   nrecord = %d;\n", nRecords_);
   Pprintf(file, "   lrec = %d;\n", lRec_);
   Pprintf(file, "   endian_key = %d;\n", endianKey);
   Pprintf(file, "   time = %f;\n", time);
   Pprintf(file, "   loop = %d;\n", loop);
   Pprintf(file, "   nfields = %d;\n", nFields_);
   Pprintf(file, "   field_names = %s;\n", fieldNames_.c_str());
   Pprintf(file, "   field_types = %s;\n", fieldTypes_.c_str());
   Pprintf(file, "   field_units = %s;\n", fieldUnits_.c_str());
   for (map<string, string>::iterator iter = otherItems_.begin();
        iter != otherItems_.end() ; ++iter)
      Pprintf(file, "   %s = %s;\n", iter->first.c_str(), iter->second.c_str());
   Pprintf(file, "}\n\n");
}
