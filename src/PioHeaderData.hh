#ifndef PIO_HEADER_DATA_HH
#define PIO_HEADER_DATA_HH

#include <string>
#include <map>
#include "Long64.hh"

struct pfile_st;


struct PioHeaderData
{
   enum DataType {ASCII, BINARY};


   void addItem(const std::string& keyword, const std::string& value);
   void addItem(const std::string& keyword, int value);
   void addItem(const std::string& keyword, unsigned value);
   void addItem(const std::string& keyword, double value);
   
   void writeHeader(pfile_st* file, int loop, double time);


   
   std::string objectName_;
   std::string className_;
   DataType dataType_;
   Long64 nRecords_;
   unsigned lRec_;
   unsigned nFields_;
   std::string fieldNames_;
   std::string fieldTypes_;
   std::string fieldUnits_;

   std::map<std::string, std::string> otherItems_;
};

#endif
