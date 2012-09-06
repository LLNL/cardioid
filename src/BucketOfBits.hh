#ifndef BUCKET_OF_BITS_HH
#define BUCKET_OF_BITS_HH

#include <string>
#include <vector>
#include <stdint.h>

class BucketOfBits
{
 public:

   enum DataType{floatType, intType, stringType};

   class Record
   {
    public:
      Record(const std::vector<DataType>& fieldTypes,
             const std::string& rawData);
   
      void getValue(unsigned fieldIndex, double& value) const;
      void getValue(unsigned fieldIndex, int& value) const;
      void getValue(unsigned fieldIndex, uint64_t& value) const;
      void getValue(unsigned fieldIndex, std::string& value) const;
      
    private:
      std::vector<DataType>    fieldTypes_;
      std::vector<std::string> rawFields_;
   };

   BucketOfBits(const std::vector<std::string>& fieldNames,
                const std::vector<std::string>& fieldTypes,
                const std::vector<std::string>& fieldUnits);
   
   unsigned nRecords() const;
   unsigned nFields() const;
   unsigned getIndex(const std::string& fieldName) const;
   const std::string& fieldName(unsigned index) const;
   const std::string& units(unsigned index) const;
   DataType dataType(unsigned index) const;
   Record getRecord(unsigned index) const;

   void addRecord(const std::string& rec);
   
 private:
   std::vector<DataType>    fieldTypes_;
   std::vector<std::string> fieldNames_;
   std::vector<std::string> fieldUnits_;
   std::vector<std::string> records_;
};

#endif
