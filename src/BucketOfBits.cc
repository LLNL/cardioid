#include "BucketOfBits.hh"

#include <cassert>
#include <cstdlib>

using namespace std;

namespace
{
   /** Returns the next word in the string buf starting from index
    * pos.  A word is defined as all characters between blocks of
    * whitespace.  The word will not include whitespace.  On return,
    * pos will contain the index one past the end of where the word
    * was found in buf (or string::npos) */
   string nextWord(const string& buf, string::size_type& pos);
}

BucketOfBits::BucketOfBits(const vector<string>& fieldNames,
                           const vector<string>& fieldTypes)
:fieldNames_(fieldNames)
{
   assert(fieldNames.size() == fieldTypes.size());

   for (unsigned ii=0; ii<fieldTypes.size(); ++ii)
   {
      if (fieldTypes[ii] == "s")
         fieldTypes_.push_back(stringType);
      else if (fieldTypes[ii] == "f")
         fieldTypes_.push_back(floatType);
      else if (fieldTypes[ii] == "u")
         fieldTypes_.push_back(intType);
      else
         assert(false);
   }
}


BucketOfBits::Record::Record(const vector<DataType>& fieldTypes,
                             const string& rawData)
: fieldTypes_(fieldTypes)
{
   string::size_type wordBegin = 0;
   rawFields_.reserve(fieldTypes_.size());
   for (unsigned ii=0; ii<fieldTypes_.size(); ++ii)
      rawFields_.push_back(nextWord(rawData, wordBegin));      
}

void BucketOfBits::Record::getValue(unsigned fieldIndex, double& value) const
{
   assert(fieldTypes_[fieldIndex] == floatType);
   value = strtod(rawFields_[fieldIndex].c_str(), NULL);
}

void BucketOfBits::Record::getValue(unsigned fieldIndex, int& value) const
{
   assert(fieldTypes_[fieldIndex] == intType);
   value = strtol(rawFields_[fieldIndex].c_str(), NULL, 10);
}

void BucketOfBits::Record::getValue(unsigned fieldIndex, unsigned long long& value) const
{
   assert(fieldTypes_[fieldIndex] == intType);
   value = strtol(rawFields_[fieldIndex].c_str(), NULL, 10);
}

void BucketOfBits::Record::getValue(unsigned fieldIndex, string& value) const
{
   assert(fieldTypes_[fieldIndex] == stringType);
   value = rawFields_[fieldIndex];
}

unsigned BucketOfBits::nRecords() const
{
   return records_.size();
}

unsigned BucketOfBits::nFields() const
{
   return fieldTypes_.size();
}

/** Returns the index of the field with the specified name.  Returns
 * nFields if the field name is not present in the bucket. */
unsigned BucketOfBits::getIndex(const std::string& fieldName) const
{
   for (unsigned ii=0; ii<fieldNames_.size(); ++ii)
      if (fieldName == fieldNames_[ii])
         return ii;
   return fieldNames_.size();
}


const std::string& BucketOfBits::fieldName(unsigned index) const
{
   return fieldNames_[index];
}

BucketOfBits::DataType BucketOfBits::dataType(unsigned index) const
{
   return fieldTypes_[index];
}

BucketOfBits::Record BucketOfBits::getRecord(unsigned index) const
{
   assert(index < records_.size());
   return Record(fieldTypes_, records_[index]);
}

void BucketOfBits::addRecord(const string& rec)
{
   records_.push_back(rec);
}


namespace
{
   string nextWord(const string& buf, string::size_type& pos)
   {
      static string whitespace(" \n\t");
      string::size_type begin = buf.find_first_not_of(whitespace, pos);
      pos = buf.find_first_of(whitespace, begin);
      return buf.substr(begin, pos);
   }
}

   
