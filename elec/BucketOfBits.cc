#include "BucketOfBits.hh"

#include <cassert>
#include <cstdlib>
#include "ioUtils.h"

using namespace std;

namespace
{
   /** Returns the start position of the next field (or string::npos) in
    *  the string buf starting from beginPos (i.e., the start of the
    *  field after the field that beginPos points to).  A field consists
    *  of zero or more leading whitespace characers followed by one or
    *  more non-whitespace characters.  Hence, the next field is located
    *  by finding the first non-whitespace character after beginPos and
    *  then finding whilespace (or the end of buf). */
   string::size_type findNextField(const string& buf, string::size_type beginPos);
}

/** All fields are expected to have a name, a type, and a unit.  If a
 * field is dimensionless, or if there are no units provided with the
 * data being places in the bucket, use "1" as the unit for the
 * corresponding field. */
BucketOfBits::BucketOfBits(const vector<string>& fieldNames,
                           const vector<string>& fieldTypes,
                           const vector<string>& fieldUnits)
:fieldNames_(fieldNames),
 fieldUnits_(fieldUnits)
{
   assert(fieldNames.size() == fieldTypes.size());
   assert(fieldNames.size() == fieldUnits.size());

   for (unsigned ii=0; ii<fieldTypes.size(); ++ii)
   {
      if (fieldTypes[ii] == "s")
         fieldTypes_.push_back(stringType);
      else if (fieldTypes[ii] == "f")
         fieldTypes_.push_back(floatType);
      else if (fieldTypes[ii] == "u")
         fieldTypes_.push_back(intType);
      else if (fieldTypes[ii] == "d")   // ewd: d used by GradientVoronoiSensor, needs to be supported
         fieldTypes_.push_back(intType);
      else if (fieldTypes[ii] == "f4")
         fieldTypes_.push_back(f4Type);
      else if (fieldTypes[ii] == "f8")
         fieldTypes_.push_back(f8Type);
      else if (fieldTypes[ii] == "u8")
         fieldTypes_.push_back(u8Type);
      else
         assert(false);
   }
}


BucketOfBits::Record::Record(const vector<DataType>& fieldTypes,
                             const string& rawData)
: fieldTypes_(fieldTypes), rawRecord_(rawData)
{
   string::size_type nextField = 0;
   offsets_.reserve(fieldTypes_.size()+1);
   offsets_.push_back(nextField);
   for (unsigned ii=0; ii<fieldTypes_.size(); ++ii)
   {
      switch (fieldTypes_[ii])
      {
        case floatType:
        case intType:
        case stringType:
         nextField = findNextField(rawRecord_, offsets_[ii]);
         break;
        case u8Type:
        case f8Type:
         nextField += 8;
         break;
        case f4Type:
         nextField += 4;
         break;
        default:
         assert(false);
      }
      offsets_.push_back(nextField);
   }
}

const string& BucketOfBits::Record::getRawData() const
{
   return rawRecord_;
}

void BucketOfBits::Record::getValue(unsigned fieldIndex, double& value) const
{
   const char* startP = rawRecord_.data()+offsets_[fieldIndex];
   switch (fieldTypes_[fieldIndex])
   {
     case floatType:
      value = strtod(startP, NULL);
      break;
     case f8Type:
      value = mkDouble((const unsigned char*)startP, "f8");
      break;
     case f4Type:
      value = mkDouble((const unsigned char*)startP, "f4");
      break;
     default:
      assert(false);
   }
}

void BucketOfBits::Record::getValue(unsigned fieldIndex, int& value) const
{
   const char* startP = rawRecord_.data()+offsets_[fieldIndex];
   switch (fieldTypes_[fieldIndex])
   {
     case intType:
      value = strtol(startP, NULL, 10);
      break;
     case u8Type:
      value = mkInt((const unsigned char*)startP, "u8");
      break;
     default:
      assert(false);
   }
}

void BucketOfBits::Record::getValue(unsigned fieldIndex, uint64_t& value) const
{
   const char* startP = rawRecord_.data()+offsets_[fieldIndex];
   switch (fieldTypes_[fieldIndex])
   {
     case intType:
      value = strtoull(startP, NULL, 10);
      break;
     case u8Type:
      value = mkInt((const unsigned char*)startP, "u8");
      break;
     default:
      assert(false);
   }
}

void BucketOfBits::Record::getValue(unsigned fieldIndex, string& value) const
{
   assert(fieldTypes_[fieldIndex] == stringType);
   string::size_type len = offsets_[fieldIndex+1] - offsets_[fieldIndex];
   value = rawRecord_.substr(offsets_[fieldIndex], len);
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

const std::string& BucketOfBits::units(unsigned index) const
{
   return fieldUnits_[index];
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

void BucketOfBits::clearRecords()
{
   records_.clear();
}

namespace
{
   /** This function won't work for mixed binary and ascii records.  Of
    * course, we don't mix types like that so it isn't a problem.  I'm
    * not even sure how you could mix types since there wout be no way
    * to tell if a whitespace was the end of an ascii field or the start
    * of a binary. */
   string::size_type findNextField(const string& buf, string::size_type beginPos)
   {
      static string whitespace(" \n\t");
      string::size_type begin = buf.find_first_not_of(whitespace, beginPos);
      string::size_type endPos = buf.find_first_of(whitespace, begin);
      return endPos;
   }
}

   
