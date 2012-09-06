#ifndef OBJECT_CC_HH
#define OBJECT_CC_HH

#include <string>
#include <vector>

#include "object.h"
#include "Long64.hh"

OBJECT* objectFind(const std::string& name, const std::string& objClass);

void objectGet(OBJECT* obj, const std::string& name, std::string& value, const std::string& defVal);
void objectGet(OBJECT* obj, const std::string& name, double& value,      const std::string& defVal);
void objectGet(OBJECT* obj, const std::string& name, bool& value,        const std::string& defVal);
void objectGet(OBJECT* obj, const std::string& name, int& value,         const std::string& defVal);
void objectGet(OBJECT* obj, const std::string& name, unsigned& value,    const std::string& defVal);
void objectGet(OBJECT* obj, const std::string& name, uint64_t& value,    const std::string& defVal);
void objectGet(OBJECT* obj, const std::string& name, int64_t&  value,    const std::string& defVal);

void objectGet(OBJECT* obj, const std::string& name, double& value, const std::string& defVal,
               const std::string& unitConvertTo);


void objectGet(OBJECT* obj, const std::string& name, std::vector<std::string>& value);
void objectGet(OBJECT* obj, const std::string& name, std::vector<int>&    value);
void objectGet(OBJECT* obj, const std::string& name, std::vector<unsigned>&    value);
void objectGet(OBJECT* obj, const std::string& name, std::vector<uint64_t>&      value);
void objectGet(OBJECT* obj, const std::string& name, std::vector<double >&     value);



#endif // #ifndef OBJECT_CC_HH
