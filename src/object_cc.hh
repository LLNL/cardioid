#ifndef OBJECT_CC_HH
#define OBJECT_CC_HH

#include <string>
#include <vector>

#include "object.h"

void objectGet(OBJECT* obj, const std::string& name, std::string& value, const std::string& defVal);
void objectGet(OBJECT* obj, const std::string& name, double& value,      const std::string& defVal);
void objectGet(OBJECT* obj, const std::string& name, int& value,         const std::string& defVal);
void objectGet(OBJECT* obj, const std::string& name, unsigned& value,    const std::string& defVal);


void objectGet(OBJECT* obj, const std::string& name, std::vector<std::string>& value);
void objectGet(OBJECT* obj, const std::string& name, std::vector<unsigned>&    value);
void objectGet(OBJECT* obj, const std::string& name, std::vector<double >&     value);



#endif // #ifndef OBJECT_CC_HH
