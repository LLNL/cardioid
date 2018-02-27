#include "ExcitationContraction.hpp"

using namespace std;

std::vector<int> XXX::getHandles(const std::vector<std::string>& varNames) const
{
   std::vector<int> varHandles(varNames.size());
   for (int ii=0; ii<varHandles.size(); ii++)
   {
      varHandles[ii] = getHandle(varNames[ii]);
   }
   return varHandles;
}


bool XXX::getOrder(const std::string type, const std::vector<std::string>& inputNames, std::vector<int>& permutation) const
{
   vector<string> varNames = getVarGroup(type);
   bool allOK = true;
   if (varNames.size() != inputNames.size()) { allOK = false; }
   permutation.resize(inputNames.size(), -1);
   for (int iname=0; iname < varNames.size(); iname++)
   {
      for (int iinput=0; iinput < inputNames.size(); iinput++)
      {
         if (inputNames[iinput] == varNames[iname])
         {
            permutation[iinput] = iname;
         }
      }
   }
   for (int iinput=0; iinput < inputNames.size(); iinput++)
   {
      if (permutation[iinput] == -1)
      {
         allOK = false;
      }
   }
   return allOK;
}
