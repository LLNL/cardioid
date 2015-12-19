#include "reactionFactory.hh"

#include <cassert>
#include <cstdlib>
#include <mpi.h>
#include <string.h>
#include "object_cc.hh"
#include "ioUtils.h"
#include "mpiUtils.h"
#include "OHaraRudy_Reaction.hh"    // OHara Rudy .
namespace  scanReaction
{
   char *getKeywordValue(char *ptr, char **keywordPtr, char **valuePtr)
   {
      char *end; 
      int size; 

      end = index(ptr,'='); 
      size = end-ptr; 
      char *keyword = (char *)malloc(size+1); 
      for (int i=0;i<size;i++) (keyword)[i] = ptr[i]; 
      keyword[size] = '\0';
      ptr = end + 1; 

      end = index(ptr,';'); 
      size = end-ptr; 
      char *value = (char *)malloc(size+1); 
      for (int i=0;i<size;i++) value[i] = ptr[i]; 
      value[size] = '\0';
      ptr = end +1;
      *keywordPtr = keyword; 
      *valuePtr = value; 
      return ptr; 
   }
   Reaction* scanOHaraRudy(OBJECT* obj, const Anatomy& anatomy)
   {
      const char *defaultCurrentNames[]={"INaCai","INaCass","INaK","INab","ICab","IKb","IpCa","INaFast","INaL","Ito","IKr","IKs","IK1","ICa"};
      OHaraRudy_Parms parms; 

      vector<string>currentNames(defaultCurrentNames,defaultCurrentNames+14); 
      int exists = object_testforkeyword(obj, "currents");
      if (exists) objectGet(obj,"currents",currentNames); 
      for (int ii=0;ii<currentNames.size();ii++)
      {
         string value; 
         string name = currentNames[ii]; 
         string defaultValue ="OHaraRudy"; 
         if (name == "INaFast")  defaultValue = "OHaraRudyMod"; 
         objectGet(obj,name,value,defaultValue); 
         printf("name=%s value=%s\n",name.c_str(),value.c_str());  fflush(stdout); 
         parms.currentNames.push_back(name); 
         parms.currentModels.push_back(value); 
      }
      Reaction *reaction = new OHaraRudy_Reaction(anatomy,parms);
         char *ptr = obj->value;
         while (*ptr != '\0')
         {
         char *keyword;
         char *valueString;
         ptr = getKeywordValue(ptr,&keyword,&valueString);
         assert(ptr != NULL);
         string  Keyword = keyword;
         int handle = reaction->getVarHandle(Keyword);
         if (handle != -1)
         {
         double value0;
         double value = atof(valueString);
         int iComp = handle >> 16;
         int iVar  = handle & 0xffff;
         value0=reaction->getValue(0,handle);
         printf("%8x %2d %2d %s %g %g\n",handle,iComp,iVar,keyword,value,value0);
         reaction->setValue(0,handle,value);
         value0=reaction->getValue(0,handle);
         printf("%8x %2d %2d %s %g %g\n",handle,iComp,iVar,keyword,value,value0);
         }
         free(keyword);
         free(valueString);
         }
      return  reaction; 
   }
}
