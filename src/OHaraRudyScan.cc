#include "reactionFactory.hh"

#include <cassert>
#include <cstdlib>
#include <mpi.h>
#include "object_cc.hh"
#include "ioUtils.h"
#include "mpiUtils.h"
#include "OHaraRudy_Reaction.hh"    // OHara Rudy .
#include "OHaraRudy.hh"    // OHara Rudy .
Reaction* scanOHaraRudy(OBJECT* obj, const Anatomy& anatomy)
{
   OHaraRudy_Parms parms; 
   int nCurrent=0; 
   std::string name=OHaraRudyCurrentNames[nCurrent];
   string defaultValue; 
   while(name != "")
   {
      string value; 
      defaultValue ="OHaraRudy"; 
      if (name == "INaFast")  defaultValue = "OHaraRudyMod"; 
      objectGet(obj,name,value,defaultValue); 
      printf("%s %s\n",name.c_str(),value.c_str()); 
      parms.currentNames.push_back(name); 
      parms.currentModels.push_back(value); 
      name = OHaraRudyCurrentNames[++nCurrent]; 
   }
   Reaction *reaction = new OHaraRudy_Reaction(anatomy,parms);
   return  reaction; 
}
