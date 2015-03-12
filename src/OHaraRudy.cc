#include "OHaraRudy.hh"
#include "OHaraRudy.h"
#include <cmath>
#include <cassert>

using namespace std;
#include "OHaraRudySetValue.hh" 
#include "OHaraRudyGetValue.hh" 
#include "OHaraRudyGetHandleMap.hh"

OHaraRudy::OHaraRudy(int cellType, Long64 gid)
{
   
   static int setup =0; 
   if ( !setup) OHaraRudySetup(); 
   setup=1; 
   cellParms_= OHaraRudyCellParms(cellType);
   OHaraRudyInitialState(cellType,&state_);
   gid_= gid; 
   defaultVoltage_ = state_.Vm;
}
/*
OHaraRudy::OHaraRudy(CELLPARMS cell, STATE state)
{
   
   cellParms_= new cell
   state_ = state; 
   OHaraRudyInitialState(cellType,&state);
   defaultVoltage_ = state_.Vm;
}
*/


double OHaraRudy::calc(double dt, double Vm, double iStim)
{
   STATE D; 
   state_.Vm = Vm;
   double dVm = OHaraRudyIntegrate(dt, iStim, &state_,cellParms_,&D); 
   return dVm;  
}
double OHaraRudy::calcS(double dt, double Vm, double iStim)
{
   STATE D; 
   state_.Vm = Vm;
   double dVm = OHaraRudyIntegrateS(dt, iStim, &state_,cellParms_,&D); 
   return dVm;  
}


double OHaraRudy::defaultVoltage()
{
   return defaultVoltage_;
}

void OHaraRudy::getCheckpointInfo(vector<string>& name,
      vector<string>& unit)
{
   const HandleMap& handleMap = getHandleMap();
   for (HandleMap::const_iterator
         iter=handleMap.begin(); iter!=handleMap.end(); ++iter)
   {
      if (iter->second.checkpoint_)
      {
         name.push_back(iter->first);
         unit.push_back(iter->second.unit_);
      }
   }
}

/* the handle representation.  Returns the value undefinedName for * unrecognized varName. */
int OHaraRudy::getVarHandle(const string& varName)
{
   return getHandleMap()[varName].handle_;
}

void OHaraRudy::getValue(const vector<int>& handle,
      vector<double>& value) const
{
   for (unsigned ii=0; ii<handle.size(); ++ii)
      value[ii] = getValue((VarHandle)handle[ii]);
}

const string OHaraRudy::getUnit(const string& varName)
{
   return getHandleMap()[varName].unit_;
}

void OHaraRudy::initConsts(int cellType)
{
   static bool initialized = false;
   if (! initialized)
   {
   initialized = true; 
   }
}

void OHaraRudy::initStates(int cellType)
{
   OHaraRudyInitialState(cellType,&state_); 
}

#ifdef SA
int main()
{
OHaraRudySetup(); 
OHaraRudy *cell = new OHaraRudy(0); 
double Vm =  -85; //cell->defaultVoltage_; 
double dV = cell->calc(0.01,Vm,0.0) ; 
printf("dV=%e\n",dV); 
}
#endif
