#include "OHaraRudy_Reaction.hh"
#include <cmath>
#include <string>
#include <map>
#include <cstdio>
#include "Anatomy.hh"
#include "OHaraRudy.h"
#include "BucketOfBits.hh"
#include "units.h"

using namespace std;
HandleMap  OHaraRudy_Reaction::handleMap_ ;

OHaraRudy_Reaction::OHaraRudy_Reaction(const int numPoints,OHaraRudy_Parms &parms)
{
   nCells_ = numPoints;
   cellType = parms.cellType;
   OHaraRudyInit(0.0,nCells_,cellType);
   makeHandleMap(); 

}

void OHaraRudy_Reaction::calc(double dt, const VectorDouble32& Vm,
      const vector<double>& iStim, VectorDouble32& dVm)
{
   assert(nCells_ == dVm.size());

   OHaraRudyPut(dt,nCells_,&Vm[0],&iStim[0]); 
   OHaraRudyCalc(); 
   OHaraRudyGet(&dVm[0]); 
}

void OHaraRudy_Reaction::initializeMembraneVoltage(VectorDouble32& Vm)
{
   assert(Vm.size() >= nCells_);
   int varHandle = OHaraRudy_Reaction::handleMap_["Vm"].handle_; 
   for (unsigned ii=0; ii<nCells_; ++ii)
   {
      Vm[ii] = OHaraRudyGetValue(ii,varHandle); 
   }
}

void OHaraRudy_Reaction::makeHandleMap()
{
   int ncomp=OHaraRudyGet_nComp(); 
   COMPONENTINFO* info=OHaraRudyGet_compInfo(); 
   for (int i=0;i<ncomp;i++)
   {
      VARINFO *varinfo = info[i].varInfo; 
      for (int j=0;j<info[i].nVar;j++) 
      {
         string name = varinfo[j].name ;
         string units = varinfo[j].units;
         unsigned int index = (i << 16) + varinfo[j].index;
         bool  checkpoint = false;
         if (  name != "Vm" ) checkpoint = true; 
         OHaraRudy_Reaction::handleMap_[name] = CheckpointVarInfo(index, checkpoint, units );
      }
   } 
}
int OHaraRudy_Reaction::getVarHandle(const string& varName) const
{
   return OHaraRudy_Reaction::handleMap_[varName].handle_; 
}
void OHaraRudy_Reaction::getCheckpointInfo(vector<string>& name, vector<string>& unit) const
{
   const HandleMap& handleMap = OHaraRudy_Reaction::handleMap_;
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
const string OHaraRudy_Reaction::getUnit(const string& varName) const 
{
   return OHaraRudy_Reaction::handleMap_[varName].unit_;
}

void OHaraRudy_Reaction::setValue(int iCell, int varHandle, double value)
{
   OHaraRudySetValue(iCell, varHandle, value);
}

double OHaraRudy_Reaction::getValue(int iCell, int varHandle) const
{
   return OHaraRudyGetValue(iCell,varHandle); 
}

void OHaraRudy_Reaction::getValue(int iCell, const vector<int>& handle, vector<double>& value) const
{
   for (int ii=0;ii<handle.size();ii++) value[ii]= OHaraRudyGetValue(iCell,handle[ii]); 

}


