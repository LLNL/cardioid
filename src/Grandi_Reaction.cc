#include "Grandi_Reaction.hh"
#include <cmath>
#include <string>
#include <map>
#include <cstdio>
#include "Anatomy.hh"
#include "Grandi.h"
#include "BucketOfBits.hh"
#include "units.h"

using namespace std;
HandleMap  Grandi_Reaction::handleMap_ ;

Grandi_Reaction::Grandi_Reaction(const Anatomy& anatomy,Grandi_Parms &parms)
{
   ttType_.resize(256, -1); 
   ttType_[100] = RA_AF;
   ttType_[101] = LA_AF;
   ttType_[102] = RA_SR;
   ttType_[103] = LA_SR;
   indexS_=-2; 
   nCells_ = anatomy.nLocal(); 
   int cellType[nCells_]; 
   for (unsigned ii=0; ii<nCells_; ++ii)
   {
      assert(anatomy.cellType(ii) >= 0 && anatomy.cellType(ii) < 256);
      cellType[ii] = ttType_[anatomy.cellType(ii)];
   }
   GrandiInit(0.0,nCells_,cellType);
   makeHandleMap(); 

}

void Grandi_Reaction::calc(double dt, const VectorDouble32& Vm,
      const vector<double>& iStim, VectorDouble32& dVm)
{
   assert(nCells_ == dVm.size());

   GrandiPut(dt,nCells_,&Vm[0],&iStim[0]); 
   GrandiCalc(); 
   GrandiGet(&dVm[0]); 
}

void Grandi_Reaction::initializeMembraneVoltage(VectorDouble32& Vm)
{
   assert(Vm.size() >= nCells_);
   int varHandle = Grandi_Reaction::handleMap_["Vm"].handle_; 
   for (unsigned ii=0; ii<nCells_; ++ii)
   {
      Vm[ii] = GrandiGetValue(ii,varHandle); 
   }
}

void Grandi_Reaction::makeHandleMap()
{
   int ncomp=GrandiGet_nComp(); 
   COMPONENTINFO* info=GrandiGet_compInfo(); 
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
         Grandi_Reaction::handleMap_[name] = CheckpointVarInfo(index, checkpoint, units );
      }
   } 
}
int Grandi_Reaction::getVarHandle(const string& varName) const
{
   return Grandi_Reaction::handleMap_[varName].handle_; 
}
void Grandi_Reaction::getCheckpointInfo(vector<string>& name, vector<string>& unit) const
{
   const HandleMap& handleMap = Grandi_Reaction::handleMap_;
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
const string Grandi_Reaction::getUnit(const string& varName) const 
{
   return Grandi_Reaction::handleMap_[varName].unit_;
}

void Grandi_Reaction::setValue(int iCell, int varHandle, double value)
{
   GrandiSetValue(iCell, varHandle, value);
}

double Grandi_Reaction::getValue(int iCell, int varHandle) const
{
   return GrandiGetValue(iCell,varHandle); 
}

void Grandi_Reaction::getValue(int iCell, const vector<int>& handle, vector<double>& value) const
{
   for (int ii=0;ii<handle.size();ii++) value[ii]= GrandiGetValue(iCell,handle[ii]); 

}
