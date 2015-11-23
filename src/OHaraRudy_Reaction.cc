#include "OHaraRudy_Reaction.hh"
#include <cmath>
#include <string>
#include <map>
#include <cstdio>
#include "Anatomy.hh"
#include "OHaraRudy.hh"
#include "OHaraRudy.h"
#include "BucketOfBits.hh"
#include "units.h"

using namespace std;
HandleMap  OHaraRudy_Reaction::handleMap_ ;

OHaraRudy_Reaction::OHaraRudy_Reaction(const Anatomy& anatomy,OHaraRudy_Parms &parms)
{
   OHaraRudyInit(0.0,anatomy.nLocal());
   makeHandleMap(); 
   ttType_.resize(256, -1); 
   ttType_[30] = ENDO_CELL;
   ttType_[31] = ENDO_CELL;
   ttType_[75] = ENDO_CELL;
   ttType_[76] = M_CELL;
   ttType_[77] = EPI_CELL;
   ttType_[100] = ENDO_CELL;
   ttType_[101] = M_CELL;
   ttType_[102] = EPI_CELL;
   indexS_=-2; 

   cells_.reserve(anatomy.nLocal());
   int ttType0 = ttType_[anatomy.cellType(0)];
   for (unsigned ii=0; ii<anatomy.nLocal(); ++ii)
   {
      assert(anatomy.cellType(ii) >= 0 && anatomy.cellType(ii) < 256);
      int ttType = ttType_[anatomy.cellType(ii)];
      Long64 gid = anatomy.gid(ii);
      if ( gid ==  (Long64)(-1)) 
         cells_.push_back(OHaraRudyDebug(ttType,parms));
      else
         cells_.push_back(OHaraRudy(ttType,parms));
   }
}

void OHaraRudy_Reaction::calc(double dt,
      const VectorDouble32& Vm,
      const vector<double>& iStim,
      VectorDouble32& dVm)
{
   assert(cells_.size() == dVm.size());

   int nCells = cells_.size();
//   OHaraRudyPut(nCells,&Vm[0],&iStim[0]); 
//   OHaraRudyCalc(); 
//   OHaraRudyGet(&dVm[0]); 
   for (int ii=0; ii<nCells; ++ii)
      dVm[ii] = cells_[ii].calc(dt, Vm[ii], iStim[ii]);
}

void OHaraRudy_Reaction::initializeMembraneVoltage(VectorDouble32& Vm)
{
   assert(Vm.size() >= cells_.size());
   for (unsigned ii=0; ii<cells_.size(); ++ii)
      Vm[ii] = cells_[ii].defaultVoltage();
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
         if (  name != "vm" ) checkpoint = true; 
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
   cells_[iCell].setValue(varHandle, value);
   OHaraRudySetValue(iCell, varHandle, value);
}

double OHaraRudy_Reaction::getValue(int iCell, int varHandle) const
{
   return cells_[iCell].getValue(varHandle);
}

void OHaraRudy_Reaction::getValue(int iCell, const vector<int>& handle, vector<double>& value) const
{

   cells_[iCell].getValue(handle, value);
}


