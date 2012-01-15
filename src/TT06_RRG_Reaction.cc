#include "TT06_RRG_Reaction.hh"
#include <cmath>
#include <string>
#include <map>
#include "Anatomy.hh"
#include "TT06_RRG.hh"
#include "BucketOfBits.hh"

using namespace std;

TT06_RRG_Reaction::TT06_RRG_Reaction(const Anatomy& anatomy)
{
   ttType_.resize(256, -1); 
   ttType_[30] = 0;
   ttType_[31] = 0;
   ttType_[75] = 0;
   ttType_[76] = 1;
   ttType_[77] = 2;
   ttType_[100] = 0;
   ttType_[101] = 1;
   ttType_[102] = 2;

   cells_.reserve(anatomy.nLocal());
   for (unsigned ii=0; ii<anatomy.nLocal(); ++ii)
   {
      assert(anatomy.cellType(ii) >= 0 && anatomy.cellType(ii) < 256);
      int ttType = ttType_[anatomy.cellType(ii)];
      cells_.push_back(TT06_RRG(ttType));
   }
}

void TT06_RRG_Reaction::calc(double dt,
                             const vector<double>& Vm,
                             const vector<double>& iStim,
                             vector<double>& dVm)
{
   assert(cells_.size() == dVm.size());

   int nCells = cells_.size();
#pragma omp parallel for
   for (int ii=0; ii<nCells; ++ii)
      dVm[ii] = cells_[ii].calc(dt, Vm[ii], iStim[ii]);
}

void TT06_RRG_Reaction::initializeMembraneVoltage(std::vector<double>& Vm)
{
   assert(Vm.size() >= cells_.size());
   for (unsigned ii=0; ii<cells_.size(); ++ii)
      Vm[ii] = cells_[ii].defaultVoltage();
}

void TT06_RRG_Reaction::loadState(const BucketOfBits& data)
{
   assert(cells_.size() == data.nRecords());

   typedef map<int, TT06_RRG::VarHandle> FieldMap;
   FieldMap fieldMap;

   for (unsigned ii=0; ii<data.nFields(); ++ii)
   {
      TT06_RRG::VarHandle handle = TT06_RRG::getVarHandle(data.fieldName(ii));
      if (handle != TT06_RRG::undefinedName)
         fieldMap[ii] = handle;
   }
   
   for (unsigned ii=0; ii<cells_.size(); ++ii)
   {
      BucketOfBits::Record iRec = data.getRecord(ii);
      for (FieldMap::const_iterator iter=fieldMap.begin();
           iter!=fieldMap.end(); ++iter)
      {
         int iField = iter->first;
         TT06_RRG::VarHandle handle = iter->second;
         double value;
         switch (data.dataType(iField))
         {
           case BucketOfBits::floatType:
            iRec.getValue(iField, value);
            break;
           case BucketOfBits::intType:
            int tmp;
            iRec.getValue(iField, tmp);
            value = double(tmp);
            break;
           default:
            assert(false);
         }
         cells_[ii].setVariable(handle, value);
      }
   }
}
