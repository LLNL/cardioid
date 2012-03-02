#include "TT06_CellML_Reaction.hh"
#include <cmath>
#include "Anatomy.hh"
#include "TT06_CellML.hh"
#include "TT06_CellML_Endo.hh"
#include "TT06_CellML_Mid.hh"
#include "TT06_CellML_Epi.hh"

using namespace std;

struct TT06_CellMLState
{
   double state[19];
};




TT06_CellML_Reaction::TT06_CellML_Reaction(const Anatomy& anatomy,
                                           IntegratorType integrator)
: nCells_(anatomy.nLocal()),
          integrator_(integrator)
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

   
   cellModel_.reserve(nCells_);
   for (int ii=0; ii<nCells_; ++ii)
   {
      assert(anatomy.cellType(ii) >= 0 && anatomy.cellType(ii) < 256);
      int ttType = ttType_[anatomy.cellType(ii)];
      switch (ttType)
      {
        case 0:
         cellModel_.push_back(new TT06_CellML_Endo());
         break;
        case 1:
         cellModel_.push_back(new TT06_CellML_Mid());
         break;
        case 2:
         cellModel_.push_back(new TT06_CellML_Epi());
         break;
        default:
         assert(false);
      }
   }

   s_.resize(nCells_);
   for (int ii=0; ii<nCells_; ++ii)
   {
      for (unsigned jj=0; jj<19; ++jj)
      {
         s_[ii].state[jj] = cellModel_[ii]->defaultState(jj);
      }
   }
   
   
}

TT06_CellML_Reaction::~TT06_CellML_Reaction()
{
   for (unsigned ii=0; ii<cellModel_.size(); ++ii)
   {
      delete cellModel_[ii];
   }
}

void TT06_CellML_Reaction::calc(double dt,
                                const vector<double>& Vm,
                                const vector<double>& iStim,
                                vector<double>& dVm)
{
   assert(nCells_ == dVm.size());

   switch (integrator_)
   {
     case forwardEuler:
      forwardEulerIntegrator(dt, Vm, iStim, dVm);
      break;
     case rushLarsen:
      rushLarsenIntegrator(dt, Vm, iStim, dVm);
      break;
     default:
      assert(false);
   }
}

void TT06_CellML_Reaction::initializeMembraneVoltage(std::vector<double>& Vm)
{
   assert(Vm.size() >= s_.size());
   for (unsigned ii=0; ii<s_.size(); ++ii)
      Vm[ii] = cellModel_[ii]->defaultState(0);
}

void TT06_CellML_Reaction::getCheckpointInfo(vector<string>& name,
                                             vector<string>& unit) const
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


/** This function maps the string representation of a variable name to
 *  the handle representation.  Returns the value undefinedName for
 *  unrecognized varName. */
int TT06_CellML_Reaction::getVarHandle(const string& varName) const
{
   return getHandleMap()[varName].handle_;
}

/** This function maps the string representation of a variable name to
 *  the handle representation.  Returns the value undefinedName for
 *  unrecognized varName. */
vector<int> TT06_CellML_Reaction::getVarHandle(const vector<string>& varName) const
{
   vector<int> handle;
   for (unsigned ii=0; ii<varName.size(); ++ii)
      handle.push_back(getHandleMap()[varName[ii]].handle_);

   return handle;
}
      
void TT06_CellML_Reaction::setValue(int iCell, int varHandle, double value)
{
//ddt   cout << "Setting var "<< varHandle << "=" << value<<endl;
   
   switch (varHandle)
   {
     case undefinedName:
      assert(false);
      break;
     case Vm:          s_[iCell].state[0] = value;      break;
     case K_i:         s_[iCell].state[1] = value;      break;
     case Na_i:        s_[iCell].state[2] = value;      break;
     case Ca_i:        s_[iCell].state[3] = value;      break;
     case Xr1_gate:    s_[iCell].state[4] = value;      break;
     case Xr2_gate:    s_[iCell].state[5] = value;      break;
     case Xs_gate:     s_[iCell].state[6] = value;      break;
     case m_gate:      s_[iCell].state[7] = value;      break;
     case h_gate:      s_[iCell].state[8] = value;      break;
     case j_gate:      s_[iCell].state[9] = value;      break;
     case Ca_ss:       s_[iCell].state[10] = value;     break;
     case d_gate:      s_[iCell].state[11] = value;     break;
     case f_gate:      s_[iCell].state[12] = value;     break;
     case f2_gate:     s_[iCell].state[13] = value;     break;
     case fCass_gate:  s_[iCell].state[14] = value;     break;
     case s_gate:      s_[iCell].state[15] = value;     break;
     case r_gate:      s_[iCell].state[16] = value;     break;
     case Ca_SR:       s_[iCell].state[17] = value;     break;
     case R_prime:     s_[iCell].state[18] = value;     break;
     case nVars:
      assert(false);
      break;
   }
}

double TT06_CellML_Reaction::getValue(int iCell, int handle) const
{
   switch (handle)
   {
     case undefinedName:
      assert(false);
      break;
     case Vm:         return s_[iCell].state[0];      break;
     case K_i:        return s_[iCell].state[1];      break;
     case Na_i:       return s_[iCell].state[2];      break;
     case Ca_i:       return s_[iCell].state[3];      break;
     case Xr1_gate:   return s_[iCell].state[4];      break;
     case Xr2_gate:   return s_[iCell].state[5];      break;
     case Xs_gate:    return s_[iCell].state[6];      break;
     case m_gate:     return s_[iCell].state[7];      break;
     case h_gate:     return s_[iCell].state[8];      break;
     case j_gate:     return s_[iCell].state[9];      break;
     case Ca_ss:      return s_[iCell].state[10];     break;
     case d_gate:     return s_[iCell].state[11];     break;
     case f_gate:     return s_[iCell].state[12];     break;
     case f2_gate:    return s_[iCell].state[13];     break;
     case fCass_gate: return s_[iCell].state[14];     break;
     case s_gate:     return s_[iCell].state[15];     break;
     case r_gate:     return s_[iCell].state[16];     break;
     case Ca_SR:      return s_[iCell].state[17];     break;
     case R_prime:    return s_[iCell].state[18];     break;
     case nVars:
      assert(false);
      break;
   }
   return 0.;
}

void TT06_CellML_Reaction::getValue(int iCell,
                                    const vector<int>& handle,
                                    vector<double>& value) const
{
   for (unsigned ii=0; ii<handle.size(); ++ii)
      value[ii] = getValue(iCell, (VarHandle)handle[ii]);
}


const string TT06_CellML_Reaction::getUnit(const string& varName) const
{
   return getHandleMap()[varName].unit_;
}


/** Remember that down in the cell models the units don't necessarily
 *  correspond to the internal units of Cardioid.  The units in this map
 *  are the units the cell model expects the variables to have. */
TT06_CellML_Reaction::HandleMap& TT06_CellML_Reaction::getHandleMap()
{
   static HandleMap handleMap;
   if (handleMap.size() == 0)
   {
      handleMap["Vm"]          = VarInfo(Vm,         false, "mV");
      handleMap["K_i"]         = VarInfo(K_i,        true,  "mM");
      handleMap["Na_i"]        = VarInfo(Na_i,       true,  "mM");
      handleMap["Ca_i"]        = VarInfo(Ca_i,       true,  "mM");
      handleMap["Xr1_gate"]    = VarInfo(Xr1_gate,   true,  "1");
      handleMap["Xr2_gate"]    = VarInfo(Xr2_gate,   true,  "1");
      handleMap["Xs_gate"]     = VarInfo(Xs_gate,    true,  "1");
      handleMap["m_gate"]      = VarInfo(m_gate,     true,  "1");
      handleMap["h_gate"]      = VarInfo(h_gate,     true,  "1");
      handleMap["j_gate"]      = VarInfo(j_gate,     true,  "1");
      handleMap["Ca_ss"]       = VarInfo(Ca_ss,      true,  "mM");
      handleMap["d_gate"]      = VarInfo(d_gate,     true,  "1");
      handleMap["f_gate"]      = VarInfo(f_gate,     true,  "1");
      handleMap["f2_gate"]     = VarInfo(f2_gate,    true,  "1");
      handleMap["fCass_gate"]  = VarInfo(fCass_gate, true,  "1");
      handleMap["s_gate"]      = VarInfo(s_gate,     true,  "1");
      handleMap["r_gate"]      = VarInfo(r_gate,     true,  "1");
      handleMap["Ca_SR"]       = VarInfo(Ca_SR,      true,  "mM");
      handleMap["R_prime"]     = VarInfo(R_prime,    true,  "1");
      assert(handleMap.size() == nVars-1);
   }
   return handleMap;
}

void TT06_CellML_Reaction::forwardEulerIntegrator(
   double dt,
   const vector<double>& Vm,
   const vector<double>& iStim,
   vector<double>& dVm)
{
#pragma omp parallel for
   for (int ii=0; ii<nCells_; ++ii)
   {
      double rates[19];
      double algebraic[70];
      dVm[ii] = cellModel_[ii]->calc(Vm[ii], iStim[ii], s_[ii].state, rates, algebraic);

      // forward euler to integrate internal state variables.
      for (unsigned jj=1; jj<19; ++jj)
         s_[ii].state[jj] += rates[jj] * dt;
   
   }
}

void TT06_CellML_Reaction::rushLarsenIntegrator(
   double dt,
   const vector<double>& Vm,
   const vector<double>& iStim,
   vector<double>& dVm)
{
#pragma omp parallel for
   for (int ii=0; ii<nCells_; ++ii)
   {
      double rates[19];
      double algebraic[70];
      dVm[ii] = cellModel_[ii]->calc(Vm[ii], iStim[ii], s_[ii].state, rates, algebraic);

      // forward euler for all states except rushLarsen for fast sodium m gate.
      for (unsigned jj=1; jj<7; ++jj)
         s_[ii].state[jj] += rates[jj] * dt;
      s_[ii].state[7] =  algebraic[3] - (algebraic[3]-s_[ii].state[7])*exp(-dt/algebraic[37]);
      for (unsigned jj=8; jj<19; ++jj)
         s_[ii].state[jj] += rates[jj] * dt;
   
   }
}

