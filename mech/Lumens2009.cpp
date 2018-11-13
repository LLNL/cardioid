#include "Lumens2009.hpp"
#include <cassert>
#include <cmath>

using namespace std;

namespace Lumens2009
{

string ThisModel::name() const
{
   return "Lumens2009";
}

enum varHandles
{

   _handle_C,
   _handle_Lsc,
   _handle_actTime,
   _handle_celltype,
   _handle_dtension_dstretchVel,
   _handle_stretch,
   _handle_stretchVel,
   _handle_tension,
   NUM_HANDLES
};

const string varNames[] = {

   "C",
   "Lsc",
   "actTime",
   "celltype",
   "dtension_dstretchVel",
   "stretch",
   "stretchVel",
   "tension",
   ""
};

const string varUnits[] = {

   "not_implemented",
   "not_implemented",
   "not_implemented",
   "not_implemented",
   "not_implemented",
   "not_implemented",
   "not_implemented",
   "not_implemented",
   ""
};

enum inputOrder
{

   _inIdx_actTime,
   _inIdx_stretch,
   _inIdx_stretchVel,
   NUM_INPUTS
};

enum outputOrder
{

   _outIdx_dtension_dstretchVel,
   _outIdx_tension,
   NUM_OUTPUTS
};

int ThisModel::getHandle(const string& varName) const
{
   for (int ii=0; ii<NUM_HANDLES; ii++)
   {
      if (varName == varNames[ii])
      {
         return ii;
      }
   }
   return -1;
}

string ThisModel::getUnit(const int varHandle) const
{
   assert(varHandle > 0 && varHandle < NUM_HANDLES);
   return varUnits[varHandle];
}

vector<string> ThisModel::getVarGroup(const std::string type) const
{
   vector<string> names;
   if (0) {}
 
   else if (type == "checkpoint")
   {
      names.reserve(2);
      names.push_back("C");
      names.push_back("Lsc");
   }
   else if (type == "flag")
   {
      names.reserve(1);
      names.push_back("celltype");
   }
   else if (type == "input")
   {
      names.reserve(3);
      names.push_back("actTime");
      names.push_back("stretch");
      names.push_back("stretchVel");
   }
   else if (type == "output")
   {
      names.reserve(2);
      names.push_back("dtension_dstretchVel");
      names.push_back("tension");
   }
   else if (type == "param")
   {
      names.reserve(0);
   }
   return names;
}


double ThisModel::get(const int varHandle) const
{
   if (0) {}

   else if (varHandle == _handle_celltype)
   {
      return celltype;
   }
   return NAN;
}
void ThisModel::set(const int varHandle, const double value)
{
   if (0) {}

   else if (varHandle == _handle_celltype)
   {
      celltype = value;
   }
   assert(0 && "Can't set a value for parameter that doesn't exist");
}

double ThisModel::get(const int varHandle, const int iCell) const
{
   if (0) {}

   else if (varHandle == _handle_C)
   {
      return _state[_old()][iCell].C;
   }
   else if (varHandle == _handle_Lsc)
   {
      return _state[_old()][iCell].Lsc;
   }
   return get(varHandle);
}

void ThisModel::set(const int varHandle, const int iCell, const double value)
{
   if (0) {}

   else if (varHandle == _handle_C)
   {
      _state[_old()][iCell].C = value;
   }
   else if (varHandle == _handle_Lsc)
   {
      _state[_old()][iCell].Lsc = value;
   }
   set(varHandle, value);
}

double ThisModel::get(const int varHandle, const int iCell, const double* inputs[]) const
{

   if (0) {}
   return get(varHandle, iCell);
}

ThisModel::ThisModel(const int numPoints)
{
   _numPoints = numPoints;
   _oldIdx = 0;
   _state[0].resize(numPoints);
   _state[1].resize(numPoints);

   //set the flags

   celltype = 2;
   //set the default parameters
   resetDefaultParameters();
}

void ThisModel::resetDefaultParameters()
{

}

void ThisModel::initialize(const double* const _inputs[])
{
   for (int _icell=0; _icell < _numPoints; _icell++)
   {
      const double actTime = _inputs[_inIdx_actTime][_icell];
      const double stretch = _inputs[_inIdx_stretch][_icell];
      const double stretchVel = _inputs[_inIdx_stretchVel][_icell];
      double C_rest = 0.0200000000000000;
      double L_s0 = 1.51000000000000;
      double Lsc_init = L_s0;
      double Lsc = Lsc_init;
      double C_init = C_rest;
      double C = C_init;
      _state[_old()][_icell].C = C;
      _state[_old()][_icell].Lsc = Lsc;
   }
}

const int internalTimestep = 20;

void ThisModel::tryTimestep(const double dt, const double* const _inputs[])
{
   for (int _icell=0; _icell < _numPoints; _icell++)
   {
      const double actTime = _inputs[_inIdx_actTime][_icell];
      const double stretchVel = _inputs[_inIdx_stretchVel][_icell];
      double stretch = _inputs[_inIdx_stretch][_icell];
      double _partial_C = 0;
      double _partial_Lsc = 0;
      double C = _state[_old()][_icell].C;
      double Lsc = _state[_old()][_icell].Lsc;
      for (int itime=0; itime<internalTimestep; itime++)
      {
         double _partial_stretch = itime*(dt/internalTimestep);
         double __melodee_temp_001 = celltype == 1;
         double sigma_act;
         double tau_d;
         double tau_r;
         double tau_sc;
         if (__melodee_temp_001)
         {
            tau_r = 48;
            tau_d = 32;
            tau_sc = 425;
         }
         else
         {
            double __melodee_temp_000 = celltype == 2;
            if (__melodee_temp_000)
            {
               tau_r = 48;
               tau_d = 32;
               tau_sc = 425;
            }
            else
            {
               tau_r = 28.1000000000000;
               tau_d = 33.8000000000000;
               tau_sc = 292.500000000000;
            }
         }
         double C_rest = 0.0200000000000000;
         double L_s0 = 1.51000000000000;
         double L_serel = 0.0400000000000000;
         double ls_unloaded = 2;
         double v_max = 0.00700000000000000;
         double ls = ls_unloaded*stretch;
         double _d_ls_wrt_stretchVel = _partial_stretch*ls_unloaded;
         double C_l = tanh(4.0*((-L_s0 + Lsc)*(-L_s0 + Lsc)));
         double _d_C_l_wrt_stretchVel = 8.0*_partial_Lsc*(-L_s0 + Lsc)*(-(tanh(4.0*((-L_s0 + Lsc)*(-L_s0 + Lsc)))*tanh(4.0*((-L_s0 + Lsc)*(-L_s0 + Lsc)))) + 1);
         double T = tau_sc*(0.3*Lsc + 0.29);
         double _d_T_wrt_stretchVel = 0.3*_partial_Lsc*tau_sc;
         double x = actTime/tau_r;
         double __melodee_temp_002 = x > 8;
         double x_001;
         if (__melodee_temp_002)
         {
            x_001 = 8;
         }
         else
         {
            x_001 = x;
         }
         double f_rise = 0.02*(x_001*x_001*x_001)*((-x_001 + 8.0)*(-x_001 + 8.0))*exp(-x_001);
         double l_snorm = (-Lsc + ls)/L_serel;
         double _d_l_snorm_wrt_stretchVel = (_d_ls_wrt_stretchVel - _partial_Lsc)/L_serel;
         double Lsc_diff = v_max*(l_snorm - 1);
         double _d_Lsc_diff_wrt_stretchVel = _d_l_snorm_wrt_stretchVel*v_max;
         double C_diff = C_l*f_rise/tau_r + (-C + C_rest)/(tau_d*(exp((T - actTime)/tau_d) + 1));
         double _d_C_diff_wrt_stretchVel = _d_C_l_wrt_stretchVel*f_rise/tau_r - _d_T_wrt_stretchVel*(-C + C_rest)*exp((T - actTime)/tau_d)/((tau_d*tau_d)*((exp((T - actTime)/tau_d) + 1)*(exp((T - actTime)/tau_d) + 1))) - _partial_C/(tau_d*(exp((T - actTime)/tau_d) + 1));
         C += C_diff*(dt/internalTimestep);
         Lsc += Lsc_diff*(dt/internalTimestep);
         _partial_C += _d_C_diff_wrt_stretchVel*(dt/internalTimestep);
         _partial_Lsc += _d_Lsc_diff_wrt_stretchVel*(dt/internalTimestep);
         stretch += stretchVel*(dt/internalTimestep);
      }
      _state[_new()][_icell].C = C;
      _state[_new()][_icell].Lsc = Lsc;
      _state[_new()][_icell]._partial_C = _partial_C;
      _state[_new()][_icell]._partial_Lsc = _partial_Lsc;
  }
}

void ThisModel::outputForNextTimestep(const double dt, const double* const _inputs[], double* const outputs[])
{
   for (int _icell=0; _icell < _numPoints; _icell++)
   {
      const double actTime = _inputs[_inIdx_actTime][_icell];
      const double stretchVel = _inputs[_inIdx_stretchVel][_icell];
      const double C = _state[_new()][_icell].C;
      const double Lsc = _state[_new()][_icell].Lsc;
      const double _partial_C = _state[_new()][_icell]._partial_C;
      const double _partial_Lsc = _state[_new()][_icell]._partial_Lsc;
      const double stretch = _inputs[_inIdx_stretch][_icell]+dt*stretchVel;
      const double _partial_stretch = dt;
      double __melodee_temp_001 = celltype == 1;
      double sigma_act;
      double tau_d;
      double tau_r;
      double tau_sc;
      if (__melodee_temp_001)
      {
         sigma_act = 100;
      }
      else
      {
         double __melodee_temp_000 = celltype == 2;
         if (__melodee_temp_000)
         {
            sigma_act = 120;
         }
         else
         {
            sigma_act = 60;
         }
      }
      double L_s0 = 1.51000000000000;
      double L_serel = 0.0400000000000000;
      double ls_unloaded = 2;
      double ls = ls_unloaded*stretch;
      double _d_ls_wrt_stretchVel = _partial_stretch*ls_unloaded;
      double l_snorm = (-Lsc + ls)/L_serel;
      double _d_l_snorm_wrt_stretchVel = (_d_ls_wrt_stretchVel - _partial_Lsc)/L_serel;
      double __melodee_temp_004 = -L_s0 + Lsc < 0;
      double tension;
      double _d_tension_wrt_stretchVel;
      if (__melodee_temp_004)
      {
         tension = 0;
         _d_tension_wrt_stretchVel = 0;
      }
      else
      {
         double __melodee_temp_003 = C < 0;
         if (__melodee_temp_003)
         {
            tension = 0;
            _d_tension_wrt_stretchVel = 0;
         }
         else
         {
            tension = C*l_snorm*sigma_act*(-L_s0 + Lsc);
            _d_tension_wrt_stretchVel = C*_d_l_snorm_wrt_stretchVel*sigma_act*(-L_s0 + Lsc) + C*_partial_Lsc*l_snorm*sigma_act + _partial_C*l_snorm*sigma_act*(-L_s0 + Lsc);
         }
      }
      double dtension_dstretchVel = _d_tension_wrt_stretchVel;
      outputs[_outIdx_dtension_dstretchVel][_icell] = dtension_dstretchVel;
      outputs[_outIdx_tension][_icell] = tension;
   }
}

void ThisModel::commitTimestep()
{
   swapOldAndNew();
}

}

