#include "TT06Dev_Reaction.hh"
#include <cmath>
#include "Anatomy.hh"
#include "TT06Func.hh"
#include "pade.hh"

using namespace std;
using namespace TT06Func;


TT06Dev_Reaction::TT06Dev_Reaction(const Anatomy& anatomy,double tolerance,int mod)
: nCells_(anatomy.nLocal())
{
   static bool initialized = false;
   if (! initialized)
   {
      initialized = true;
      TT06Func::initCnst();
   }
   {
      double V0 = -100.0; 
      double V1 =  50.0; 
      double deltaV = 0.1; 
      int maxCost=128; 
      int maxTerms=64; 
      PADE **fit=TT06Func::makeFit(tolerance,V0,V1,deltaV,mod); 
      for (int i=0;fit[i]!=NULL;i++) padeCalc(fit[i],maxTerms,maxTerms,maxCost); 
      TT06Func::writeFit(fit); 
   }
   dtForFit_=0.0; 
   ttType_.resize(256, -1); 
   ttType_[30] = 0;
   ttType_[31] = 0;
   ttType_[75] = 0;
   ttType_[76] = 1;
   ttType_[77] = 2;
   ttType_[100] = 0;
   ttType_[101] = 1;
   ttType_[102] = 2;
   s_.resize(nCells_);
   for (unsigned ii=0; ii<nCells_; ++ii)
   {
      assert(anatomy.cellType(ii) >= 0 && anatomy.cellType(ii) < 256);
      int cellType = ttType_[anatomy.cellType(ii)];
      TT06Func::initState(&(s_[ii]),cellType); 
   }
   
}
TT06Dev_Reaction::~TT06Dev_Reaction()
{
}
void TT06Dev_Reaction::writeStateDev(int loop)
{
int  map[] = { dVK_i , Na_i , Ca_i , Xr1_gate , Xr2_gate , Xs_gate , m_gate , h_gate , j_gate , Ca_ss , d_gate , f_gate , f2_gate , fCass_gate , s_gate , r_gate , Ca_SR , R_prime , jL_gate};
	for (int i=0;i<nStateVar;i++) 
	{
		int k = map[i]; 
		double state = s_[0].state[k];
		printf("%d %24.14le\n",i,state); 
	}
}

void TT06Dev_Reaction::calc(double dt, const vector<double>& Vm, const vector<double>& iStim, vector<double>& dVm)
{
   TT06Func::updateNonGate(dt, nCells_,&Vm[0], &(s_[0]), &dVm[0]);
   TT06Func::updateGate(dt, nCells_,&Vm[0], &(s_[0]));
   //double c9 = TT06Func::get_c9(); 
   //for (unsigned ii=0; ii<nCells_; ++ii) s_[ii].state[TT06Func::K_i] += dt*((-dVm[ii]+iStim[ii])*c9) ;
}
void TT06Dev_Reaction::updateNonGate(double dt, const vector<double>& Vm, vector<double>& dVR)
{
    TT06Func::updateNonGate(dt, nCells_,&Vm[0], &(s_[0]), &dVR[0]);
}
void TT06Dev_Reaction::updateGate(double dt, const vector<double>& Vm)
{
   TT06Func::updateGate(dt, nCells_,&Vm[0], &(s_[0]));
}


void TT06Dev_Reaction::initializeMembraneVoltage(std::vector<double>& Vm)
{
   assert(Vm.size() >= s_.size());
   for (unsigned ii=0; ii<s_.size(); ++ii)
      Vm[ii] = TT06Func::defaultVoltage(s_[ii].cellType);
}
