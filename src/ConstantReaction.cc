#include "Anatomy.hh"
#include "ConstantReaction.hh"
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

ConstantReaction::ConstantReaction(const Anatomy& anatomy,
                 const vector<double>& eta,
                 const SymmetricTensor& sigma1,
                 const SymmetricTensor& sigma2,
                 const SymmetricTensor& sigma3,
                 const double alpha,
                 const double beta,
                 const double gamma,
                 const int printRate)
: anatomy_(anatomy)
{
   cellModel_=new ConstantModel(eta,
                 sigma1, sigma2, sigma3,
                 alpha, beta, gamma);

   printRate_=printRate;
}

void ConstantReaction::calc(double dt,
                            const Managed<ArrayView<double>> Vm_,
                            const Managed<ArrayView<double>> iStim_,
                            Managed<ArrayView<double>> dVm_)
{
   assert( cellModel_!=0 );

   ArrayView<double> dVm(dVm_);
   //double sum=0.;
   for (unsigned ii=0; ii<anatomy_.nLocal(); ++ii)
   {
      Tuple globalTuple = anatomy_.globalTuple(ii);
      const double x=((double)globalTuple.x()+0.5)*anatomy_.dx();
      const double y=((double)globalTuple.y()+0.5)*anatomy_.dy();
      const double z=((double)globalTuple.z()+0.5)*anatomy_.dz();
      dVm[ii] = cellModel_->calc(x,y,z);
       
      assert( dVm[ii]==dVm[ii] ); // test for nan
   }
}

void ConstantReaction::initializeMembraneVoltage(ArrayView<double> Vm)
{
#if 0
   for (unsigned ii=0; ii<Vm.size(); ++ii)
      Vm[ii] = 0.;
#else
   // initialize with exact solution
   const double alpha=cellModel_->getAlpha();
   const double beta =cellModel_->getBeta();
   const double gamma=cellModel_->getGamma();
   for (unsigned ii=0; ii<anatomy_.nLocal(); ++ii)
   {
      Tuple globalTuple = anatomy_.globalTuple(ii);
      const double x=((double)globalTuple.x()+0.5)*anatomy_.dx();
      const double y=((double)globalTuple.y()+0.5)*anatomy_.dy();
      const double z=((double)globalTuple.z()+0.5)*anatomy_.dz();
      Vm[ii] = alpha*x+beta*y+gamma*z;
   }
#endif
}

// simple linear function (just for testing sensor!)
double ConstantReaction::getValue(int iCell, int handle) const
{
   assert( handle==-1 );
   
   const double alpha=cellModel_->getAlpha();
   const double beta =cellModel_->getBeta();
   const double gamma=cellModel_->getGamma();

   Tuple globalTuple = anatomy_.globalTuple(iCell);
   const double x=((double)globalTuple.x()+0.5)*anatomy_.dx();
   const double y=((double)globalTuple.y()+0.5)*anatomy_.dy();
   const double z=((double)globalTuple.z()+0.5)*anatomy_.dz();

   return alpha*x+beta*y+gamma*z;
}

int ConstantReaction::getVarHandle(const string& varName) const
{
   return -1;
}
