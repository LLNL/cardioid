#include "Anatomy.hh"
#include "ConstantReaction.hh"
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

ConstantReaction::ConstantReaction(const Anatomy& anatomy,
                 const std::vector<double>& eta,
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
                            ro_larray_ptr<double> Vm,
                            ro_larray_ptr<double> iStim,
                            wo_larray_ptr<double> dVm)
{
   assert( cellModel_!=0 );
   ContextRegion region(CPU);
   Vm.use();
   iStim.use();
   dVm.use();
   
   //double sum=0.;
   for (unsigned ii=0; ii<anatomy_.nLocal(); ++ii)
   {
      Tuple globalTuple = anatomy_.globalTuple(ii);
      const double x=((double)globalTuple.x()+0.5)*anatomy_.dx();
      const double y=((double)globalTuple.y()+0.5)*anatomy_.dy();
      const double z=((double)globalTuple.z()+0.5)*anatomy_.dz();
      double update = cellModel_->calc(x,y,z);
      assert(update ==update); //test for nan
      dVm[ii] = update;
   }
}

void ConstantReaction::initializeMembraneVoltage(wo_larray_ptr<double> Vm)
{
   ContextRegion region(CPU);
   Vm.use();
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

int ConstantReaction::getVarHandle(const std::string& varName) const
{
   return -1;
}
