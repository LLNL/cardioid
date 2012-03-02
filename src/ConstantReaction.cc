#include "Anatomy.hh"
#include "ConstantReaction.hh"
#include <mpi.h>
#include <iostream>
#include <iomanip>

using namespace std;

ConstantReaction::ConstantReaction(const Anatomy& anatomy,
                 const double eta[3],
                 const SymmetricTensor& sigma1,
                 const SymmetricTensor& sigma2,
                 const SymmetricTensor& sigma3,
                 const double alpha,
                 const double beta,
                 const double gamma)
: anatomy_(anatomy)
{
   cellModel_=new ConstantModel(eta,
                 sigma1, sigma2, sigma3,
                 alpha, beta, gamma);
}

void ConstantReaction::calc(double dt,
                         const vector<double>& Vm,
                         const vector<double>& iStim,
                         vector<double>& dVm)
{
   assert( cellModel_!=0 );
   
   //double sum=0.;
   for (unsigned ii=0; ii<anatomy_.nLocal(); ++ii)
   {
      Tuple globalTuple = anatomy_.globalTuple(ii);
      const double x=((double)globalTuple.x()+0.5)*anatomy_.dx();
      const double y=((double)globalTuple.y()+0.5)*anatomy_.dy();
      const double z=((double)globalTuple.z()+0.5)*anatomy_.dz();
      dVm[ii] = cellModel_->calc(x,y,z);
      //sum+=dVm[ii];
      
      assert( dVm[ii]==dVm[ii] ); // test for nan
   }

   //double gsum=0.;
   //MPI_Allreduce(&sum, &gsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);   

   //for (unsigned ii=0; ii<anatomy_.nLocal(); ++ii)
   //{
   //   dVm[ii]-=gsum;
   //}
   
}

void ConstantReaction::initializeMembraneVoltage(std::vector<double>& Vm)
{
#if 0
   for (unsigned ii=0; ii<Vm.size(); ++ii)
      Vm[ii] = 0.;
#else
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

