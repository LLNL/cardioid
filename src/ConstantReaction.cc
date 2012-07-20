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
                         const VectorDouble32& Vm,
                         const vector<double>& iStim,
                         VectorDouble32& dVm)
{
   assert( cellModel_!=0 );
   
   static int count = 0;
   
   if( printRate_>0 )
   if( count % printRate_ == 0)
      compareWithExactSol(Vm);
   count++;

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

void ConstantReaction::compareWithExactSol(const VectorDouble32& Vm)const
{
   const double alpha=cellModel_->getAlpha();
   const double beta =cellModel_->getBeta();
   const double gamma=cellModel_->getGamma();

   double avg1=0.;
   double avg2=0.;
   for (unsigned ii=0; ii<anatomy_.nLocal(); ++ii)
   {
      Tuple globalTuple = anatomy_.globalTuple(ii);
      const double x=((double)globalTuple.x()+0.5)*anatomy_.dx();
      const double y=((double)globalTuple.y()+0.5)*anatomy_.dy();
      const double z=((double)globalTuple.z()+0.5)*anatomy_.dz();
      avg1+=alpha*x+beta*y+gamma*z;
      avg2+=Vm[ii];
   }
   double tmp=0.;
   MPI_Allreduce(&avg1,&tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   avg1=tmp/anatomy_.nGlobal();
   MPI_Allreduce(&avg2,&tmp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   avg2=tmp/anatomy_.nGlobal();


   double diff2=0.;
   for (unsigned ii=0; ii<anatomy_.nLocal(); ++ii)
   {
      Tuple globalTuple = anatomy_.globalTuple(ii);
      const double x=((double)globalTuple.x()+0.5)*anatomy_.dx();
      const double y=((double)globalTuple.y()+0.5)*anatomy_.dy();
      const double z=((double)globalTuple.z()+0.5)*anatomy_.dz();
      const double diff = Vm[ii]-avg2 - (alpha*x+beta*y+gamma*z-avg1);
      diff2+=diff*diff;
   }
   double gdiff2=0.;
   MPI_Allreduce(&diff2,&gdiff2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
   gdiff2=sqrt(gdiff2/anatomy_.nGlobal());
   if( myRank==0 )
      cout<<"Diff. with exact solution = "<<gdiff2<<endl;
}

void ConstantReaction::initializeMembraneVoltage(VectorDouble32& Vm)
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
