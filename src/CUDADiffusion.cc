#include "CUDADiffusion.hh"
#include "DiffusionUtils.hh"
#include "SymmetricTensor.hh"
#include <vector>
#include <map>
#include <iostream>
#include <stdio.h>

//update is called before calc
//dVm is combined with Vm from reaction without any rearrangement
//
//Hence
//'update' must rearrange the dVmBlock
//activate map_dVm

using namespace std;

#define CONTAINS(container, item) (container.find(item) != container.end())

typedef double Real;

extern "C"
void call_cuda_kernels(const Real *VmRaw, Real *dVmRaw, const Real *sigmaRaw, int nx, int ny, int nz, Real *dVmOut, const int *lookup,int nCells);

void copy_to_block(double* blockCPU, const int* lookupCPU, const double* sourceCPU, const int begin, const int end);

CUDADiffusion::CUDADiffusion(const Anatomy& anatomy, int simLoopType)
: anatomy_(anatomy), simLoopType_(simLoopType),
  localGrid_(DiffusionUtils::findBoundingBox(anatomy, false))
{
   //move data into position
   nx_ = localGrid_.nx();
   ny_ = localGrid_.ny();
   nz_ = localGrid_.nz();
   VmBlock_.setup(vector<double>(nx_*ny_*nz_));
   dVmBlock_.setup(vector<double>(nx_*ny_*nz_));
   cellLookup_.setup(vector<int>(anatomy.size()));

   nLocal_ = anatomy.nLocal();
   nRemote_ = anatomy.nRemote();
   nCells_ = anatomy.size();
   
   //initialize the block index as well.
   //cellFromRed_.setup(vector<int>(anatomy.size()));
   //blockFromRed_.setup(vector<int>(anatomy.size()));
   //vector<int>& blockFromRedVec(blockFromRed_.modifyOnHost());
   //vector<int>& cellFromRedVec(cellFromRed_.modifyOnHost());
   map<int,int> cellFromBlock;
   vector<int>& cellLookupVec(cellLookup_.modifyOnHost());
//   vector<int>& blockLookupVec(blockLookup_.modifyOnHost());

   for (int icell=0; icell<nCells_; icell++)
   {
      //index to coordinate 
      Tuple ll = localGrid_.localTuple(anatomy.globalTuple(icell));
      int iblock = ll.z() + nz_ * ( ll.y() + ny_ * ll.x() );
      cellFromBlock[iblock] = icell;
      cellLookupVec[icell] = iblock;
   }

   const int offsets[3] = {ny_*nz_, nz_ , 1};
   const double areas[3] =
      {
         anatomy.dy()*anatomy.dz(),
         anatomy.dx()*anatomy.dz(),
         anatomy.dx()*anatomy.dy()
      };
   const double disc[3] = {anatomy.dx(), anatomy.dy(), anatomy.dz()};
   //initialize the array of diffusion coefficients
   sigmaFaceNormal_.setup(vector<double>(nx_*ny_*nz_*9));
   vector<double>& sigmaFaceNormalVec(sigmaFaceNormal_.modifyOnHost());

   for (int icell=0; icell<nCells_; icell++)
   {
      //int icell = cellFromRedVec[ired];
      //int iblock = blockFromRedVec[ired];
      Tuple ll = localGrid_.localTuple(anatomy.globalTuple(icell));
      int xx[3] ={ ll.x() , ll.y() , ll.z() };
      int iblock = ll.z()+nz_*(ll.y()+ny_*ll.x());

      const SymmetricTensor& sss(anatomy.conductivity(icell));
      double thisSigma[3][3] = {{sss.a11, sss.a12, sss.a13},
                                {sss.a12, sss.a22, sss.a23},
                                {sss.a13, sss.a23, sss.a33}};

      for (int idim=0; idim<3; idim++)
      {
         int otherBlock = iblock-offsets[idim];
         if (CONTAINS(cellFromBlock, otherBlock)) //if my left cell block exists
         {
            for (int jdim=0; jdim<3; jdim++)
            {
               int sigmaIndex = xx[2] + nz_ * ( xx[1] + ny_ * ( xx[0] + nx_ * ( jdim + 3 * idim)));
               double sigmaValue = thisSigma[idim][jdim]*areas[idim]/disc[jdim];
               if (idim != jdim)
               {
                  bool canComputeGradient=
                        CONTAINS(cellFromBlock, iblock              -offsets[jdim])
                     && CONTAINS(cellFromBlock, iblock              +offsets[jdim])
                     && CONTAINS(cellFromBlock, iblock-offsets[idim]-offsets[jdim])
                     && CONTAINS(cellFromBlock, iblock-offsets[idim]+offsets[jdim])
                     ;
                  if (! canComputeGradient)
                  {
                     sigmaValue = 0;
                  }
               }
               sigmaFaceNormalVec[sigmaIndex] = sigmaValue;
            }
         }
         else
         {
            //cellLookupVec[lookup] = -1;
            for (int jdim=0; jdim<3; jdim++)
            {
               int sigmaIndex = xx[2] + nz_ * ( xx[1] + ny_ * ( xx[0] + nx_ * ( jdim + 3 * idim)));
               sigmaFaceNormalVec[sigmaIndex] = 0;
            }
         }
      }
   }
}

void CUDADiffusion::updateLocalVoltage(const double* VmLocal)
{
   vector<double>& VmBlockVec(VmBlock_.modifyOnDevice());
   double* VmBlockVecRaw=&VmBlockVec[0];
   const vector<int>& cellLookupVec(cellLookup_.readOnDevice());
   const int* lookupRaw=&cellLookupVec[0];
   copy_to_block(VmBlockVecRaw, lookupRaw, VmLocal, 0, nLocal_);
}

void CUDADiffusion::updateRemoteVoltage(const double* VmRemote)
{
   vector<double>& VmBlockVec(VmBlock_.modifyOnDevice());
   double* VmBlockVecRaw=&VmBlockVec[0];
   const vector<int>& cellLookupVec(cellLookup_.readOnDevice());
   const int* lookupRaw=&cellLookupVec[0];
   copy_to_block(VmBlockVecRaw, lookupRaw, VmRemote, nLocal_, nCells_);
}

void CUDADiffusion::calc(VectorDouble32& dVm){
   actualCalc(*this, dVm);
}

void actualCalc(CUDADiffusion& self, VectorDouble32& dVm)
{
   int self_nLocal_ = self.nLocal_;
   int self_nCells_ = self.nCells_;
   int nx = self.nx_;
   int ny = self.ny_;
   int nz = self.nz_;
   double* dVmOut=&dVm[0];
   if (self.simLoopType_ == 0) // it is a hard coded of enum LoopType {omp, pdr}  
   {
      #pragma omp target teams distribute parallel for firstprivate(self_nLocal_)
      for (int icell=0; icell<self_nLocal_; icell++)
      {
         dVmOut[icell] = 0;
      }
   }
   
   vector<double>& dVmBlockVec(self.dVmBlock_.modifyOnDevice());
   const vector<double>& VmBlockVec(self.VmBlock_.readOnDevice());
   const vector<double>& sigmaFaceNormalVec(self.sigmaFaceNormal_.readOnDevice());
   const vector<int>& cellLookupVec(self.cellLookup_.readOnDevice());
   
   const double* VmRaw=&VmBlockVec[0];
   double* dVmRaw=&dVmBlockVec[0];
   const double* sigmaRaw=&sigmaFaceNormalVec[0];
   const int*    lookupRaw=&cellLookupVec[0];
   double* temp=dVmOut;
   
   {
      call_cuda_kernels(VmRaw,
                        dVmRaw,
                        sigmaRaw,
                        nx,ny,nz,
                        dVmOut,
                        lookupRaw,
                        self_nCells_);
   }

}

unsigned* CUDADiffusion::blockIndex() {return NULL;}
double* CUDADiffusion::VmBlock() {return NULL;}
double* CUDADiffusion::dVmBlock() {return NULL;}


