#include "CUDADiffusion.hh"
#include "DiffusionUtils.hh"
#include "SymmetricTensor.hh"
#include <vector>
#include <map>
#include <iostream>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime_api.h>

//update is called before calc
//dVm is combined with Vm from reaction without any rearrangement
//
//Hence
//'update' must rearrange the dVmBlock
//activate map_dVm
//
#define CUDA_VERIFY(x) do { cudaError_t error = x; if (error != cudaSuccess) { cout << error << endl; assert(error == cudaSuccess && #x ); } } while(0)

using namespace std;

#define CONTAINS(container, item) (container.find(item) != container.end())

typedef double Real;

void call_cuda_kernels(ro_mgarray_ptr<Real> VmRaw, rw_mgarray_ptr<Real> dVmRaw, ro_mgarray_ptr<Real> sigmaRaw, int nx, int ny, int nz, wo_mgarray_ptr<Real> dVmOut, ro_mgarray_ptr<int> lookup,int nCells);

void copy_to_block(wo_mgarray_ptr<double> blockCPU, ro_mgarray_ptr<int> lookupCPU, ro_mgarray_ptr<double> sourceCPU, const int begin, const int end);

CUDADiffusion::CUDADiffusion(const Anatomy& anatomy, int simLoopType, double newDiffusionScale)
: anatomy_(anatomy), simLoopType_(simLoopType), Diffusion(newDiffusionScale),
  localGrid_(DiffusionUtils::findBoundingBox(anatomy, false))
{
   //move data into position
   nx_ = localGrid_.nx();
   ny_ = localGrid_.ny();
   nz_ = localGrid_.nz();
   VmBlock_.resize(nx_*ny_*nz_);
   dVmBlock_.resize(nx_*ny_*nz_);
   cellLookup_.resize(anatomy.size());

   nLocal_ = anatomy.nLocal();
   nRemote_ = anatomy.nRemote();
   nCells_ = anatomy.size();
   
   //initialize the block index as well.
   //cellFromRed_.setup(PinnedVector<int>(anatomy.size()));
   //blockFromRed_.setup(PinnedVector<int>(anatomy.size()));
   //vector<int>& blockFromRedVec(blockFromRed_.modifyOnHost());
   //vector<int>& cellFromRedVec(cellFromRed_.modifyOnHost());
   map<int,int> cellFromBlock;
   wo_array_ptr<int> cellLookupVec = cellLookup_.writeonly(CPU);
   for (int icell=0; icell<nCells_; icell++)
   {
      //index to coordinate 
      Tuple ll = localGrid_.localTuple(anatomy.globalTuple(icell));
      int iblock = ll.z() + nz_ * ( ll.y() + ny_ * ll.x() );
      cellFromBlock[iblock] = icell;
      cellLookupVec[icell] = iblock;
   }

   //auto cellLookupAccess=cellLookup_.writeonly();
   //copy(cellLookupVec.begin(), cellLookupVec.end(), cellLookupAccess.begin());

   const int offsets[3] = {ny_*nz_, nz_ , 1};
   const double areas[3] =
      {
         anatomy.dy()*anatomy.dz(),
         anatomy.dx()*anatomy.dz(),
         anatomy.dx()*anatomy.dy()
      };
   const double disc[3] = {anatomy.dx(), anatomy.dy(), anatomy.dz()};
   const double gridCellVolume = disc[0]*disc[1]*disc[2];
   //initialize the array of diffusion coefficients
   sigmaFaceNormal_.resize(nx_*ny_*nz_*9);
   wo_array_ptr<double> sigmaFaceNormalVec = sigmaFaceNormal_.writeonly(CPU);

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
               double sigmaValue = thisSigma[idim][jdim]*areas[idim]/disc[jdim]*diffusionScale()/gridCellVolume;
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
   {
      wo_array_ptr<double> VmBlock = VmBlock_.writeonly(GPU);
      CUDA_VERIFY(cudaMemset(VmBlock.raw(), 0, VmBlock.size()*sizeof(double)));
   }
}

void CUDADiffusion::updateLocalVoltage(ro_mgarray_ptr<double> VmLocal)
{
   copy_to_block(VmBlock_, cellLookup_, VmLocal, 0, nLocal_);
}

void CUDADiffusion::updateRemoteVoltage(ro_mgarray_ptr<double> VmRemote)
{
   copy_to_block(VmBlock_, cellLookup_, VmRemote, nLocal_, nCells_);
}

void CUDADiffusion::calc(rw_mgarray_ptr<double> dVm)
{
   if (simLoopType_ == 0) // it is a hard coded of enum LoopType {omp, pdr}  
   {
      wo_array_ptr<double> dVm_to_zero = dVm.useOn(GPU);
      CUDA_VERIFY(cudaMemset(dVm_to_zero.raw(), 0, dVm_to_zero.size()*sizeof(double)));
   }

   call_cuda_kernels(VmBlock_,dVmBlock_,sigmaFaceNormal_,nx_,ny_,nz_,
                     dVm,cellLookup_,nLocal_);
}

unsigned* CUDADiffusion::blockIndex() {return NULL;}
double* CUDADiffusion::VmBlock() {return NULL;}
double* CUDADiffusion::dVmBlock() {return NULL;}


