#include "GPUDiffusion.hh"
#include "DiffusionUtils.hh"
#include "SymmetricTensor.hh"
#include <vector>
#include <map>

using namespace std;

#define CONTAINS(container, item) (container.find(item) != container.end())

GPUDiffusion::GPUDiffusion(const Anatomy& anatomy, int simLoopType)
: simLoopType_(simLoopType),
  localGrid_(DiffusionUtils::findBoundingBox(anatomy, false))
{
   //move data into position
   nx_ = localGrid_.nx();
   ny_ = localGrid_.ny();
   nz_ = localGrid_.nz();
   VmBlock_.setup(vector<double>(nx_*ny_*nz_));

   nLocal_ = anatomy.nLocal();
   nRemote_ = anatomy.nRemote();
   nCells_ = anatomy.size();
   
   //initialize the block index as well.
   cellFromRed_.setup(vector<int>(anatomy.size()));
   blockFromRed_.setup(vector<int>(anatomy.size()));
   vector<int>& blockFromRedVec(blockFromRed_.modifyOnHost());
   vector<int>& cellFromRedVec(cellFromRed_.modifyOnHost());
   map<int,int> cellFromBlock;

   nRed_ = 0;
   int nBlack = 0;
   nRedLocal_ = 0;
   for (int icell=0; icell<nCells_; icell++)
   {
      Tuple ll = localGrid_.localTuple(anatomy.globalTuple(icell));
      bool isRed = ((ll.x()+ll.y()+ll.z()) % 2);
      int ired;
      if (isRed) {
         ired = nRed_++;
         if (icell < nLocal_) { nRedLocal_++; }
      }
      else
      {
         ired = nCells_-(++nBlack);
      }
      cellFromRedVec[ired] = icell;
      int iblock = ll.x()+nx_*(ll.y()+ny_*ll.z());
      blockFromRedVec[ired] = iblock;
      cellFromBlock[iblock] = icell;
   }
   assert(nRed_+nBlack == nCells_);

   const int offsets[3] = {1, nx_, nx_*ny_};
   const double areas[3] =
      {
         anatomy.dy()*anatomy.dz(),
         anatomy.dx()*anatomy.dz(),
         anatomy.dx()*anatomy.dy()
      };
   const double disc[3] = {anatomy.dx(), anatomy.dy(), anatomy.dz()};
   //initialize the array of diffusion coefficients
   cellLookup_.setup(vector<int>(anatomy.size()*3));
   vector<int>& cellLookupVec(cellLookup_.modifyOnHost());
   sigmaFaceNormal_.setup(vector<double>(anatomy.size()*9));
   vector<double>& sigmaFaceNormalVec(sigmaFaceNormal_.modifyOnHost());
   for (int ired=0; ired<nCells_; ired++)
   {
      int icell = cellFromRedVec[ired];
      int iblock = blockFromRedVec[ired];

      const SymmetricTensor& sss(anatomy.conductivity(icell));
      double thisSigma[3][3] = {{sss.a11, sss.a12, sss.a13},
                                {sss.a12, sss.a22, sss.a23},
                                {sss.a13, sss.a23, sss.a33}};
      for (int idim=0; idim<3; idim++)
      {
         int otherBlock = iblock-offsets[idim];
         int lookup = ired + idim*nCells_;
         if (CONTAINS(cellFromBlock, otherBlock))
         {
            cellLookupVec[lookup] = cellFromBlock[otherBlock];
            for (int jdim=0; jdim<3; jdim++)
            {
               int sigmaIndex = jdim +3*(idim +3*ired);
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
            cellLookupVec[lookup] = -1;
            for (int jdim=0; jdim<3; jdim++)
            {
               int sigmaIndex = jdim +3*(idim +3*ired);
               sigmaFaceNormalVec[sigmaIndex] = 0;
            }
         }
      }
   }
}

void GPUDiffusion::updateLocalVoltage(const double* VmLocal)
{
   vector<double>& VmBlockVec(VmBlock_.modifyOnDevice());
   const vector<int>& blockFromRedVec(blockFromRed_.readOnDevice());
   const vector<int>& cellFromRedVec(cellFromRed_.readOnDevice());
   
   double* VmBlockVecRaw=&VmBlockVec[0];
   const int* blockFromRedVecRaw=&blockFromRedVec[0];
   const int* cellFromRedVecRaw=&cellFromRedVec[0];
   #pragma omp target teams distribute parallel for
   for (int ired=0; ired<nRedLocal_; ired++)
   {
      int icell = cellFromRedVecRaw[ired];
      VmBlockVecRaw[blockFromRedVecRaw[ired]] = VmLocal[icell];
   }

   #pragma omp target teams distribute parallel for
   for (int ired=nRedLocal_+nRemote_; ired<nCells_; ired++)
   {
      int icell = cellFromRedVecRaw[ired];
      VmBlockVecRaw[blockFromRedVecRaw[ired]] = VmLocal[icell];
   }
   
}

void GPUDiffusion::updateRemoteVoltage(const double* VmRemote)
{
   vector<double>& VmBlockVec(VmBlock_.modifyOnDevice());
   const vector<int>& blockFromRedVec(blockFromRed_.readOnDevice());
   const vector<int>& cellFromRedVec(cellFromRed_.readOnDevice());
   
   double* VmBlockVecRaw=&VmBlockVec[0];
   const int* blockFromRedVecRaw=&blockFromRedVec[0];
   const int* cellFromRedVecRaw=&cellFromRedVec[0];   
   #pragma omp target teams distribute parallel for
   for (int ired=nRedLocal_; ired<nRedLocal_+nRemote_; ired++)
   {
      int icell = cellFromRedVecRaw[ired];
      VmBlockVecRaw[blockFromRedVecRaw[ired]] = VmRemote[icell-nLocal_];
   }
}

void GPUDiffusion::calc(VectorDouble32& dVm)
{
   double* dVmRaw=&dVm[0];
   if (simLoopType_ == 0) // it is a hard coded of enum LoopType {omp, pdr}  
   {
      #pragma omp target teams distribute parallel for
      for (int icell=0; icell<nLocal_; icell++)
      {
         dVmRaw[icell] = 0;
      }
   }
   
   const vector<double>& VmBlockVec(VmBlock_.readOnDevice());
   const vector<double>& sigmaFaceNormalVec(sigmaFaceNormal_.readOnDevice());
   const vector<int>& cellLookupVec(cellLookup_.readOnDevice());
   const vector<int>& blockFromRedVec(blockFromRed_.readOnDevice());
   const vector<int>& cellFromRedVec(cellFromRed_.readOnDevice());
   
   const double* VmBlockVecRaw=&VmBlockVec[0];
   const double* sigmaFaceNormalVecRaw=&sigmaFaceNormalVec[0];
   const int* cellLookupVecRaw=&cellLookupVec[0];
   const int* blockFromRedVecRaw=&blockFromRedVec[0];
   const int* cellFromRedVecRaw=&cellFromRedVec[0];    

   const int extents[3] = {0, nRed_, nCells_};
   const int offset[3] = {1, nx_, nx_*ny_};
  
   for (int idim=0; idim<3; idim++)
   {
      for (int red=0; red<=1; red++)
      {
         #pragma omp target teams distribute parallel for
         for (int ired=extents[red]; ired<extents[red+1]; ired++)
         { 
            int other = cellLookupVecRaw[ired + idim*nCells_];
            if (other < 0) { continue; }
         
            int thisIndex = blockFromRedVecRaw[ired];
            double flux = 0;
            for (int jdim=0; jdim<3; jdim++)
            {
               double thisGradient;
               if (idim==jdim)
               {
                  thisGradient = VmBlockVecRaw[thisIndex] - VmBlockVecRaw[thisIndex-offset[idim]];
               }
               else
               {
                  thisGradient = 0.25*(
                       VmBlockVecRaw[thisIndex +offset[jdim]              ]
                     + VmBlockVecRaw[thisIndex +offset[jdim] -offset[idim]]
                     - VmBlockVecRaw[thisIndex -offset[jdim]              ]
                     - VmBlockVecRaw[thisIndex -offset[jdim] -offset[idim]]
                  );
               }    
               double thisSigma = sigmaFaceNormalVecRaw[jdim +3*(idim + 3*ired)];
               flux += thisSigma*thisGradient;
            }
            dVmRaw[other] -= flux;
            int icell = cellFromRedVecRaw[ired];
            if (icell < nLocal_) //if this is a local cell
            {
               dVmRaw[icell] += flux;
            }
         }
      }
   }
}

unsigned* GPUDiffusion::blockIndex() {return NULL;}
double* GPUDiffusion::VmBlock() {return NULL;}
double* GPUDiffusion::dVmBlock() {return NULL;}


