#include "OpenmpGpuRedblackDiffusion.hh"
#include "DiffusionUtils.hh"
#include "SymmetricTensor.hh"
#include <vector>
#include <map>

using namespace std;

#define CONTAINS(container, item) (container.find(item) != container.end())

/* Here's the design constraints for this piece of code:
 *
 * - When we do the hot loop, we're updating adjacent cells at the
 *   same time.  Therefore, we can never be doing operations on the
 *   adjacent cells without creating a race condition.
 *
 * - We'd like the local and the remote coordinates to be contiguous
 *   as much as possible.
 *
 * Here, I have 3 orderings:
 * cell ordering - the ordering the external arrays (where we have no control)
 * red/black ordering - redLocal redRemote blackRemote blackLocal
 * block ordering - ordering of the points in the block
 *
 * In the constructor I'm constructing the mappings between these
 * orderings. I've set up the red ordering so that all the reds are
 * first, followed by all the blacks.  The remotes are in a contiguous
 * block in the middle.
 * 
 * The hot loop runs on the red black ordering.  
 *
 */

OpenmpGpuRedblackDiffusion::OpenmpGpuRedblackDiffusion(const Anatomy& anatomy, int simLoopType)
: simLoopType_(simLoopType),
  localGrid_(DiffusionUtils::findBoundingBox(anatomy, false))
{
   //move data into position
   nx_ = localGrid_.nx();
   ny_ = localGrid_.ny();
   nz_ = localGrid_.nz();
   VmBlock_.setup(PinnedVector<double>(nx_*ny_*nz_));

   nLocal_ = anatomy.nLocal();
   nRemote_ = anatomy.nRemote();
   nCells_ = anatomy.size();
   
   //initialize the block index as well.
   cellFromRed_.setup(PinnedVector<int>(anatomy.size()));
   blockFromRed_.setup(PinnedVector<int>(anatomy.size()));
   ArrayView<int> blockFromRedVec = blockFromRed_;
   ArrayView<int> cellFromRedVec = cellFromRed_;
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
   cellLookup_.setup(PinnedVector<int>(anatomy.size()*3));
   ArrayView<int> cellLookupVec = cellLookup_;
   sigmaFaceNormal_.setup(PinnedVector<double>(anatomy.size()*9));
   ArrayView<double> sigmaFaceNormalVec = sigmaFaceNormal_;
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

void OpenmpGpuRedblackDiffusion::updateLocalVoltage(const Managed<ArrayView<double>> VmLocal_managed)
{
   ConstArrayView<double> VmLocal = VmLocal_managed;
   ArrayView<double> VmBlockVec = VmBlock_; 
   ConstArrayView<int> blockFromRedVec = blockFromRed_;
   ConstArrayView<int> cellFromRedVec = cellFromRed_;
   
   double* VmBlockVecRaw=&VmBlockVec[0];
   const int* blockFromRedVecRaw=&blockFromRedVec[0];
   const int* cellFromRedVecRaw=&cellFromRedVec[0];
   //#pragma omp target teams distribute parallel for
   for (int ired=0; ired<nRedLocal_; ired++)
   {
      int icell = cellFromRedVecRaw[ired];
      VmBlockVecRaw[blockFromRedVecRaw[ired]] = VmLocal[icell];
   }

   //#pragma omp target teams distribute parallel for
   for (int ired=nRedLocal_+nRemote_; ired<nCells_; ired++)
   {
      int icell = cellFromRedVecRaw[ired];
      VmBlockVecRaw[blockFromRedVecRaw[ired]] = VmLocal[icell];
   }
   
}

void OpenmpGpuRedblackDiffusion::updateRemoteVoltage(const Managed<ArrayView<double>> VmRemote_managed)
{
   ConstArrayView<double> VmRemote = VmRemote_managed;
   ArrayView<double> VmBlockVec = VmBlock_; 
   ConstArrayView<int> blockFromRedVec = blockFromRed_;
   ConstArrayView<int> cellFromRedVec = cellFromRed_;
   
   double* VmBlockVecRaw=&VmBlockVec[0];
   const int* blockFromRedVecRaw=&blockFromRedVec[0];
   const int* cellFromRedVecRaw=&cellFromRedVec[0];   
   //#pragma omp target teams distribute parallel for
   for (int ired=nRedLocal_; ired<nRedLocal_+nRemote_; ired++)
   {
      int icell = cellFromRedVecRaw[ired];
      VmBlockVecRaw[blockFromRedVecRaw[ired]] = VmRemote[icell-nLocal_];
   }
}

void OpenmpGpuRedblackDiffusion::calc(Managed<ArrayView<double>> dVm_managed)
{
   actualCalc(*this, dVm_managed);
}

void actualCalc(OpenmpGpuRedblackDiffusion& self, ArrayView<double> dVm)
{
   int self_nLocal_ = self.nLocal_;
   int self_nCells_ = self.nCells_;
   int self_nRed_ = self.nRed_;
   int self_nx_ = self.nx_;
   int self_ny_ = self.ny_;
   double* dVmRaw=&dVm[0];
   if (self.simLoopType_ == 0) // it is a hard coded of enum LoopType {omp, pdr}  
   {
      //#pragma omp target teams distribute parallel for firstprivate(self_nLocal_)
      for (int icell=0; icell<self_nLocal_; icell++)
      {
         dVmRaw[icell] = 0;
      }
   }
   
   ConstArrayView<double> VmBlockVec = self.VmBlock_;
   ConstArrayView<double> sigmaFaceNormalVec = self.sigmaFaceNormal_;
   ConstArrayView<int> cellLookupVec = self.cellLookup_;
   ConstArrayView<int> blockFromRedVec = self.blockFromRed_;
   ConstArrayView<int> cellFromRedVec = self.cellFromRed_;
   
   const double* VmBlockVecRaw=&VmBlockVec[0];
   const double* sigmaFaceNormalVecRaw=&sigmaFaceNormalVec[0];
   const int* cellLookupVecRaw=&cellLookupVec[0];
   const int* blockFromRedVecRaw=&blockFromRedVec[0];
   const int* cellFromRedVecRaw=&cellFromRedVec[0];    

   for (int idim=0; idim<3; idim++)
   {
      for (int red=0; red<=1; red++)
      {
         const int extents[3] = {0, self_nRed_, self_nCells_};
         int lower = extents[red];
         int upper = extents[red+1];

         #pragma target teams omp distribute parallel for firstprivate(idim, red, self_nCells_, self_nx_, self_ny_, self_nLocal_, lower, upper)
         for (int ired=lower; ired<upper; ired++)
         { 
            const int offset[3] = {1, self_nx_, self_nx_*self_ny_};
            int other = cellLookupVecRaw[ired + idim*self_nCells_];
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
            if (icell < self_nLocal_) //if this is a local cell
            {
               dVmRaw[icell] += flux;
            }
         }
      }
   }
}

unsigned* OpenmpGpuRedblackDiffusion::blockIndex() {return NULL;}
double* OpenmpGpuRedblackDiffusion::VmBlock() {return NULL;}
double* OpenmpGpuRedblackDiffusion::dVmBlock() {return NULL;}


