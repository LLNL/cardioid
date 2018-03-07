#include "OpenmpGpuFlatDiffusion.hh"
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

OpenmpGpuFlatDiffusion::OpenmpGpuFlatDiffusion(const Anatomy& anatomy, int simLoopType)
: simLoopType_(simLoopType),
  localGrid_(DiffusionUtils::findBoundingBox(anatomy, false))
{
   //move data into position
   nx_ = localGrid_.nx();
   ny_ = localGrid_.ny();
   nz_ = localGrid_.nz();
   VmBlock_.setup(PinnedVector<double>(nx_*ny_*nz_));
   dVmBlock_.setup(PinnedVector<double>(nx_*ny_*nz_));
   skip_reorder_ = false; // for now always set to false
   
   nLocal_ = anatomy.nLocal();
   nRemote_ = anatomy.nRemote();
   nCells_ = anatomy.size();
   
   //initialize the block index as well.
   blockFromCell_.setup(PinnedVector<int>(anatomy.size()));   
   ArrayView<int> blockFromCellVec = blockFromCell_;
   map<int,int> cellFromBlock;

   for (int icell=0; icell<nCells_; icell++)
   {
      Tuple ll = localGrid_.localTuple(anatomy.globalTuple(icell));
      int iblock = ll.x()+nx_*(ll.y()+ny_*ll.z());

      blockFromCellVec[icell] = iblock;
      //FIXME
      //check to see if the ordering is the same as the block order,
      //if this is true, then we can avoid a bunch of indexing in the diffusion compute
      skip_reorder_= false;

      cellFromBlock[iblock] = icell;
   }

   const int offsets[3] = {1, nx_, nx_*ny_};
   const double areas[3] =
      {
         anatomy.dy()*anatomy.dz(),
         anatomy.dx()*anatomy.dz(),
         anatomy.dx()*anatomy.dy()
      };
   const double disc[3] = {anatomy.dx(), anatomy.dy(), anatomy.dz()};

   //initialize the array of diffusion coefficients
   int sigmaNx=nx_-1;
   int sigmaNy=ny_-1;
   int sigmaNz=nz_-1;
   sigmaFaceNormal_.setup(PinnedVector<double>(sigmaNx*sigmaNy*sigmaNz*9));
   ArrayView<double> sigmaFaceNormalVec = sigmaFaceNormal_;
   for (int iz=0; iz<nz_-1; iz++)
   {
      for (int iy=0; iy<ny_-1; iy++)
      {
         for (int ix=0; ix<nx_-1; ix++)
         {
            int iblock = ix+1 + nx_*(iy+1 + ny_*(iz+1));
            int isigma = ix + sigmaNx*(iy + sigmaNy*iz);

            if (CONTAINS(cellFromBlock, iblock))
            {
               int icell = cellFromBlock[iblock];
               const SymmetricTensor& sss(anatomy.conductivity(icell));
               double thisSigma[3][3] = {{sss.a11, sss.a12, sss.a13},
                                         {sss.a12, sss.a22, sss.a23},
                                         {sss.a13, sss.a23, sss.a33}};
               for (int idim=0; idim<3; idim++)
               {
                  int otherBlock = iblock-offsets[idim];
                  if (CONTAINS(cellFromBlock, otherBlock))
                  {
                     for (int jdim=0; jdim<3; jdim++)
                     {
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
                        sigmaFaceNormalVec[jdim +3*(idim + 3*isigma)] = sigmaValue;
                     }
                  }
                  else
                  {
                     for (int jdim=0; jdim<3; jdim++)
                     {
                        sigmaFaceNormalVec[jdim +3*(idim + 3*isigma)] = 0;
                     }
                  }
               }
            }
            else
            {
               for (int idim=0; idim<3; idim++)
               {
                  for (int jdim=0; jdim<3; jdim++)
                  {
                     sigmaFaceNormalVec[jdim +3*(idim + 3*isigma)] = 0;
                  }
               }
            }
         }
      }
   }
}

void OpenmpGpuFlatDiffusion::updateLocalVoltage(const Managed<ArrayView<double>> VmLocal_managed)
{
   if (skip_reorder_)
   {
      //FIXME
   }
   else
   {
      ConstArrayView<double> VmLocal = VmLocal_managed;
      ArrayView<double> VmBlockVec = VmBlock_;
      double* VmBlockRaw=&VmBlockVec[0];

      ConstArrayView<int> blockFromCellVec = blockFromCell_;
      const int* blockFromCellRaw=&blockFromCellVec[0];

      //#pragma omp target teams distribute parallel for
      for (int icell=0; icell<nLocal_; icell++)
      {
         VmBlockRaw[blockFromCellRaw[icell]] = VmLocal[icell];
      }
   }
}

void OpenmpGpuFlatDiffusion::updateRemoteVoltage(const Managed<ArrayView<double>> VmRemote_managed)
{
   if (skip_reorder_)
   {
      //FIXME
   }
   else
   {
      ConstArrayView<double> VmRemote = VmRemote_managed;
      ArrayView<double> VmBlockVec = VmBlock_;
      double* VmBlockRaw=&VmBlockVec[0];
      
      ConstArrayView<int> blockFromCellVec = blockFromCell_;
      const int* blockFromCellRaw=&blockFromCellVec[0];
      
      //#pragma omp target teams distribute parallel for
      for (int icell=nLocal_; icell<nCells_; icell++)
      {
         VmBlockRaw[blockFromCellRaw[icell]] = VmRemote[icell-nLocal_];
      }
   }
}

void OpenmpGpuFlatDiffusion::calc(Managed<ArrayView<double>> dVm_managed)
{
   actualCalc(*this, dVm_managed);
}

void actualCalc(OpenmpGpuFlatDiffusion& self, ArrayView<double> dVm)
{
   int self_nx_ = self.nx_;
   int self_ny_ = self.ny_;
   int self_nz_ = self.nz_;
   int self_nLocal_ = self.nLocal_;
   
   double* dVmRaw=&dVm[0];
   ArrayView<double> dVmBlockVec = self.dVmBlock_;
   double* dVmBlockRaw=&dVmBlockVec[0];

   if (self.skip_reorder_)
   {
      //FIXME
   }
   else
   {
      if (self.simLoopType_ == 0) // it is a hard coded of enum LoopType {omp, pdr}  
      {
         //#pragma omp target teams distribute parallel for firstprivate(self_nLocal_)
         for (int icell=0; icell<self_nLocal_; icell++)
         {
            dVmRaw[icell] = 0;
         }
      }
 

      //#pragma omp target teams distribute parallel for firstprivate(self_nx_,self_ny_,self_nz_)
      for (int ii=0; ii<(self_nx_-2)*(self_ny_-2)*(self_nz_-2); ii++)
      {
         int iz = ii / ((self_nx_-2)*(self_ny_-2));
         int rz = ii % ((self_nx_-2)*(self_ny_-2));
         int iy = rz / (self_nx_-2);
         int ix = rz % (self_nx_-2);
         
         ix++;
         iy++;
         iz++;
         
         dVmBlockRaw[ix + self_nx_*(iy + self_ny_*iz)] = 0;
      }
   }
   
   ConstArrayView<double> VmBlockVec(self.VmBlock_.readOnDevice());
   const double* VmBlockRaw=&VmBlockVec[0];

   ConstArrayView<double> sigmaFaceNormalVec(self.sigmaFaceNormal_.readOnDevice());
   const double* sigmaFaceNormalRaw=&sigmaFaceNormalVec[0];

   int iterateMax = (self.nx_-1)*(self.ny_-1)*(self.nz_-1);
   
   #pragma unroll(3)
   for (int idim=0; idim<3; idim++)
   {
      #pragma unroll(2)
      for (int red=0; red<=1; red++)
      {
         //#pragma omp target teams distribute parallel for firstprivate(self_nx_, self_ny_, iterateMax)
         for (int ired=red; ired<iterateMax; ired+=2)
         {
            const int offset[3] = {1, self_nx_, self_nx_*self_ny_};

            int sigmaNx=self_nx_-1;
            int sigmaNy=self_ny_-1;
      
            int iz = ired / (sigmaNy*sigmaNx);
            int rzo = ired % (sigmaNy*sigmaNx);
            
            int iy_guess = rzo / sigmaNx;

            int rz = rzo ^
               (
                  (
                     (!sigmaNx & (iy_guess^iz))
                     |
                     ( sigmaNx & !(sigmaNy) & iz)
                  )
                  %2
               );
            int iy = rz / sigmaNx;
            int ix = rz % sigmaNx;

            int iblock = ix+1 + self_nx_*(iy+1 + self_ny_*(iz+1));
            int isigma = ix + sigmaNx*(iy + sigmaNy*iz);
            double flux = 0;
            for (int jdim=0; jdim<3; jdim++)
            {
               double thisGradient;
               if (idim==jdim)
               {
                  thisGradient = VmBlockRaw[iblock] - VmBlockRaw[iblock-offset[idim]];
               }
               else
               {
                  thisGradient = 0.25*(
                       VmBlockRaw[iblock +offset[jdim]              ]
                     + VmBlockRaw[iblock +offset[jdim] -offset[idim]]
                     - VmBlockRaw[iblock -offset[jdim]              ]
                     - VmBlockRaw[iblock -offset[jdim] -offset[idim]]
                  );
               }    
               double thisSigma = sigmaFaceNormalRaw[jdim +3*(idim + 3*isigma)];
               flux += thisSigma*thisGradient;
            }
            dVmBlockRaw[iblock] -= flux;
            dVmBlockRaw[iblock-offset[idim]] += flux;
         }
      }
   }

   if (self.skip_reorder_)
   {
      //FIXME
   }
   else
   {
      ConstArrayView<int> blockFromCellVec = self.blockFromCell_;
      const int* blockFromCellRaw=&blockFromCellVec[0];
      //#pragma omp target teams distribute parallel for
      for (int icell=0; icell<self_nLocal_; icell++)
      {
         dVmRaw[icell] += dVmBlockRaw[blockFromCellRaw[icell]];
      }
   }
}

unsigned* OpenmpGpuFlatDiffusion::blockIndex() {return NULL;}
double* OpenmpGpuFlatDiffusion::VmBlock() {return NULL;}
double* OpenmpGpuFlatDiffusion::dVmBlock() {return NULL;}


