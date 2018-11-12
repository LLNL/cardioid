#include "Grid3DStencil.hh"
#include "IndexToTuple.hh"
#include "TupleToIndex.hh"

using std::vector;

Grid3DStencil::Grid3DStencil(Long64 gid, int nx, int ny, int nz)
{
   const bool faces = true;
   const bool edges = true;
   const bool corners = false;
   nbrGids_.reserve(18);
   
   IndexToTuple gidToTuple(nx, ny, nz);
   TupleToIndex tupleToGid(nx, ny, nz);
   
   Tuple tt = gidToTuple(gid);
   if (faces)
   {
      if (tt.x() < nx-1)
      {
         Long64 newGid = tupleToGid(tt.x()+1, tt.y(), tt.z());
         nbrGids_.push_back(newGid);
      }
      if (tt.x() > 0 )
      {
         Long64 newGid = tupleToGid(tt.x()-1, tt.y(), tt.z());
         nbrGids_.push_back(newGid);
      }
      if (tt.y() < ny-1)
      {
         Long64 newGid = tupleToGid(tt.x(), tt.y()+1, tt.z());
         nbrGids_.push_back(newGid);
      }
      if (tt.y() > 0 )
      {
         Long64 newGid = tupleToGid(tt.x(), tt.y()-1, tt.z());
         nbrGids_.push_back(newGid);
      }
      if (tt.z() < nz-1)
      {
         Long64 newGid = tupleToGid(tt.x(), tt.y(), tt.z()+1);
         nbrGids_.push_back(newGid);
      }
      if (tt.z() > 0 )
      {
         Long64 newGid = tupleToGid(tt.x(), tt.y(), tt.z()-1);
         nbrGids_.push_back(newGid);
      }    
   }
   if (edges)
   {
      if (tt.x() < nx-1 && tt.y() < ny-1)
      {
         Long64 newGid = tupleToGid(tt.x()+1, tt.y()+1, tt.z());
         nbrGids_.push_back(newGid);
      }
      if (tt.x() > 0 && tt.y() < ny-1)
      {
         Long64 newGid = tupleToGid(tt.x()-1, tt.y()+1, tt.z());
         nbrGids_.push_back(newGid);
      }
      if (tt.x() < nx-1 && tt.y() > 0)
      {
         Long64 newGid = tupleToGid(tt.x()+1, tt.y()-1, tt.z());
         nbrGids_.push_back(newGid);
      }
      if (tt.x() > 0 && tt.y() > 0)
      {
         Long64 newGid = tupleToGid(tt.x()-1, tt.y()-1, tt.z());
         nbrGids_.push_back(newGid);
      }
      if (tt.x() < nx-1 && tt.z() < nz-1)
      {
         Long64 newGid = tupleToGid(tt.x()+1, tt.y(), tt.z()+1);
         nbrGids_.push_back(newGid);
      }
      if (tt.x() > 0 && tt.z() < nz-1)
      {
         Long64 newGid = tupleToGid(tt.x()-1, tt.y(), tt.z()+1);
         nbrGids_.push_back(newGid);
      }
      if (tt.x() < nx-1 && tt.z() > 0)
      {
         Long64 newGid = tupleToGid(tt.x()+1, tt.y(), tt.z()-1);
         nbrGids_.push_back(newGid);
      }
      if (tt.x() > 0 && tt.z() > 0)
      {
         Long64 newGid = tupleToGid(tt.x()-1, tt.y(), tt.z()-1);
         nbrGids_.push_back(newGid);
      }
      if (tt.y() < ny-1 && tt.z() < nz-1)
      {
         Long64 newGid = tupleToGid(tt.x(), tt.y()+1, tt.z()+1);
         nbrGids_.push_back(newGid);
      }
      if (tt.y() > 0 && tt.z() < nz-1)
      {
         Long64 newGid = tupleToGid(tt.x(), tt.y()-1, tt.z()+1);
         nbrGids_.push_back(newGid);
      }
      if (tt.y() < ny-1 && tt.z() > 0)
      {
         Long64 newGid = tupleToGid(tt.x(), tt.y()+1, tt.z()-1);
         nbrGids_.push_back(newGid);
      }
      if (tt.y() > 0 && tt.z() > 0)
      {
         Long64 newGid = tupleToGid(tt.x(), tt.y()-1, tt.z()-1);
         nbrGids_.push_back(newGid);
      }
   }
   if (corners)
   {
      if (tt.x() < nx-1 && tt.y() < ny-1 && tt.z() < nz-1)
      {
         Long64 newGid = tupleToGid(tt.x()+1, tt.y()+1, tt.z()+1);
         nbrGids_.push_back(newGid);
      }
      if (tt.x() > 0 && tt.y() < ny-1 && tt.z() < nz-1)
      {
         Long64 newGid = tupleToGid(tt.x()-1, tt.y()+1, tt.z()+1);
         nbrGids_.push_back(newGid);
      }
      if (tt.x() > 0 && tt.y() > 0 && tt.z() < nz-1)
      {
         Long64 newGid = tupleToGid(tt.x()-1, tt.y()-1, tt.z()+1);
         nbrGids_.push_back(newGid);
      }
      if (tt.x() > 0 && tt.y() < ny-1 && tt.z() > 0)
      {
         Long64 newGid = tupleToGid(tt.x()-1, tt.y()+1, tt.z()-1);
         nbrGids_.push_back(newGid);
      }
      if (tt.x() < nx-1 && tt.y() > 0 && tt.z() < nz-1)
      {
         Long64 newGid = tupleToGid(tt.x()+1, tt.y()-1, tt.z()+1);
         nbrGids_.push_back(newGid);
      }
      if (tt.x() < nx-1 && tt.y() > 0 && tt.z() > 0)
      {
         Long64 newGid = tupleToGid(tt.x()+1, tt.y()-1, tt.z()-1);
         nbrGids_.push_back(newGid);
      }
      if (tt.x() < nx-1 && tt.y() < ny-1 && tt.z() > 0)
      {
         Long64 newGid = tupleToGid(tt.x()+1, tt.y()+1, tt.z()-1);
         nbrGids_.push_back(newGid);
      }
      if (tt.x() > 0 && tt.y() > 0 && tt.z() > 0)
      {
         Long64 newGid = tupleToGid(tt.x()-1, tt.y()-1, tt.z()-1);
         nbrGids_.push_back(newGid);
      }
   }
}
