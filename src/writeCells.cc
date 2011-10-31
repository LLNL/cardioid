#include "writeCells.hh"

#include <cassert>
#include <iostream>
#include <algorithm>

#include "Simulate.hh"
#include "Anatomy.hh"
#include "pio.h"
#include "IndexToVector.hh"

using std::vector;
using namespace std;


/** Only write records for the local cells */
void writeCells(const Simulate& sim,
		const std::string& filename)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   const Anatomy& anatomy = sim.anatomy_;
   
   const vector<AnatomyCell>& cells = anatomy.cellArray();
   IndexToVector indexToVector(anatomy.nx(), anatomy.ny(), anatomy.nz());

   int halfNx = anatomy.nx()/2;
   int halfNy = anatomy.ny()/2;
   int halfNz = anatomy.nz()/2;
   

   Long64 nLocal = anatomy.nLocal();
   Long64 nGlobal;
   MPI_Allreduce(&nLocal, &nGlobal, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
   
   PFILE* file = Popen(filename.c_str(), "w", MPI_COMM_WORLD);

   char fmt[] = "%5d %5d %5d %4u %8u %18.12f";
   int lrec = 55;
   int nfields = 6; 


   if (myRank == 0)
   {
      // write header
      int nfiles;
      Pget(file,"ngroup",&nfiles);
      Pprintf(file, "cellViz FILEHEADER {\n");
      Pprintf(file, "  lrec = %d;\n", lrec);
      Pprintf(file, "  datatype = FIXRECORDASCII;\n");
      Pprintf(file, "  nrecords = %llu;\n", nGlobal);
      Pprintf(file, "  nfields = %d;\n", nfields);
      Pprintf(file, "  field_names = rx ry rz cellType domain Vm;\n");
      Pprintf(file, "  field_types = u u u u u f;\n" );
      Pprintf(file, "  nfiles = %u;\n", nfiles);
      Pprintf(file, "  h = %4u  0    0\n", anatomy.nx());
      Pprintf(file, "        0    %4u  0\n", anatomy.ny());
      Pprintf(file, "        0    0    %4u;\n", anatomy.nz());
      Pprintf(file, "}\n\n");
   }
   
   char line[lrec+1];
   int lastDest = -1;
   int count = 0;
   for (unsigned ii=0; ii<nLocal; ++ii)
   {
      if (cells[ii].dest_ != lastDest)
      {
	 count =0;
	 lastDest = cells[ii].dest_;
      }

      Vector v = indexToVector(cells[ii].gid_);
      int ix = int(v.x()) - halfNx;
      int iy = int(v.y()) - halfNy;
      int iz = int(v.z()) - halfNz;

      int l = snprintf(line, lrec, fmt,
		       ix, iy, iz, cells[ii].cellType_, cells[ii].dest_,
		       sim.VmArray_[ii]);
      
      for (; l < lrec - 1; l++) line[l] = (char)' ';
      line[l++] = (char)'\n';
      assert (l==lrec);
      Pwrite(line, lrec, 1, file);
   }
   
   Pclose(file);
}

/** Writes all the cells in the cell array. */
void writeCells(const vector<AnatomyCell>& cells,
		int nx, int ny, int nz,
		const std::string& filename)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   IndexToVector indexToVector(nx, ny, nz);

   int halfNx = nx/2;
   int halfNy = ny/2;
   int halfNz = nz/2;
   

   Long64 nLocal = cells.size();
   Long64 nGlobal;
   MPI_Allreduce(&nLocal, &nGlobal, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
   
   PFILE* file = Popen(filename.c_str(), "w", MPI_COMM_WORLD);

   char fmt[] = "%5d %5d %5d %4u %8u";
   int lrec = 35;
   int nfields = 5; 


   if (myRank == 0)
   {
      // write header
      int nfiles;
      Pget(file,"ngroup",&nfiles);
      Pprintf(file, "cellViz FILEHEADER {\n");
      Pprintf(file, "  lrec = %d;\n", lrec);
      Pprintf(file, "  datatype = FIXRECORDASCII;\n");
      Pprintf(file, "  nrecords = %llu;\n", nGlobal);
      Pprintf(file, "  nfields = %d;\n", nfields);
      Pprintf(file, "  field_names = rx ry rz cellType domain;\n");
      Pprintf(file, "  field_types = u u u u u;\n" );
      Pprintf(file, "  nfiles = %u;\n", nfiles);
      Pprintf(file, "  h = %4u  0    0\n", nx);
      Pprintf(file, "        0    %4u  0\n", ny);
      Pprintf(file, "        0    0    %4u;\n", nz);
      Pprintf(file, "}\n\n");
   }
   
   char line[lrec+1];
   int lastDest = -1;
   int count = 0;
   for (unsigned ii=0; ii<nLocal; ++ii)
   {
      if (cells[ii].dest_ != lastDest)
      {
	 count =0;
	 lastDest = cells[ii].dest_;
      }

      Vector v = indexToVector(cells[ii].gid_);
      int ix = int(v.x()) - halfNx;
      int iy = int(v.y()) - halfNy;
      int iz = int(v.z()) - halfNz;

      int l = snprintf(line, lrec, fmt,
		       ix, iy, iz, cells[ii].cellType_, cells[ii].dest_);
      
      for (; l < lrec - 1; l++) line[l] = (char)' ';
      line[l++] = (char)'\n';
      assert (l==lrec);
      Pwrite(line, lrec, 1, file);
   }
   
   Pclose(file);
}   
