#include "writeCells.hh"

#include <cassert>
#include <iostream>
#include <algorithm>

#include "AnatomyReader.hh"
#include "pio.h"
#include "IndexToVector.hh"

using std::vector;
using namespace std;



void writeCells(const AnatomyReader& anatomy,
		const std::string& filename)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   const vector<AnatomyCell>& cells = anatomy._anatomy;
   IndexToVector indexToVector(anatomy._nx,
			       anatomy._ny,
			       anatomy._nz);

   int halfNx = anatomy._nx/2;
   int halfNy = anatomy._ny/2;
   int halfNz = anatomy._nz/2;
   

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
      Pprintf(file, "  h = %4u  0    0\n", anatomy._nx);
      Pprintf(file, "        0    %4u  0\n", anatomy._ny);
      Pprintf(file, "        0    0    %4u;\n", anatomy._nz);
      Pprintf(file, "}\n\n");
   }
   
   char line[lrec+1];
   int lastDest = -1;
   int count = 0;
   for (unsigned ii=0; ii<cells.size(); ++ii)
   {
      if (cells[ii]._dest != lastDest)
      {
	 count =0;
	 lastDest = cells[ii]._dest;
      }

//       if (count++ > 1000)
// 	 continue;
      

      THREE_VECTOR v = indexToVector(cells[ii]._gid);
      int ix = int(v.x) - halfNx;
      int iy = int(v.y) - halfNy;
      int iz = int(v.z) - halfNz;
//       int ix = int(v.x);
//       int iy = int(v.y);
//       int iz = int(v.z);

      int l = snprintf(line, lrec, fmt,
		   ix, iy, iz, cells[ii]._cellType, cells[ii]._dest);
      
      for (; l < lrec - 1; l++) line[l] = (char)' ';
      line[l++] = (char)'\n';
      assert (l==lrec);
      Pwrite(line, lrec, 1, file);
   }
   
   Pclose(file);
}
