#include "writeCells.hh"

#include <cassert>
#include <iostream>
#include <algorithm>
#include <set>

#include "Simulate.hh"
#include "Anatomy.hh"
#include "PioHeaderData.hh"
#include "pio.h"
#include "IndexToVector.hh"
#include <algorithm>

using namespace std;


/** Writes all the cells in the cell array. */
void writeCells(const vector<AnatomyCell>& cells,
                int nx, int ny, int nz,
                const std::string& filename)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   Long64 nLocal = cells.size();
   Long64 nGlobal;
   MPI_Allreduce(&nLocal, &nGlobal, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
   
   PFILE* file = Popen(filename.c_str(), "w", MPI_COMM_WORLD);

   char fmt[] = "%12llu %4u %8u";
   int lRec = 27;
   int nFields = 3; 

   PioHeaderData header;
   header.objectName_ = "cellViz";
   header.className_  = "FILEHEADER";
   header.dataType_   = PioHeaderData::ASCII;
   header.nRecords_   = nGlobal;
   header.lRec_       = lRec;
   header.nFields_    = nFields;
   header.fieldNames_ = "gid cellType domain";
   header.fieldTypes_ = "u u u";
   header.fieldUnits_ = "1 1 1";
   header.addItem("nx", nx);
   header.addItem("ny", ny);
   header.addItem("nz", nz);
   
   if (myRank == 0)
      header.writeHeader(file, 0, 0);

   char line[lRec+1];
   for (unsigned ii=0; ii<nLocal; ++ii)
   {
      int l = snprintf(line, lRec, fmt,
                       cells[ii].gid_, cells[ii].cellType_, cells[ii].dest_);
      
      for (; l < lRec - 1; l++) line[l] = (char)' ';
      line[l++] = (char)'\n';
      assert (l==lRec);
      Pwrite(line, lRec, 1, file);
   }
   
   Pclose(file);
}   
