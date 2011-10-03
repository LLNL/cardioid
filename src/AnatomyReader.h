#ifndef ANATOMY_READER_H
#define ANATOMY_READER_H

#include <vector>
#include <string>
#include <mpi.h>

#include "AnatomyCell.h"

struct pfile_st;


class AnatomyReader
{
  public:

   AnatomyReader(const std::string& filename, MPI_Comm comm);
   
  private:

   void asciiReader(pfile_st* file);
   void binaryReader(pfile_st* file);


   int _nx; // nCells in x-direction in global anatomy grid
   int _ny; // nCells in y-direction in global anatomy grid
   int _nz; // nCells in z-direction in global anatomy grid
   std::vector<std::string> _cellTypes;
   std::vector<AnatomyCell> _anatomy;
   
};

#endif
