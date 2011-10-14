#ifndef ANATOMY_READER_H
#define ANATOMY_READER_H

#include <vector>
#include <string>
#include <mpi.h>

#include "AnatomyCell.hh"

struct pfile_st;
class Simulate;

class AnatomyReader
{
  public:

   AnatomyReader(const std::string& filename, MPI_Comm comm,
   Simulate& sim);


  private:

   void asciiReader(pfile_st* file);
   void binaryReader(pfile_st* file);

//   std::vector<std::string> _cellTypes;
   std::vector<AnatomyCell>& _anatomy;
};

#endif
