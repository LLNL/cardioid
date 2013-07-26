#include "readCellList.hh"
#include <mpi.h>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h> 

using namespace std;

/** Initialize cellVec with gids listed in file filename */
void readCellList(const string filename, vector<Long64>& cellVec)
{
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   // try to open file
   int openfail;
   ifstream input;
   if (myRank == 0)
   {
      input.open(filename.c_str(),ifstream::in);
      if (!input.is_open())
      {
         cerr << "Could not open cell list file " << filename << endl;
         exit(1);
      }
      
   }

   if (myRank == 0)
   {
      while (!input.eof())
      {
         string query;
         if ( !getline ( input, query ) )
             break;
         istringstream ss ( query );
         Long64 igid;
         while( ss >> igid )
         {
            assert( igid>=0 );
            cellVec.push_back(igid);
         }
      }
   }   
   int nCells = cellVec.size();
   MPI_Bcast(&nCells, 1, MPI_INT, 0, MPI_COMM_WORLD);
   cellVec.resize(nCells);
   MPI_Bcast(&cellVec[0], nCells, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
}
