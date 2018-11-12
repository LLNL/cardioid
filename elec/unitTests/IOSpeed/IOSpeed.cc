#include <mpi.h>

#include "heap.h"
#include "pio.h"
#include "hardwareInfo.h"
#include "mpiUtils.h"
#include "ioUtils.h"

#include <iostream>
#include <vector>
#include <utility>
#include <sstream>

using namespace std;

MPI_Comm COMM_LOCAL;

void ioTest(int kBytes, int nFiles)
{
   Pio_setNumWriteFiles(nFiles);
   vector<char> data(kBytes*1024);
   stringstream dirname;
   dirname << kBytes <<"kPerTask_"<<nFiles<<"files";


   MPI_Comm comm = MPI_COMM_WORLD;
   int myRank;
   MPI_Comm_rank(comm, &myRank);
   if (myRank == 0)
      DirTestCreate(dirname.str().c_str());
   string filename = dirname.str() + "/data";
   PFILE* file = Popen(filename.c_str(), "w", comm);

   Pwrite(&data[0], 1, kBytes*1024, file);
   
   Pclose(file);
}




int main(int argc, char** argv)
{
   MPI_Init(&argc, &argv);
   COMM_LOCAL = MPI_COMM_WORLD;
   heap_start(50);
   
   int myRank;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

   int nIoTasks = hi_nIoTasks(MPI_COMM_WORLD);
   if (myRank == 0)
      cout << "nIoTasks = " << nIoTasks << endl;
   

   vector<pair<int, int> > testList; //kb/task, nFiles
   // testList.push_back(make_pair(10, 512));
   // testList.push_back(make_pair(10, 1024));
   // testList.push_back(make_pair(100, 512));
   // testList.push_back(make_pair(100, 1024));
   // testList.push_back(make_pair(1000, 512));
   // testList.push_back(make_pair(1000, 1024));
   testList.push_back(make_pair(10000, 384));
   testList.push_back(make_pair(10000, 512));
   testList.push_back(make_pair(10000, 768));
   testList.push_back(make_pair(10000, 1024));
   testList.push_back(make_pair(10000, 1450));
   testList.push_back(make_pair(10000, 1500));
   testList.push_back(make_pair(10000, 1536));
   

   
   
   for (unsigned ii=0; ii<testList.size(); ++ii)
   {
      ostringstream startMsg;
      ostringstream endMsg;
      int kBytesPerTask = testList[ii].first;
      int nFiles = testList[ii].second;
      startMsg << "Starting I/O test for " <<kBytesPerTask << "k per task, "
               << nFiles << " files";
      endMsg << "Finished I/O test for " <<kBytesPerTask << "k per task, "
               << nFiles << " files";
      timestampBarrier(startMsg.str().c_str(), MPI_COMM_WORLD);
      ioTest(kBytesPerTask, nFiles);
      timestampBarrier(endMsg.str().c_str(), MPI_COMM_WORLD);
   }
   
   
   MPI_Finalize();
}
