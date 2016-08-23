#include <cstdlib>
#include <iostream>
#include <mpi.h>
#include <vector>
#include <set>
#include "CommTable.hh"
#include "HaloExchange.hh"

using namespace std;
#define timebase(x) asm volatile ("mftb %0" : "=r"(x) : )

//number of neighbors to talk to
const int Nneighbor=20;

//number of doubles to send per node
const int MessageSize = 70;

//variation range
const int VarRange = 37;

//locality range
const int LocalRange = 3;

//#define RANDOM

extern "C"{
uint32_t coord_list(int *myid, int *totN, int* mapping);
}


int Nnode;
int main(int argc, char** argv)
{
   uint64_t t1,t2,t3,t4;
   int npes, mype;
   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD, &npes);
   MPI_Comm_rank(MPI_COMM_WORLD, &mype);

  //get my ID
  int myID=mype;
  int sum;

  Nnode=npes;

  assert(Nnode > Nneighbor);

  //randomize 
  srand((myID+3213)*2627);
 
  set<int> sendTask_s;
  vector<int> sendOffset;

  
  //assign nodes
  #ifdef RANDOM
  while(sendTask_s.size()<Nneighbor)
  {
    sendTask_s.insert((rand()%(Nnode-1) +1+ myID )%Nnode);
  }
  vector<int> sendTask(sendTask_s.begin(),sendTask_s.end());
  #else
  int* mapping = new int[Nnode*5];
  int myid_spi;
  int totN_spi;
  int dist;
  vector<int> sendTask;
  coord_list(&myid_spi,&totN_spi,mapping);
//  cout << "mpi_rank=" << mype << " spi_rank=" << myid_spi << endl;
//  cout << "mpi_totN=" << npes << " spi_totN=" << totN_spi << endl;
  for(int ii=0;ii<totN_spi;ii++)
  {
    dist =  abs(mapping[myid_spi*5+0] - mapping[ii*5+0]) +
            abs(mapping[myid_spi*5+1] - mapping[ii*5+1]) +
            abs(mapping[myid_spi*5+2] - mapping[ii*5+2]) +
            abs(mapping[myid_spi*5+3] - mapping[ii*5+3]) +
            abs(mapping[myid_spi*5+4] - mapping[ii*5+4]) ;

    if (dist==0)
    {
      assert(ii==myid_spi);
      continue;
    }
    if ( dist <= LocalRange ) sendTask.push_back(ii);
    cout << dist << " ";
  } 
  cout << endl;
  assert(sendTask.size() >= Nneighbor);
  while(sendTask.size() > Nneighbor )
  {
    sendTask.erase(sendTask.begin()+rand()%sendTask.size());
  }
  assert(sendTask.size() == Nneighbor);
  delete [] mapping;
  #endif


  cout << "my dest is :";
  for(int ii=0 ; ii < sendTask.size() ; ++ii )
    cout << sendTask[ii] << " ";
  cout<<endl;

  sum=0;
  for(int ii=0 ; ii < sendTask.size() ; ++ii )
  {
    sendOffset.push_back(sum);
    int new_size=rand()%VarRange + MessageSize;
    new_size -= new_size%4; //why align by 4?
    sum+=new_size;
  }
  sendOffset.push_back(sum);
    
  //create sendMap
  vector<int> sendMap;
  sendMap.resize(sum);

  cout<<"my total send size is " << sum << " doubles" << endl;

  for(int ii=0;ii<sendMap.size();ii++) sendMap[ii]=ii;

  //create Data
  vector<double> Data;
  for(int ii=0;ii<sum;ii++) Data.push_back(mype*1000+ii); 
  int nLocal = sum;

  //create commTable
  CommTable my_cmt(sendTask,sendOffset,MPI_COMM_WORLD);
  
  cout<<"my total recv size is " << my_cmt.recvSize() << " doubles" << endl;

  HaloExchange<double> spi_HE(sendMap,&my_cmt);

  cout << " total time = " << (double)t3/10.0/1600.0  << " us" << endl;
  
  MPI_Finalize();

  return 0;
}
