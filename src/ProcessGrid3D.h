#ifndef PROCESSGRID3D_H
#define PROCESSGRID3D_H

#include <vector>
#ifdef USE_MPI
#include <mpi.h> 
#else
typedef int MPI_Comm;
#endif

class ProcessGrid3D {
  private:

  MPI_Comm comm_;
  int np0_,np1_,np2_;
  int ip0_,ip1_,ip2_;
  int mype_, npes_;
  bool active_;
  
  public:

  ProcessGrid3D(MPI_Comm comm, int np0, int np1, int np2);
  ~ProcessGrid3D();

  int mype(void) { return mype_; }
  int npes(void) { return npes_; }
  int np0(void) {return np0_;}    // x dimension of process grid
  int np1(void) {return np1_;}    // y dimension of process grid
  int np2(void) {return np2_;}    // z dimension of process grid
  int ip0(void) {return ip0_;}    // x value of current process
  int ip1(void) {return ip1_;}    // y value of current process
  int ip2(void) {return ip2_;}    // z value of current process
  bool active(void) {return active_;} 
  int pe(int i, int j, int k) { return i + j*np0_ + k*np0_*np1_;}  // process number 
  MPI_Comm subcomm(int nsubx, int nsuby, int nsubz, int istart, int jstart, int kstart);
  int dsend(int destpe, int dind, std::vector<double>& vec, std::vector<int>& id);
  int drecv(int srcpe, int dind, std::vector<double>& vec, std::vector<int>& id);
};
#endif
