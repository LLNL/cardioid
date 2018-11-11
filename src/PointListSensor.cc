#include "PointListSensor.hh"
#include "Anatomy.hh"
#include "ioUtils.h"
#include "Simulate.hh"

#include <mpi.h>
#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

PointListSensor::PointListSensor(const SensorParms& sp, const PointListSensorParms& p, 
                                 const Anatomy& anatomy,
                                 const PotentialData& vdata)
: Sensor(sp),
  printDerivs_(p.printDerivs),
  vdata_(vdata)
{
  int myRank;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &myRank);

  const int plistsize = p.pointList.size();
  assert(plistsize > 0);

  vector<string> outfiles_loc;  // filenames of output files owned by this task
  vector<int> gidfound(plistsize,0);
  
  // loop through grid points on this task, check against pointlist
  for (unsigned ii=0; ii<plistsize; ++ii)
  {
    const Long64& gid = p.pointList[ii];
    for (unsigned jj=0; jj<anatomy.nLocal(); ++jj)
    {
      if (anatomy.gid(jj) == gid)
      {
        localCells_.push_back(gid);
        sensorind_.push_back(jj);
        ostringstream ossnum;
        ossnum.width(5);
        ossnum.fill('0');
        ossnum << ii;
        string filename = p.dirname + "/" + p.filename + "." + ossnum.str();
        outfiles_loc.push_back(filename);
        gidfound[ii] = 1;
      }
    }
  }

  //ewd:  error checking, print out warning when no cell of this gid is found (i.e. type is probably zero)
  vector<int> gidsum(plistsize,0);
  MPI_Allreduce(&gidfound[0], &gidsum[0], plistsize, MPI_INT, MPI_SUM, comm);
  if (myRank == 0)
  {
    for (unsigned ii=0; ii<plistsize; ++ii)
    {
      if (gidsum[ii] == 0)
        cout << "Warning:  PointListSensor could not find non-zero cell type with gid " << p.pointList[ii] << "!  Skipping." << endl;
      else if (gidsum[ii] > 1)
        cout << "ERROR:  PointListSensor found multiple processors which own cell gid = " << p.pointList[ii] << "!" << endl;
    }
  }
      
  if (myRank == 0)
     DirTestCreate(p.dirname.c_str());
  MPI_Barrier(comm); // none shall pass before task 0 creates directory
  
  // loop through local files, initialize ofstream, print header
  if (outfiles_loc.size() > 0)
  {
    for (unsigned ii=0; ii<outfiles_loc.size(); ++ii)
    {
      ofstream* fout_ii = new ofstream;
      fout_loc_.push_back(fout_ii);
      fout_loc_[ii]->open(outfiles_loc[ii].c_str(),ofstream::out);
      fout_loc_[ii]->setf(ios::scientific,ios::floatfield);
      if (printDerivs_)
        (*fout_loc_[ii]) << "#    time   V_m  dVm_r  dVm_d   for grid point " << localCells_[ii] << endl;
      else
        (*fout_loc_[ii]) << "#    time   V_m   for grid point " << localCells_[ii] << endl;
    }
  }
}

PointListSensor::~PointListSensor()
{
  for (unsigned ii=0; ii<fout_loc_.size(); ++ii)
  {
    fout_loc_[ii]->close();
    delete fout_loc_[ii];
  }
}

void PointListSensor::print(double time, int /*loop*/)
{
   if (printDerivs_)
      printDerivs(time);
   else
      print(time);
}


void PointListSensor::print(double time)
{
   ro_array_ptr<double> VmArray = vdata_.VmTransport_.useOn(CPU);
   for (unsigned ii=0; ii<fout_loc_.size(); ++ii)
   {
      int ind = sensorind_[ii];
      (*fout_loc_[ii]) << setprecision(10) << " " << time << "     " << VmArray[ind] << endl;
   }
}

void PointListSensor::printDerivs(double time)
{
   ro_array_ptr<double> VmArray = vdata_.VmTransport_.useOn(CPU);
   ro_array_ptr<double> dVmReaction = vdata_.dVmReactionTransport_.useOn(CPU);
   ro_array_ptr<double> dVmDiffusion = vdata_.dVmDiffusionTransport_.useOn(CPU);
   for (unsigned ii=0; ii<fout_loc_.size(); ++ii)
   {
      int ind = sensorind_[ii];
      (*fout_loc_[ii]) << setprecision(10) << " " << time << "     " << VmArray[ind] 
                       << "   " << dVmReaction[ind] 
                       << "   " << dVmDiffusion[ind] << endl;
   }
}
