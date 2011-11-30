#include "PointListSensor.hh"
#include "Anatomy.hh"
#include <mpi.h>
#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

PointListSensor::PointListSensor(const PointListSensorParms& p, const Anatomy& anatomy)
    : startTime_(p.startTime), endTime_(p.endTime),
      filebase_(p.filebase), printRate_(p.printRate),
      printDerivs_(p.printDerivs)
{
  const int plistsize = p.pointlist.size();
  assert(plistsize > 0);

  vector<string> outfiles_loc;  // filenames of output files owned by this task
  vector<int> gidfound(plistsize,0);
  
  // loop through grid points on this task, check against pointlist
  for (unsigned ii=0; ii<plistsize; ++ii)
  {
    const unsigned gid = p.pointlist[ii];
    //for (unsigned jj=0; jj<anatomy.size(); ++jj)
    for (unsigned jj=0; jj<anatomy.nLocal(); ++jj)
    {
      if (anatomy.gid(jj) == gid)
      {
        pointlist_loc_.push_back(gid);
        sensorind_.push_back(jj);
        ostringstream ossnum;
        ossnum.width(5);
        ossnum.fill('0');
        ossnum << ii;
        string filename = filebase_ + "." + ossnum.str();
        outfiles_loc.push_back(filename);
        gidfound[ii] = 1;
      }
    }
  }

  //ewd:  error checking, print out warning when no cell of this gid is found (i.e. type is probably zero)
  vector<int> gidsum(plistsize,0);
  MPI_Allreduce(&gidfound[0], &gidsum[0], plistsize, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  int myRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  if (myRank == 0)
  {
    for (unsigned ii=0; ii<plistsize; ++ii)
    {
      if (gidsum[ii] == 0)
        cout << "Warning:  PointListSensor could not find non-zero cell type with gid " << p.pointlist[ii] << "!  Skipping." << endl;
      else if (gidsum[ii] > 1)
        cout << "ERROR:  PointListSensor found multiple processors which own cell gid = " << p.pointlist[ii] << "!" << endl;
    }
  }
      


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
        (*fout_loc_[ii]) << "#    time   V_m  dVm_r  dVm_d  dVm_e   for grid point " << pointlist_loc_[ii] << endl;
      else
        (*fout_loc_[ii]) << "#    time   V_m   for grid point " << pointlist_loc_[ii] << endl;
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

void PointListSensor::print(double time, std::vector<double>& Vm)
{
  if (time >= startTime_ && (endTime_ <= 0.0 || time <= endTime_))
  {
    for (unsigned ii=0; ii<fout_loc_.size(); ++ii)
    {
      int ind = sensorind_[ii];
      (*fout_loc_[ii]) << setprecision(10) << " " << time << "     " << Vm[ind] << endl;
    }
  }
}

void PointListSensor::print(double time, vector<double>& Vm, vector<double>& dVm_r, vector<double>& dVm_d, vector<double>& dVm_e)
{
  if (time >= startTime_ && (endTime_ <= 0.0 || time <= endTime_))
  {
    for (unsigned ii=0; ii<fout_loc_.size(); ++ii)
    {
      int ind = sensorind_[ii];
      (*fout_loc_[ii]) << setprecision(10) << " " << time << "     " << Vm[ind] << "   " << dVm_r[ind] << "   " << dVm_d[ind] << "   " << dVm_e[ind] << endl;
    }
  }
}
