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
  assert(p.pointlist.size() > 0);

  vector<string> outfiles_loc;  // filenames of output files owned by this task
  
  // loop through grid points on this task, check against pointlist
  for (unsigned ii=0; ii<p.pointlist.size(); ++ii)
  {
    const unsigned gid = p.pointlist[ii];
    for (unsigned jj=0; jj<anatomy.size(); ++jj)
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
      }
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
