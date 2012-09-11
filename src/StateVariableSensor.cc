#include "StateVariableSensor.hh"
#include "Anatomy.hh"
#include "ioUtils.h"
#include "Simulate.hh"
#include "Reaction.hh"
#include "GridPoint.hh"

#include <mpi.h>
#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

StateVariableSensor::StateVariableSensor(const SensorParms& sp, const StateVariableSensorParms& p, 
                                         const Simulate& sim)
    : Sensor(sp), startTime_(p.startTime), endTime_(p.endTime),
      gidCenter_(p.gidCenter), radius_(p.radius), sim_(sim)
{
  int myRank;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &myRank);


  // if fieldList[0] == "all", compute fieldNames from sim.anatomy
  assert(p.fieldList.size() >= 0);
  if (p.fieldList[0] == "all")
  {  
     vector<string> fieldUnits_;
     sim.reaction_->getCheckpointInfo(fieldNames_, fieldUnits_);
  }
  else
  {
     fieldNames_ = p.fieldList;        
  }
  handles_ = sim.reaction_->getVarHandle(fieldNames_);

  // loop through grid points on this task, determine if they are within radius of gidCenter
  GridPoint ctrPt(gidCenter_,sim.nx_,sim.ny_,sim.nz_);
  vector<string> outfiles_loc;  // filenames of output files owned by this task
  for (unsigned jj=0; jj<sim.anatomy_.nLocal(); ++jj)
  {
     Long64 gid = sim.anatomy_.gid(jj);
     GridPoint gpt(gid,sim.nx_,sim.ny_,sim.nz_);
     double locRadSq = (gpt.x-ctrPt.x)*(gpt.x-ctrPt.x) + (gpt.y-ctrPt.y)*(gpt.y-ctrPt.y) + (gpt.z-ctrPt.z)*(gpt.z-ctrPt.z);
     if ( gid == gidCenter_ || locRadSq <= radius_*radius_)
     {
        localCells_.push_back(gid);
        sensorind_.push_back(jj);

        // save output filename for this gid
        ostringstream ossnum;
        ossnum.width(10);
        ossnum.fill('0');
        ossnum << gid;
        string filename = p.dirname + "/" + ossnum.str() + ".sv.dat";
        outfiles_loc.push_back(filename);
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
      GridPoint gpt(localCells_[ii],sim.nx_,sim.ny_,sim.nz_);
      (*fout_loc_[ii]) << "#    gid " << localCells_[ii] << " (" << gpt.x << "," << gpt.y << "," << gpt.z << ")" << endl;
      (*fout_loc_[ii]) << "#    time    Vm    ";
      for (int jj=0; jj<fieldNames_.size(); jj++)
         (*fout_loc_[ii]) << fieldNames_[jj] << " ";
      (*fout_loc_[ii]) << endl;
    }
  }
}

StateVariableSensor::~StateVariableSensor()
{
  for (unsigned ii=0; ii<fout_loc_.size(); ++ii)
  {
    fout_loc_[ii]->close();
    delete fout_loc_[ii];
  }
}

void StateVariableSensor::print(double time, int /*loop*/)
{
   print(time);
}


void StateVariableSensor::print(double time)
{
  if (time >= startTime_ && (endTime_ <= 0.0 || time <= endTime_))
  {
     for (unsigned ii=0; ii<fout_loc_.size(); ++ii)
     {
        int kk = sensorind_[ii];
        (*fout_loc_[ii]) << setprecision(10) << " " << time << "  " << sim_.vdata_.VmArray_[kk] << "  ";
        vector<double> value(handles_.size(), 0.0);
        sim_.reaction_->getValue(ii, handles_, value);
        for (unsigned kk=0; kk<value.size(); ++kk)
           (*fout_loc_[ii]) << setprecision(10) << value[kk] << "  ";
        (*fout_loc_[ii]) << endl;
     }
  }
}

