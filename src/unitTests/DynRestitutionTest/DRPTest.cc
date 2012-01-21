#include <cstdio>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>

#include "Anatomy.hh"
#include "Reaction.hh"
#include "TT04_CellML_Reaction.hh"
#include "TT04Dev_Reaction.hh"
#include "TT06_CellML_Reaction.hh"
#include "TT06Dev_Reaction.hh"
#include "TT06_RRG_Reaction.hh"
#include "ReactionFHN.hh"

using namespace std;


Anatomy buildAnatomy(int cellType)
{
  AnatomyCell c;
  c.gid_ = 0;
  c.cellType_ = 100+cellType;
  c.theta_ = 0;
  c.phi_ = 0;
  Anatomy a;
  a.setGridSize(1, 1, 1);
  a.cellArray().push_back(c);
  a.nLocal() = 1;
  a.nRemote() = 0;
  a.dx() = 0.2;
  a.dy() = 0.2;
  a.dz() = 0.2;
  return a;
}

Reaction* factory(const string& name, const Anatomy& anatomy, double tolerance,int mod)
{
   if (name == "cellml_tt04")    return new TT04_CellML_Reaction(anatomy, TT04_CellML_Reaction::rushLarsen);
   if (name == "cellml_tt04_fe") return new TT04_CellML_Reaction(anatomy, TT04_CellML_Reaction::forwardEuler);
   if (name == "cellml_tt06")    return new TT06_CellML_Reaction(anatomy, TT06_CellML_Reaction::rushLarsen);
   if (name == "cellml_tt06_fe") return new TT06_CellML_Reaction(anatomy, TT06_CellML_Reaction::forwardEuler);
   if (name == "tt06rrg")        return new TT06_RRG_Reaction(anatomy);
   if (name == "fhn")            return new ReactionFHN(anatomy);
   if (name == "tt04dev")        return new TT04Dev_Reaction(anatomy);
   if (name == "tt06dev")        return new TT06Dev_Reaction(anatomy,tolerance,mod);
   assert(false);
   return 0;
}


int main(int argc, char *argv[])
{

  const double apdThresh = 0.90;  // fraction of repolarization used to define apd time
  const double stimOffset = 0.0;  // delay onset of stim at each beat, msec
  const int printRate = 50;       // frequency of voltage printing
  
  if (argc < 11)
  {
    cout << "program arguments:" << endl;
    cout << "argv[1] - method name (see list below)" << endl;
    cout << "argv[2] - amplitude of stimulus -52.0" << endl;
    cout << "argv[3] - length of stimulus [ms] 1 ms" << endl;
    cout << "argv[4] - time step [ms] 2e-4 ms"<< endl;
    cout << "argv[5] - initial basic cycle length [ms] 4000 ms" << endl;
    cout << "argv[6] - final basic cycle length [ms] 500 ms" << endl;
    cout << "argv[7] - beats per restitution block  50" << endl;
    cout << "argv[8] - extra time between blocks [ms]  0 ms" << endl; 
    cout << "argv[9] - cell position [endo=0; mid=1; epi=2]" << endl;
    cout << "argv[10] - tolerance for pade approximations" << endl;
    cout << "argv[11] - mod (used by tt06dev only)" << endl;
      cout <<endl;
      cout << "Supported cell models:" <<endl;
      cout << "----------------------" <<endl;
      cout << "   cellml_tt04      TT04 from CellML.  Rush-Larsen integrator" << endl;
      cout << "   cellml_tt04_fe   TT04 from CellML.  Forward Euler integrator" << endl;
      cout << "   cellml_tt06      TT06 from CellML.  Rush-Larsen integrator" << endl;
      cout << "   cellml_tt06_fe   TT06 from CellML.  Forward Euler integrator" << endl;
      cout << "   tt06rrg          TT06 as modified by Rice et al." << endl;
      cout << "   fhn              FitzHugh-Nagumo" << endl;
      cout << "   tt04dev          Developmental version of TT04" << endl;
      cout << "   tt06dev          Developmental version of TT06" << endl;
    return 0;
  }

  // if true, define baseline voltage from initial value only
  // if false, reset baseline voltage to value just before next stimulus pulse
  const bool fixedBaseline = true;  
                                    
  string method =            (argv[1]);
  double stimMagnitude =     atof(argv[2]);
  double stimLength =        atof(argv[3]);
  double dt =                atof(argv[4]);
  double initCycleLength =   atof(argv[5]);
  double finalCycleLength =  atof(argv[6]);
  int nBeatsPerBlock =       atoi(argv[7]);
  double timeBetween =       atof(argv[8]);
  int cellPosition =         atoi(argv[9]);
  double tolerance =        atof(argv[10]);
  int mod = 0;
  if (argc > 11) 
     mod =                atoi(argv[11]);
  
  Anatomy anatomy = buildAnatomy(cellPosition);
  Reaction* cellModel = factory(method, anatomy,tolerance,mod);
     
  cout << "# method: " << method
       << "\tstimMagnitude " << stimMagnitude
       << "\tstimStart " << stimLength
       << "\tstimLength " << stimLength
       << "\tdt " << dt
       << "\tBCL (basic cycle length) range: "
       << initCycleLength << " - " << finalCycleLength
       << "\t nBeatsPerBlock " << nBeatsPerBlock
       << "\t timeBetween " << timeBetween
       << "\tcell type " << cellPosition
       << endl;

  unsigned nLocal = anatomy.nLocal();
  vector<double> Vm(nLocal, -86.2);
  vector<double> dVmReaction(nLocal, 0);
  vector<double> iStim(nLocal, 0.0);
  double baseline = Vm[0];
  double lastVm = 0.0;
  double time = 0.0;
  ofstream vmOutfile("vm_vs_t.dat");
  ofstream apdOutfile("apd_vs_di.dat");
  vmOutfile << "# time Vm iStim dVmReaction" << endl;
  apdOutfile << "# DI time   APD time    simulation time t    Vm(t)" << endl;
  cout << "# DI time   APD time  time t    Vm(t)   cycleLength" << endl;
  vmOutfile.setf(ios::scientific,ios::floatfield);
  apdOutfile.setf(ios::scientific,ios::floatfield);
  cout.setf(ios::scientific,ios::floatfield);
  
  double cycleLength = initCycleLength;
  while (cycleLength > finalCycleLength)
  {
    // restitution block of nBeats
    bool apdThreshSeek = false;
    bool apdThreshHit = false;
    double apdVm, apd, di;
    double stimOnTime = -1.E+30;
    double stimOffTime = +1.E+30;
    double apdThreshTime = +1.E+30;
    for (int ibeat = 0; ibeat < nBeatsPerBlock; ibeat++)
    {
      int nsteps = cycleLength/dt;
      for (int istep = 0; istep < nsteps; istep++)
      {
        double tloc = istep*dt;
        double stim = 0.0;
        if (tloc >= stimOffset && tloc <= stimOffset+stimLength)
          stim = stimMagnitude;

        if (stim != iStim[0] && iStim[0] == 0.0)
        {
          stimOnTime = time;
          if (!fixedBaseline)
            baseline = lastVm;  // use the voltage right before stimulus was applied as our baseline
        }
        else if (stim != iStim[0] && stim == 0.0)
        {
          stimOffTime = time;
          apdThreshSeek = true;
        }

        if (Vm[0] <= baseline*apdThresh && apdThreshSeek)
        {
          apdThreshHit = true;
          apdThreshSeek = false;
          apdThreshTime = time;
          apdVm = Vm[0];
          //cout << "APD Threshold hit:  " << Vm[0] << "  " << baseline*apdThresh << ", baseline = " << baseline << ", time = " << time << ", stimOffTime = " << stimOffTime << endl;
        }

        // when next stimulus is applied, print out computed APD and DI times
        if (stimOnTime > apdThreshTime && apdThreshHit)
        {
          apdThreshHit = false;
          di = stimOnTime - apdThreshTime;
          apd = apdThreshTime-stimOffTime;     
          apdOutfile << " " << di << "   " << apd << "   " << apdThreshTime << "   " << apdVm << "   " << cycleLength << endl;
        }

        if (istep%printRate == 0) 
          vmOutfile << setprecision(10) << time << "  " << Vm[0] << "  " << iStim[0] << "  " << dVmReaction[0] << endl;

        // advance cell model
        iStim.assign(nLocal, stim);
        lastVm = Vm[0];
        cellModel->calc(dt, Vm, iStim, dVmReaction);
        for (unsigned ii=0; ii<nLocal; ++ii)
          Vm[ii] += dt*(dVmReaction[ii] - iStim[ii]);
        time += dt;
      }
    }

    // delay period after last beat
    for (int istep = 0; istep < (int)(timeBetween/dt); istep++)
    {
      if (Vm[0] <= baseline*apdThresh && apdThreshSeek)
      {
        apdThreshHit = true;
        apdThreshSeek = false;
        apdThreshTime = time;
        apdVm = Vm[0];
        //cout << "APD Threshold hit:  " << Vm[0] << "  " << baseline*apdThresh << ", baseline = " << baseline << ", time = " << time << ", stimOffTime = " << stimOffTime << endl;
      }

      // when next stimulus is applied, print out computed APD and DI times
      if (stimOnTime > apdThreshTime && apdThreshHit)
      {
        apdThreshHit = false;
        di = stimOnTime - apdThreshTime;
        apd = apdThreshTime-stimOffTime;     
        apdOutfile << " " << di << "   " << apd << "   " << apdThreshTime << "   " << apdVm << endl;
      }

      if (istep%printRate == 0) 
        vmOutfile << setprecision(10) << time << "  " << Vm[0] << "  " << iStim[0] << "  " << dVmReaction[0] << endl;
      
      double tloc = istep*dt;
      iStim.assign(nLocal, 0.0);
      cellModel->calc(dt, Vm, iStim, dVmReaction);
      for (unsigned ii=0; ii<nLocal; ++ii)
        Vm[ii] += dt*dVmReaction[ii];
      time += dt;
    }

    // print last computed apd vs. di for this cycle
    cout << " " << di << "   " << apd << "   " << apdThreshTime << "   " << apdVm << "   " << cycleLength << endl;
    
    
    // decrease cycle length
    if (cycleLength > 1000.)
      cycleLength -= 1000.;
    else if (cycleLength > 500.)
      cycleLength -= 250.;
    else if (cycleLength > 400.)
      cycleLength -= 50.;
    else if (cycleLength > 300.)
      cycleLength -= 10.;
    else if (cycleLength > 250.)
      cycleLength -= 5.;
    else if (cycleLength > 50.)
      cycleLength -= 1.;
    else
    {
      cout << "Restitution protocol finished.";
      return 0;
    }
  }
  vmOutfile.close();
  apdOutfile.close();

  return 0;
}
