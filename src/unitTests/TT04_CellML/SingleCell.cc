#include <cstdio>
#include <iostream>
#include <cassert>
#include <cmath>

using namespace std;

#include "TT04_CellML.hh"
#include "TT04_CellML_Endo.hh"
#include "TT04_CellML_Mid.hh"
#include "TT04_CellML_Epi.hh"




int main(int argc, char *argv[])
{
   if (argc < 8)
   {
      cout << "program agruments:" << endl;
      cout << "argv[1] - ic file" << endl;
      cout << "argv[2] - amplitude of stimulus -52.0" << endl;
      cout << "argv[3] - length of stimulus [time steps]    100 (* 0.02 = 2 ms)" << endl;
      cout << "argv[4] - frequency of stimulus [time steps] 50000 (* 0.02 = 1000 ms, ie. 1Hz)" << endl;
      cout << "argv[5] - simulation time [ms]  1000.0 (= 1s)" << endl;
      cout << "argv[6] - time step [ms]        0.02"<< endl;
      cout << "argv[7] - print Vm every N time steps   50" << endl;
      cout << "argv[8] - cell position [endo=0; mid=1; epi=2]" << endl;
      return 0;
   }

   char *icFile = (argv[1]);
   double stimMagnitude = (double) atof(argv[2]);
   int stimLength = atoi(argv[3]);
   int stim_cycleLength = atoi(argv[4]);
   double tend = (double) atof(argv[5]);
   double dt = (double) atof(argv[6]);
   int printVMtimesteps = atoi(argv[7]);
//  double IstimLim = (double) atof(argv[8]);
   int cellPosition = atoi(argv[8]);
   int i = printVMtimesteps;                


   double tstart = 0.0;
   //double tcurrent = 0.0;
  
   bool stimulate = false;
   int stim_event = 10000; 	// first stimulation event
   int time_cnt = 1;
   int stim_cnt = 0;


   double Istim = 0.0;
   double Vm = -86.2;
   double Iion = 0.0;  
  
  
   cout << "# ev File: " << icFile
	<< "\tstimMagnitude " << stimMagnitude
	<< "\tstimLength " << stimLength
	<< "\tstim_cycleLength " << stim_cycleLength
	<< "\ttend " << tend
	<< "\tdt " << dt
	<< endl;
                                            
   TT04_CellML* cellModel;
   switch (cellPosition)
   {
     case 0:
      cellModel = new TT04_CellML_Endo();
      break;
     case 1:
      cellModel = new TT04_CellML_Mid();
      break;
     case 2:
      cellModel = new TT04_CellML_Epi();
      break;
     default:
      assert(false);
   }
   
         
	 
   printf("# Starting the computation time loop\n");
   for(double tcurrent = tstart; tcurrent <= tend; tcurrent += dt,i++)
   {
      if (stimulate) { stim_cnt++; }
//   if (Istim > IstimLim) { cout << "High current Istim = " << Istim << endl; Istim = IstimLim; }
      if (time_cnt == stim_event)
      {
//      cout << "Trigger stimulus at t = " << tcurrent << endl;
	 stimulate = true;
         stim_event += stim_cycleLength;      //update next stim_event
      }

      if (stimulate)
      {
         if (stim_cnt < stimLength)
         { 
            //if ((tcurrent*1000) > 500.0)
            Istim = stimMagnitude;
         }
         else 
         {
            Istim = 0;
            stim_cnt = 0;
            stimulate = false;
         }
      }
	  

      double dVExternal = -Istim;
      double dVReaction = cellModel->calc(dt, Vm, Istim);
       
      Vm += (dt *  (dVExternal + dVReaction));               // In CellML it is just dVm/dt = (-1)*(Istim + Iion)
       
      
      
      if ((i == printVMtimesteps)) 
      {
	 printf("%f %le %le %le\n",tcurrent,Vm,Istim,dVReaction); 
	 fflush (stdout);
	 i = 0;
      }
      time_cnt++;
   } // end time loop

   return 0;

}
