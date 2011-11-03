#include <cstdio>
#include <unistd.h>
#include <iostream>

#include <cmath>

using namespace std;

#include "IBM_TT04.hh"

#ifndef PI
#define PI 3.14159265358979323846
#endif


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
   double totalCalcTime = 0.0;
  
   bool stimulate = false;
   int stim_event = 100; 	// first stimulation event
   int time_cnt = 1;
   int stim_cnt = 0;


   double Istim = 0.0;
   double Vm = -86.2;
   double Iion = 0.0;  
  
  
   cout << "ev File: " << icFile
	<< "\tstimMagnitude " << stimMagnitude
	<< "\tstimLength " << stimLength
	<< "\tstim_cycleLength " << stim_cycleLength
	<< "\ttend " << tend
	<< "\tdt " << dt
	<< endl;
                                            
   IBM_TT04 *pemIBM;
   pemIBM = new IBM_TT04(icFile, cellPosition);  
         
	 
   double MinusOne_over_Cm = (-1)/((pemIBM)->getCm());	                                              
   printf("Starting the computation time loop\n");
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
      double dVReaction = (pemIBM)->Calc(dt, (Vm), (Istim));
       
//       Vm += (dt * (MinusOne_over_Cm) * (Istim + Iion)); // Eq. 1 tenTusscher et al. AJP Heart Circ Physiol 2004
      Vm += (dt *  (dVExternal + dVReaction));               // In CellML it is just dVm/dt = (-1)*(Istim + Iion)
       
      
      
      if ((i == printVMtimesteps)) 
      {
//	    printf("t %d Vm %lf stim %lf ion %lf K1 %lf to %lf Kr %lf Ks %lf CaL %lf NaK %lf Na %lf NaCa %lf bCa %lf pK %lf pCa %lf leak %lf up %lf rel %lf\n",
	 printf("%f %le %le %le\n",tcurrent,Vm,Istim,Iion); 
/*
  printf("%d %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",
  (int) tcurrent, Vm, Istim, (pemIBM)->getIion(), (pemIBM)->getIK1(), (pemIBM)->getIto(), (pemIBM)->getIKr(), (pemIBM)->getIKs(),
  (pemIBM)->getICaL(), (pemIBM)->getINaK(), (pemIBM)->getINa(), (pemIBM)->getIbNa(), (pemIBM)->getINaCa(), (pemIBM)->getIbCa(),
  (pemIBM)->getIpK(), (pemIBM)->getIpCa(), (pemIBM)->getIleak(), (pemIBM)->getIup(), (pemIBM)->getIrel());
*/
	 fflush (stdout);
	 i = 0;
      }
      time_cnt++;
   } // end time loop

   return 0;

}
/*      
	inline double getV(){return (y_TT04[tt04_V]);};
	inline double getNai(){return (y_TT04[tt04_Nai]);};
	inline double getKi(){return (y_TT04[tt04_Ki]);};
	inline double getCai(){return (y_TT04[tt04_Cai]);};
	inline double getXr1(){return (y_TT04[tt04_xr1]);};
	inline double getXr2(){return (y_TT04[tt04_xr2]);};
	inline double getXs(){return (y_TT04[tt04_xs]);};
	inline double getm(){return (y_TT04[tt04_m]);};
	inline double geth(){return (y_TT04[tt04_h]);};
	inline double getj(){return (y_TT04[tt04_j]);};
	inline double getd(){return (y_TT04[tt04_d]);};
	inline double getf(){return (y_TT04[tt04_f]);};
	inline double getfCa(){return (y_TT04[tt04_fCa]);};
	inline double gets(){return (y_TT04[tt04_s]);};
	inline double getr(){return (y_TT04[tt04_r]);};
	inline double getCaSR(){return (y_TT04[tt04_CaSR]);};
	inline double getg(){return (y_TT04[tt04_g]);};	      
*/	      
