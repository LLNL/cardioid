#include "Reaction.hh"
#include "Interpolation.hh"
#include "object.h"
#include "reactionFactory.hh"
#include <vector>
#include <sstream>

#ifdef USE_CUDA
# include "lazy_array.hh"
# include <nvrtc.h>
# include <cuda.h>
#else //USE_CUDA
# if 0
#  include <simdops/resetArch.hpp>
# endif
# include <simdops/simdops.hpp>
# include "VectorDouble32.hh"
#endif //USE_CUDA

REACTION_FACTORY(Grandi)(OBJECT* obj, const double dt, const int numPoints, const ThreadTeam& group);    

namespace Grandi
{

#ifndef USE_CUDA
   struct State
   {

      double CaM[SIMDOPS_FLOAT64V_WIDTH];
      double Cai[SIMDOPS_FLOAT64V_WIDTH];
      double Caj[SIMDOPS_FLOAT64V_WIDTH];
      double Casl[SIMDOPS_FLOAT64V_WIDTH];
      double Casr[SIMDOPS_FLOAT64V_WIDTH];
      double Ki[SIMDOPS_FLOAT64V_WIDTH];
      double Myc[SIMDOPS_FLOAT64V_WIDTH];
      double Mym[SIMDOPS_FLOAT64V_WIDTH];
      double NaBj[SIMDOPS_FLOAT64V_WIDTH];
      double NaBsl[SIMDOPS_FLOAT64V_WIDTH];
      double Nai[SIMDOPS_FLOAT64V_WIDTH];
      double Naj[SIMDOPS_FLOAT64V_WIDTH];
      double Nasl[SIMDOPS_FLOAT64V_WIDTH];
      double RyRi[SIMDOPS_FLOAT64V_WIDTH];
      double RyRo[SIMDOPS_FLOAT64V_WIDTH];
      double RyRr[SIMDOPS_FLOAT64V_WIDTH];
      double SLHj[SIMDOPS_FLOAT64V_WIDTH];
      double SLHsl[SIMDOPS_FLOAT64V_WIDTH];
      double SLLj[SIMDOPS_FLOAT64V_WIDTH];
      double SLLsl[SIMDOPS_FLOAT64V_WIDTH];
      double SRB[SIMDOPS_FLOAT64V_WIDTH];
      double TnCHc[SIMDOPS_FLOAT64V_WIDTH];
      double TnCHm[SIMDOPS_FLOAT64V_WIDTH];
      double TnCL[SIMDOPS_FLOAT64V_WIDTH];
      double d[SIMDOPS_FLOAT64V_WIDTH];
      double f[SIMDOPS_FLOAT64V_WIDTH];
      double fcaBj[SIMDOPS_FLOAT64V_WIDTH];
      double fcaBsl[SIMDOPS_FLOAT64V_WIDTH];
      double h[SIMDOPS_FLOAT64V_WIDTH];
      double hL[SIMDOPS_FLOAT64V_WIDTH];
      double j[SIMDOPS_FLOAT64V_WIDTH];
      double m[SIMDOPS_FLOAT64V_WIDTH];
      double mL[SIMDOPS_FLOAT64V_WIDTH];
      double xkr[SIMDOPS_FLOAT64V_WIDTH];
      double xks[SIMDOPS_FLOAT64V_WIDTH];
      double xkur[SIMDOPS_FLOAT64V_WIDTH];
      double xtf[SIMDOPS_FLOAT64V_WIDTH];
      double ykur[SIMDOPS_FLOAT64V_WIDTH];
      double ytf[SIMDOPS_FLOAT64V_WIDTH];
   };
#endif //USE_CUDA

   class ThisReaction : public Reaction
   {
    public:
      ThisReaction(const int numPoints, const double __dt);
      std::string methodName() const;
      
      void createInterpolants(const double _dt);
      //void updateNonGate(double dt, const VectorDouble32&Vm, VectorDouble32&dVR);
      //void updateGate   (double dt, const VectorDouble32&Vm) ;
      virtual void getCheckpointInfo(std::vector<std::string>& fieldNames,
                                     std::vector<std::string>& fieldUnits) const;
      virtual int getVarHandle(const std::string& varName) const;
      virtual void setValue(int iCell, int varHandle, double value);
      virtual double getValue(int iCell, int varHandle) const;
      virtual double getValue(int iCell, int varHandle, double V) const;
      virtual const std::string getUnit(const std::string& varName) const;

    private:
      unsigned nCells_;
      double __cachedDt;

    public:
      //PARAMETERS
      double AF;
      double ISO;
      double RA;
    public:
      void calc(double dt,
                ro_mgarray_ptr<double> Vm_m,
                ro_mgarray_ptr<double> iStim_m,
                wo_mgarray_ptr<double> dVm_m);
      void initializeMembraneVoltage(wo_mgarray_ptr<double> Vm);
      virtual ~ThisReaction();
#ifdef USE_CUDA
      void constructKernel();

      lazy_array<double> stateTransport_;
      std::string _program_code;
      nvrtcProgram _program;
      std::vector<char> _ptx;
      CUmodule _module;
      CUfunction _kernel;
      int blockSize_;
#else //USE_CUDA

      std::vector<State, AlignedAllocator<State> > state_;
#endif

      //BGQ_HACKFIX, compiler bug with zero length arrays
      Interpolation _interpolant[45+1];
      FRIEND_FACTORY(Grandi)(OBJECT* obj, const double dt, const int numPoints, const ThreadTeam& group);
   };
}


