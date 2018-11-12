#include "Reaction.hh"
#include "Interpolation.hh"
#include "object.h"
#include "reactionFactory.hh"
#include <vector>
#include <sstream>

#ifdef USE_CUDA
# include "TransportCoordinator.hh"
# include <nvrtc.h>
# include <cuda.h>
#else //USE_CUDA
# if 0
#  include <simdops/resetArch.hpp>
# endif
# include <simdops/simdops.hpp>
#endif //USE_CUDA

REACTION_FACTORY(sundnes_et_al_2016_FHN)(OBJECT* obj, const double dt, const int numPoints, const ThreadTeam& group);    

namespace sundnes_et_al_2016_FHN
{

#ifndef USE_CUDA
   struct State
   {

      double W[SIMDOPS_FLOAT64V_WIDTH];
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
      double Vpeak;
      double Vrest;
      double Vthresh;
#ifdef USE_CUDA
    public:
      void calc(double dt,
                const Managed<ArrayView<double>> Vm_m,
                const Managed<ArrayView<double>> iStim_m,
                Managed<ArrayView<double>> dVm_m);
      void initializeMembraneVoltage(ArrayView<double> Vm);
      virtual ~ThisReaction();
      void constructKernel();

      TransportCoordinator<PinnedVector<double> > stateTransport_;
      std::string _program_code;
      nvrtcProgram _program;
      std::vector<char> _ptx;
      CUmodule _module;
      CUfunction _kernel;
      int blockSize_;
#else //USE_CUDA
    public:
      void calc(double dt,
                const VectorDouble32& Vm,
                const std::vector<double>& iStim,
                VectorDouble32& dVm);
      void initializeMembraneVoltage(VectorDouble32& Vm);

      std::vector<State, AlignedAllocator<State> > state_;
#endif

      //BGQ_HACKFIX, compiler bug with zero length arrays
      Interpolation _interpolant[0+1];
      FRIEND_FACTORY(sundnes_et_al_2016_FHN)(OBJECT* obj, const double dt, const int numPoints, const ThreadTeam& group);
   };
}


