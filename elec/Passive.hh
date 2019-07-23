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

REACTION_FACTORY(Passive)(OBJECT* obj, const double dt, const int numPoints, const ThreadTeam& group);    

namespace Passive
{

#ifndef USE_CUDA
   struct State
   {

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
      double E_R;
      double G;
    public:
      void calc(double dt,
                ro_mgarray_ptr<int> indexArray,
                ro_mgarray_ptr<double> Vm_m,
                ro_mgarray_ptr<double> iStim_m,
                wo_mgarray_ptr<double> dVm_m);
      void initializeMembraneVoltage(ro_mgarray_ptr<int> indexArray, wo_mgarray_ptr<double> Vm);
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
      Interpolation _interpolant[0+1];
      FRIEND_FACTORY(Passive)(OBJECT* obj, const double dt, const int numPoints, const ThreadTeam& group);
   };
}


