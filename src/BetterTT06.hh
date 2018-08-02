#include "Reaction.hh"
#include "Interpolation.hh"
#include "TransportCoordinator.hh"
#include "object.h"
#include "reactionFactory.hh"
#include <vector>
#include <sstream>
#include <nvrtc.h>
#include <cuda.h>

REACTION_FACTORY(BetterTT06 )(OBJECT* obj, const double dt, const int numPoints, const ThreadTeam& group);

namespace BetterTT06
{
   class ThisReaction : public Reaction
   {
    public:
      ThisReaction(const int numPoints, const double __dt);
      std::string methodName() const;
      
      void createInterpolants(const double _dt);
      void calc(double dt,
                const Managed<ArrayView<double>> Vm_m,
                const Managed<ArrayView<double>> iStim_m,
                Managed<ArrayView<double>> dVm_m);
      //void updateNonGate(double dt, const VectorDouble32&Vm, VectorDouble32&dVR);
      //void updateGate   (double dt, const VectorDouble32&Vm) ;
      void initializeMembraneVoltage(ArrayView<double> Vm);
      virtual void getCheckpointInfo(std::vector<std::string>& fieldNames,
                                     std::vector<std::string>& fieldUnits) const;
      virtual int getVarHandle(const std::string& varName) const;
      virtual void setValue(int iCell, int varHandle, double value);
      virtual double getValue(int iCell, int varHandle) const;
      virtual double getValue(int iCell, int varHandle, double V) const;
      virtual const std::string getUnit(const std::string& varName) const;

      virtual ~ThisReaction();
    public:

    private:
      void constructKernel();

      //PARAMETERS
      unsigned nCells_;
      TransportCoordinator<PinnedVector<double> > stateTransport_;
      double __cachedDt;
      std::string _program_code;
      nvrtcProgram _program;
      std::vector<char> _ptx;
      CUmodule _module;
      CUfunction _kernel;
      int blockSize_;

      Interpolation _interpolant[31];
      FRIEND_FACTORY(BetterTT06)(OBJECT* obj, const double dt, const int numPoints, const ThreadTeam& group);
   };
}


