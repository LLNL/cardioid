#include "Reaction.hh"
#include "Interpolation.hh"
#include "TransportCoordinator.hh"
#include "object.h"
#include <vector>
#include <sstream>
#include <nvrtc.h>
#include <cuda.h>

namespace scanReaction 
{
    Reaction* scanBetterTT06(OBJECT* obj, const int numPoints, const double __dt);
}

namespace BetterTT06
{
   class ThisReaction : public Reaction
   {
    public:
      ThisReaction(const int numPoints, const double __dt);
      std::string methodName() const;
      
      void createInterpolants(const double _dt);
      void calc(double dt,
                const VectorDouble32& Vm,
                const std::vector<double>& iStim,
                VectorDouble32& dVm);
      //void updateNonGate(double dt, const VectorDouble32&Vm, VectorDouble32&dVR);
      //void updateGate   (double dt, const VectorDouble32&Vm) ;
      void initializeMembraneVoltage(VectorDouble32& Vm);
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
      TransportCoordinator<std::vector<double> > stateTransport_;
      double __cachedDt;
      std::string _program_code;
      nvrtcProgram _program;
      std::vector<char> _ptx;
      CUmodule _module;
      CUfunction _kernel;
      int blockSize_;

      Interpolation _interpolant[31];
      friend Reaction* scanReaction::scanBetterTT06(OBJECT* obj, const int numPoints, const double __dt);
   };
}


