#ifndef TEST_REACTION_HH
#define TEST_REACTION_HH

#include "Reaction.hh"

struct TestReactionParms
{
   double initialVoltage;
   double delta;
};

/** Note that this class does not have many of the features we might
 *  want in a fully featured Test class.  Those would include (for
 *  example):
 *
 *  - ability to set the random seed
 *  - construction of random numbers such that the initial voltages are
 *    the same regardless of the number of tasks on which the problem is
 *    run.
 *
 *  We know how to do such things, I just haven't taken the time in this
 *  case.
 */
class TestReaction : public Reaction
{
 public:
   
   TestReaction(const TestReactionParms& parms)
   : V0_(parms.initialVoltage),
     delta_(parms.delta)
   {};

   std::string methodName() const {return "test";}

   void calc(double dt,
             const Managed<ArrayView<double>> Vm,
             const Managed<ArrayView<double>> iStim,
             Managed<ArrayView<double>> dVm) {}
   void initializeMembraneVoltage(ArrayView<double> Vm)
   {
      for (unsigned ii=0; ii< Vm.size(); ++ii)
         Vm[ii] = V0_ + delta_ * (2*drand48() - 1.0);
   };

   /** Functions needed for checkpoint/restart */
   void getCheckpointInfo(std::vector<std::string>& fieldNames,
                          std::vector<std::string>& fieldUnits) const{};
   int getVarHandle(const std::string& varName) const {return -1;};
   void setValue(int iCell, int varHandle, double value){};
   double getValue(int iCell, int varHandle) const{return 0.0;};
   void getValue(int iCell,
                 const std::vector<int>& handle,
                 std::vector<double>& value) const{};
   const std::string getUnit(const std::string& varName) const{return "1";};


 private:

   double V0_;
   double delta_;

};

#endif
