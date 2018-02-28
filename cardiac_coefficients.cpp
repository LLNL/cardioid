#include "mfem.hpp"
#include "cardiac_physics.hpp"
#include "cardiac_coefficients.hpp"
#include <memory>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
using namespace mfem;

double QuadratureFunctionCoefficient::Eval(ElementTransformation &T,
                                           const IntegrationPoint &ip)
{
   mfem_error ("QuadratureFunctionCoefficient::Eval (...)\n"
               "   is not implemented for this class.");
   return 0.0;

}

double QuadratureFunctionCoefficient::EvalQ(ElementTransformation &T,
                                            const int num_ip)
{
   int elem_no = T.ElementNo;
   return QuadF->GetElementValue(elem_no, num_ip);
}

ActiveTensionFunction::ActiveTensionFunction(QuadratureSpace *qs, FiniteElementSpace *f)
   : QuadratureFunction(qs), QuadS(qs), fes(f), tester(nCells)
{
   nCells = this->Size();

   stretch.resize(nCells);
   stretchVel.resize(nCells);
   tension.resize(nCells);
   dtension.resize(nCells);
   actTime.resize(nCells);
   nextStretch.resize(nCells);

   bool allOK = false;

   outputNames.push_back("tension");
   outputNames.push_back("dtension_dstretchVel");

   allOK = tester.getOrder("output",outputNames,outOrder);
   MFEM_VERIFY(allOK, "Cell model output not set!");;

   outArrays[outOrder[0]] = &tension[0];
   outArrays[outOrder[1]] = &dtension[0];
   
   inputNames.push_back("stretch");
   inputNames.push_back("stretchVel");
   inputNames.push_back("actTime");

   allOK = tester.getOrder("input",inputNames,inOrder);
   MFEM_VERIFY(allOK, "Cell model input not set!");

   inArrays[inOrder[0]] = &stretch[0];
   inArrays[inOrder[1]] = &stretchVel[0];
   inArrays[inOrder[2]] = &actTime[0];

}

void ActiveTensionFunction::Initialize()
{

   for (int icell=0; icell<nCells; icell++)
   {
      stretch[icell] = 1;
      stretchVel[icell] = 1;
      actTime[icell] = 0;
   }
   
   tester.initialize(inArrays);

   //this->SetData(tension.data());   
}

void ActiveTensionFunction::CalcStretch(const Vector &x, const double dt)
{
   for (int icell=0; icell<nCells; icell++) {
      nextStretch[icell] = 1;
      stretchVel[icell] = (nextStretch[icell]-stretch[icell])/dt;
   }

}

void ActiveTensionFunction::TryStep(const Vector &x, const double dt)
{
   CalcStretch(x, dt);
   tester.tryTimestep(dt, inArrays);
   tester.outputForNextTimestep(dt, inArrays, outArrays);
   tester.commitTimestep();
}

void ActiveTensionFunction::CommitStep(const double dt)
{
   tester.outputForNextTimestep(dt, inArrays, outArrays);
   tester.commitTimestep();
   swap(stretch,nextStretch);
   
}

ActiveTensionFunction & ActiveTensionFunction::operator=(double value)
{
   Vector::operator=(value);
   return *this;
}
