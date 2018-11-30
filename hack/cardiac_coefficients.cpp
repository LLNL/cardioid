#include "mfem.hpp"
#include "cardiac_physics.hpp"
#include "cardiac_coefficients.hpp"
#include <memory>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
using namespace mfem;

/// Standard coefficient evaluation is not valid
double QuadratureFunctionCoefficient::Eval(ElementTransformation &T,
                                           const IntegrationPoint &ip)
{
   mfem_error ("QuadratureFunctionCoefficient::Eval (...)\n"
               "   is not implemented for this class.");
   return 0.0;

}

/// Evaluate the function coefficient at a specific quadrature point
double QuadratureFunctionCoefficient::EvalQ(ElementTransformation &T,
                                            const int num_ip)
{
   int elem_no = T.ElementNo;
   return QuadF->GetElementValue(elem_no, num_ip);
}

ActiveTensionFunction::ActiveTensionFunction(QuadratureSpace *qs, FiniteElementSpace *f, VectorCoefficient &fib)
   : QuadratureFunction(qs), QuadS(qs), fes(f), tester(qs->GetSize()), Q(&fib)
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
      stretch[icell] = 1.0;
      stretchVel[icell] = 0.0;
      actTime[icell] = 0.0;
   }
   
   tester.initialize(inArrays);

   this->SetData(tension.data());   
}

void ActiveTensionFunction::CalcStretch(const Vector &x, const double dt)
{
   Vector xs_true(x.GetData(), fes->GetTrueVSize());
   Vector xs(fes->GetVSize());
   fes->GetProlongationMatrix()->Mult(xs_true, xs);

   ElementTransformation *T;
   const FiniteElement *fe;
   Vector el_x;
   Array<int> vdofs;
   Vector fib(3);
   Vector fib_out(3);
   
   DenseMatrix J0i(3);
   DenseMatrix DSh_u;
   DenseMatrix DS_u;
   DenseMatrix J(3);
   DenseMatrix PMatI_u;

   int dof;
   
   for (int i = 0; i < fes->GetNE(); ++i) {
      T = fes->GetElementTransformation(i);
      fes->GetElementVDofs(i, vdofs);
      fe = fes->GetFE(i);
      xs.GetSubVector(vdofs, el_x);

      dof = fe->GetDof();
      PMatI_u.UseExternalData(el_x.GetData(), dof, 3);
      DSh_u.SetSize(dof, 3);
      DS_u.SetSize(dof, 3);
            
      int intorder = 2*fe->GetOrder() + 3; 
      const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(), intorder);
      
      for (int i_num=0; i_num<ir.GetNPoints(); i_num++) {

         const IntegrationPoint &ip = ir.IntPoint(i_num);
         T->SetIntPoint(&ip);
         CalcInverse(T->Jacobian(), J0i);

         fe->CalcDShape(ip, DSh_u);
         Mult(DSh_u, J0i, DS_u);
         MultAtB(PMatI_u, DS_u, J);
         Q->Eval(fib, *T, ip);
         fib /= fib.Norml2();
         J.Mult(fib, fib_out);         
         nextStretch[(this)->GetElementOffset(T->ElementNo) + i_num] =  fib_out.Norml2();                  
      }
   }
      
   for (int icell=0; icell<nCells; icell++) {
      stretchVel[icell] = (nextStretch[icell]-stretch[icell])/dt;
   }

}

void ActiveTensionFunction::TryStep(const Vector &x, const double dt)
{
   CalcStretch(x, dt);
   tester.tryTimestep(dt, inArrays);
   tester.outputForNextTimestep(dt, inArrays, outArrays);
}

void ActiveTensionFunction::CommitStep(const double dt)
{
   tester.commitTimestep();
   swap(stretch,nextStretch);
   
}

ActiveTensionFunction & ActiveTensionFunction::operator=(double value)
{
   Vector::operator=(value);
   return *this;
}
