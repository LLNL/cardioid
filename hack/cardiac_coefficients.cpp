#include "mfem.hpp"
#include "cardiac_coefficients.hpp"
#include "reactionFactory.hh"
#include "ThreadServer.hh"
#include "object.h"
#include "object_cc.hh"
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
   //int elem_no = T.ElementNo;
   //return QuadF->GetElementValue(elem_no, num_ip);
}

ReactionFunction::ReactionFunction(QuadratureSpace *qs, FiniteElementSpace *f, const double new_dt, const std::string new_objectName, const ThreadTeam& group)
: QuadratureFunction(qs), QuadS(qs), fes(f), dt(new_dt), objectName(new_objectName), threadGroup(group)
{
   nCells = this->Size();

   Vm.resize(nCells);
   dVm.resize(nCells);
   iStim.resize(nCells);

   double* VmRaw = Vm.writeonly(CPU).raw();
   VmQuad.SetSpace(QuadS, VmRaw, 1);
   
   reaction = std::shared_ptr<Reaction>(reactionFactory(objectName, dt, nCells, threadGroup));
}

void ReactionFunction::Initialize()
{
   initializeMembraneState(reaction.get(), objectName, Vm);
   this->SetData(dVm.readwrite(CPU).raw());
}

void ReactionFunction::CalcVm(const Vector &x)
{
   Vector xs_true(x.GetData(), fes->GetTrueVSize());
   Vector xs(fes->GetVSize());
   fes->GetProlongationMatrix()->Mult(xs_true, xs);

   ElementTransformation *T;
   const FiniteElement *fe;
   Vector el_x;
   Array<int> vdofs;
   Vector basisFunctions;
   
   int dof;

   wo_array_ptr<double> VmRaw = Vm.writeonly(CPU);
   for (int i = 0; i < fes->GetNE(); ++i) {
      T = fes->GetElementTransformation(i);
      fes->GetElementVDofs(i, vdofs);
      fe = fes->GetFE(i);
      xs.GetSubVector(vdofs, el_x);

      dof = fe->GetDof();
      basisFunctions.SetSize(dof);

      const IntegrationRule &ir = QuadS->GetElementIntRule(i);
      Vector localVm;
      VmQuad.GetElementValues(i, localVm);
      for (int i_num=0; i_num<ir.GetNPoints(); i_num++) {
         const IntegrationPoint &ip = ir.IntPoint(i_num);
         T->SetIntPoint(&ip);
         fe->CalcShape(ip, basisFunctions);
         localVm[i_num] = el_x*basisFunctions;
      }
   }      
}

void ReactionFunction::Calc(const Vector &x)
{
   CalcVm(x);
   reaction->calc(dt, Vm, iStim, dVm);
   dVm.readonly(CPU);
}
