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

ReactionFunction::ReactionFunction(QuadratureSpace *qs, FiniteElementSpace *f, const double new_dt, const std::string new_objectName, const ThreadTeam& group)
: QuadratureFunction(qs), QuadS(qs), fes(f), dt(new_dt), objectName(new_objectName), threadGroup(group)
{

      
   nCells = this->Size();

   int vectorLength = nCells;
   Vm.resize(vectorLength);
   dVm.resize(vectorLength);
   iStim.resize(vectorLength);

   double* VmRaw = Vm.writeonly(CPU).raw();
   VmQuad.SetSpace(QuadS, VmRaw, 1);

   std::vector<int> cellTypes;
   for (int i = 0; i < fes->GetNE(); ++i)
   {
      ElementTransformation *T = fes->GetElementTransformation(i);
      //This is a hack.  There's no way to get access to the offsets() array
      //in Quadrature Space without declaring ourselves to be a friend class.
      //This is broken and I hope it is fixed in 4.0
      Vector localVm;
      VmQuad.GetElementValues(i, localVm);
      for ( int j=0; j<localVm.Size(); j++)
      {
         cellTypes.push_back(T->Attribute);
      }
   }
      
   reaction.addReaction(objectName);
   reaction.create(dt, cellTypes, threadGroup);
}

void ReactionFunction::Initialize()
{
   reaction.initializeMembraneState(Vm);
   this->SetData(dVm.readwrite(CPU).raw());
}

void ReactionFunction::CalcVm(const Vector &xs)
{
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
   reaction.calc(dt, Vm, iStim, dVm);
   dVm.readonly(CPU);
}

QuadratureIntegrator::QuadratureIntegrator(QuadratureFunction* p_newQuadFunction, const double new_scale)
: p_quadFunction(p_newQuadFunction),
  scale(new_scale)
{}

void QuadratureIntegrator::AssembleRHSElementVect(const FiniteElement &el, ElementTransformation &Tr, Vector &result)
{
   int dof_u = el.GetDof();
   result.SetSize(dof_u);
   Vector shape(dof_u);
   
   const IntegrationRule &ir = p_quadFunction->GetSpace()->GetElementIntRule(Tr.ElementNo);

   result = 0.0;
   Vector localQuad;
   p_quadFunction->GetElementValues(Tr.ElementNo, localQuad);
   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      el.CalcShape(ip, shape);
      shape *= ip.weight*Tr.Weight()*scale*localQuad[i];
      result += shape;
   }
}

double StimulusCollection::Eval(ElementTransformation& T, const IntegrationPoint &ip)
{
   double x[3];
   Vector transip(x, 3);
 
   T.Transform(ip, transip);

   double result = 0;
   for (std::size_t istim=0; istim<stim_.size(); istim++)
   {
      result += dt_*stim_[istim].eval(time_, T.ElementNo, transip);
   }
   return result;
}

void StimulusCollection::add(Stimulus newStim)
{
   stim_.push_back(newStim);
}
