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

ReactionWrapper::ReactionWrapper(const double new_dt, const std::vector<std::string> new_objectNames, const ThreadTeam& group, const std::vector<int>& cellTypes)
: dt(new_dt), threadGroup(group)
{
      
   nCells = cellTypes.size();

   Vm.resize(nCells);
   dVm.resize(nCells);
   iStim.resize(nCells);

   for (auto objectName : new_objectNames) {
      reaction.addReaction(objectName);
   }
   reaction.create(dt, cellTypes, threadGroup);

   rw_array_ptr<double> Vm_ptr = Vm.readwrite(CPU);
   Vm_vector.SetDataAndSize(Vm_ptr.raw(), nCells);
   rw_array_ptr<double> Iion_ptr = dVm.readwrite(CPU);
   Iion_vector.SetDataAndSize(Iion_ptr.raw(), nCells);
}

void ReactionWrapper::Initialize()
{
   reaction.initializeMembraneState(Vm);
}

void ReactionWrapper::Calc()
{
   reaction.calc(dt, Vm, iStim, dVm);
}

Vector& ReactionWrapper::getVmReadwrite()
{
   rw_array_ptr<double> Vm_ptr = Vm.useOn(CPU);
   return Vm_vector;
}

Vector& ReactionWrapper::getIionReadwrite()
{
   rw_array_ptr<double> Iion_ptr = dVm.useOn(CPU);
   return Iion_vector;
}

const Vector& ReactionWrapper::getVmReadonly() const
{
   ro_array_ptr<double> Vm_ptr = Vm.useOn(CPU);
   return Vm_vector;
}

const Vector& ReactionWrapper::getIionReadonly() const
{
   ro_array_ptr<double> Iion_ptr = dVm.useOn(CPU);
   return Iion_vector;
}

ReactionFunction::ReactionFunction(QuadratureSpace *qs, FiniteElementSpace *f,ReactionWrapper* rw)
: QuadratureFunction(qs),reactionWrapper(rw), QuadS(qs), fes(f)
{
   double* VmRaw = rw->getVmReadwrite().GetData();
   VmQuad.SetSpace(QuadS, VmRaw, 1);
}

void ReactionFunction::Calc(const Vector &xs)
{
   ElementTransformation *T;
   const FiniteElement *fe;
   Vector el_x;
   Array<int> vdofs;
   Vector basisFunctions;
   
   int dof;

   wo_array_ptr<double> VmRaw = reactionWrapper->Vm.writeonly(CPU);
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
   reactionWrapper->Calc();
   reactionWrapper->getIionReadonly();
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
