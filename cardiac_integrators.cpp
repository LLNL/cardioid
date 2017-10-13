#include "mfem.hpp"
#include "cardiac_integrators.hpp"

namespace mfem
{

void CardiacModel::EvalP(const DenseMatrix &J, const double pres, const Vector &fiber, DenseMatrix &P) const
{
   int dim = J.Width();
   double dJ = J.Det();

   E.SetSize(dim);
   C.SetSize(dim);
   JT.SetSize(dim);
   I.SetSize(dim);
   P.SetSize(dim);
   PK2.SetSize(dim);
   FinvT.SetSize(dim);

   I = 0.0;
   for (int d=0; d<dim; d++) {
      I(d,d) = 1.0;
   }

   JT = J;
   JT.Transpose();

   Mult(JT, J, C);

   CalcInverseTranspose(J, FinvT);

   Add(0.5, C, -0.5, I, E);

   // Determine orthonormal fiber coordinate system
   Vector fib = fiber;
   fib /= fib.Norml2();

   Vector orth1(3);
   Vector orth2(3);
   
   Vector test(3);
   test = 0.0;
   test(0) = 1.0;
   
   cross(fib, test, orth1);

   if (orth1.Norml2() < 1.0e-8) {
      test = 0.0;
      test(1) = 1.0;
      cross(fib, test, orth1);
   }
   
   orth1 /= orth1.Norml2();
   cross(fib, orth1, orth2);
   orth2 /= orth2.Norml2();

   // Compute transformation tensor to fiber coordinates
   DenseMatrix orth(3); 

   for (int i=0; i<3; i++) {
      orth(i,0) = fib(i);
      orth(i,1) = orth1(i);
      orth(i,2) = orth2(i);
   }

   DenseMatrix orth_transpose(3);
   orth_transpose.Transpose(orth);

   //Transform E into fiber coordinates
   // E' = Q^T E Q
   DenseMatrix dummy(3);
   Mult(orth_transpose, E, dummy);
   Mult(dummy, orth, E);

   // Apply cardiac constitutive model
   double Q = b_ff * E(0,0) * E(0,0) +
      b_ss * E(1,1) * E(1,1) +
      b_nn * E(2,2) * E(2,2) +
      b_fs * (E(0,1) * E(0,1) + E(1,0) * E(1,0)) +
      b_fn * (E(0,2) * E(0,2) + E(2,0) * E(2,0)) +
      b_ns * (E(1,2) * E(1,2) + E(2,1) * E(2,1));
      
   double fact = (C_1/2.0) * exp(Q);

   PK2(0,0) = fact * b_ff * E(0,0);
   PK2(1,1) = fact * b_ss * E(1,1);
   PK2(2,2) = fact * b_nn * E(2,2);

   PK2(0,1) = fact * b_fs * E(0,1);
   PK2(1,0) = fact * b_fs * E(1,0);

   PK2(0,2) = fact * b_fn * E(0,2);
   PK2(2,0) = fact * b_fn * E(2,0);

   PK2(1,2) = fact * b_ns * E(1,2);
   PK2(2,1) = fact * b_ns * E(2,1);

   // Transform back from fiber coordinates
   // PK2 = Q PK2' Q^T 
   Mult(orth, PK2, dummy);
   Mult(dummy, orth_transpose, PK2);

   Mult(J, PK2, P); 
   P.Add(-1.0*pres*dJ, FinvT);
 
}

double CardiacModel::EvalW(const DenseMatrix &J) const
{
   mfem_error("CardiacModel::EvalW"
              " is not overloaded!");

   return 0.0;

}

void CardiacModel::EvalP(const DenseMatrix &J, DenseMatrix &P) const
{
   mfem_error("CardiacModel::EvalP"
              " is not overloaded!");
}

void CardiacModel::AssembleH(
   const DenseMatrix &J, const DenseMatrix &DS, const double weight,
   DenseMatrix &A) const
{
   mfem_error("CardiacModel::AssembleH"
              " is not overloaded!");
}

void CardiacNLFIntegrator::AssembleElementVector(Array<const FiniteElement *> &el,
                                                 ElementTransformation &Tr,
                                                 Array<Vector *> &elfun, 
                                                 Array<Vector *> &elvec)
{
   int dof_u = el[0]->GetDof();
   int dim_u = el[0]->GetDim();

   int dof_p = el[1]->GetDof();
   int dim_p = el[1]->GetDim();
      
   if (dim_u != dim_p) {
      mfem_error("CardiacNLFIntegrator::AssembleElementVector"
                 " dimensions of FE spaces not compatible");
   }

   DSh_u.SetSize(dof_u, dim_u);
   DS_u.SetSize(dof_u, dim_u);
   J0i.SetSize(dim_u);
   J.SetSize(dim_u);
   P.SetSize(dim_u);
   PMatI_u.UseExternalData(elfun[0]->GetData(), dof_u, dim_u);
   elvec[0]->SetSize(dof_u*dim_u);
   PMatO_u.UseExternalData(elvec[0]->GetData(), dof_u, dim_u);

   Sh_p.SetSize(dof_p);
   elvec[1]->SetSize(dof_p);

   int intorder = 2*el[0]->GetOrder() + 3; 
   const IntegrationRule &ir = IntRules.Get(el[0]->GetGeomType(), intorder);

   *elvec[0] = 0.0;
   *elvec[1] = 0.0;

   double pres;

   Vector fiber;

   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      Tr.SetIntPoint(&ip);
      CalcInverse(Tr.Jacobian(), J0i);

      el[0]->CalcDShape(ip, DSh_u);
      Mult(DSh_u, J0i, DS_u);
      MultAtB(PMatI_u, DS_u, J);

      el[1]->CalcShape(ip, Sh_p);
      pres = Sh_p * *elfun[1];

      Q->Eval(fiber, Tr, ip);
      model->EvalP(J, pres, fiber, P);

      P *= ip.weight*Tr.Weight();
      AddMultABt(DS_u, P, PMatO_u);

      elvec[1]->Add(ip.weight * Tr.Weight() * (J.Det() - 1.0), Sh_p);

   }

}

CardiacNLFIntegrator::~CardiacNLFIntegrator()
{ }

void BodyForceNLFIntegrator::AssembleElementVector(Array<const FiniteElement *> &el,
                                                         ElementTransformation &Tr,
                                                         Array<Vector *> &elfun, 
                                                         Array<Vector *> &elvec)
{
   int dof_u = el[0]->GetDof();
   int dim_u = el[0]->GetDim();   

   shape.SetSize(dof_u);
   elvec[0]->SetSize(dof_u*dim_u);

   elvec[1]->SetSize(el[1]->GetDof());


   int intorder = 2*el[0]->GetOrder() + 3; 
   const IntegrationRule &ir = IntRules.Get(el[0]->GetGeomType(), intorder);

   *elvec[0] = 0.0;
   *elvec[1] = 0.0;

   double val, cf;

   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);

      Tr.SetIntPoint(&ip);
      val = Tr.Weight();

      el[0]->CalcShape(ip, shape);
      function.Eval(Qvec, Tr, ip);

      for (int k = 0; k < dim_u; k++) {
         cf = -1.0 * val * Qvec(k);

         for (int s = 0; s < dof_u; s++) {
            (*elvec[0])(dof_u*k + s) += ip.weight * cf * shape(s);
         }
      }
   }

}

BodyForceNLFIntegrator::~BodyForceNLFIntegrator()
{ }

void ActiveTensionNLFIntegrator::AssembleElementVector(Array<const FiniteElement *> &el,
                                                         ElementTransformation &Tr,
                                                         Array<Vector *> &elfun, 
                                                         Array<Vector *> &elvec)
{
   int dof_u = el[0]->GetDof();
   int dim_u = el[0]->GetDim();

   DSh_u.SetSize(dof_u, dim_u);
   DS_u.SetSize(dof_u, dim_u);
   J0i.SetSize(dim_u);
   J.SetSize(dim_u);
   P.SetSize(dim_u);
   AStress.SetSize(dim_u);
   AStress_current.SetSize(dim_u);
   PMatI_u.UseExternalData(elfun[0]->GetData(), dof_u, dim_u);
   elvec[0]->SetSize(dof_u*dim_u);
   PMatO_u.UseExternalData(elvec[0]->GetData(), dof_u, dim_u);

   int intorder = 2*el[0]->GetOrder() + 3; 
   const IntegrationRule &ir = IntRules.Get(el[0]->GetGeomType(), intorder);

   elvec[1]->SetSize(el[1]->GetDof());

   *(elvec[0]) = 0.0;
   *(elvec[1]) = 0.0;

   for (int i = 0; i < ir.GetNPoints(); i++)
   {

      const IntegrationPoint &ip = ir.IntPoint(i);
      Tr.SetIntPoint(&ip);
      CalcInverse(Tr.Jacobian(), J0i);

      el[0]->CalcDShape(ip, DSh_u);
      Mult(DSh_u, J0i, DS_u);
      MultAtB(PMatI_u, DS_u, J);

      function.Eval(AStress, Tr, ip);
      Mult(J, AStress, AStress_current);

      AStress_current *= ip.weight*Tr.Weight();
      AddMultABt(DS_u, AStress_current, PMatO_u);
   }
}

ActiveTensionNLFIntegrator::~ActiveTensionNLFIntegrator()
{ }


void PressureBoundaryNLFIntegrator::AssembleRHSElementVector(Array<const FiniteElement *> &el,
                                                             FaceElementTransformations &Tr,
                                                             Array<Vector *> &elfun, 
                                                             Array<Vector *> &elvec)
{
   int dim = el[0]->GetDim();
   int dof_u = el[0]->GetDof();
   int dof_p = el[1]->GetDof();

   shape.SetSize (dof_u);
   nor.SetSize (dim);
   fnor.SetSize (dim);
   elvec[0]->SetSize (dim*dof_u);
   elvec[1]->SetSize (dof_p);

   DSh_u.SetSize(dof_u, dim);
   DS_u.SetSize(dof_u, dim);
   J0i.SetSize(dim);
   J.SetSize(dim);
   Jinv.SetSize(dim);

   PMatI_u.UseExternalData(elfun[0]->GetData(), dof_u, dim);

   int intorder = 2*el[0]->GetOrder() + 3; 
   const IntegrationRule &ir = IntRules.Get(Tr.FaceGeom, intorder);

   *(elvec[0]) = 0.0;
   *(elvec[1]) = 0.0;

   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      IntegrationPoint eip;
      Tr.Loc1.Transform(ip, eip);

      Tr.Face->SetIntPoint(&ip);
      CalcOrtho(Tr.Face->Jacobian(), nor);

      //Normalize vector
      double norm = nor.Norml2();
      nor /= norm;

      Tr.Elem1->SetIntPoint(&eip);
      CalcInverse(Tr.Elem1->Jacobian(), J0i);

      el[0]->CalcDShape(eip, DSh_u);
      Mult(DSh_u, J0i, DS_u);
      MultAtB(PMatI_u, DS_u, J);

      CalcInverse(J, Jinv);
      
      Jinv.MultTranspose(nor, fnor);

      el[0]->CalcShape (eip, shape);
      fnor *= -1.0 * function.Eval(*Tr.Face, ip) * J.Det(); 
      for (int j = 0; j < dof_u; j++)
         for (int k = 0; k < dim; k++)
         {
            (*elvec[0])(dof_u*k+j) += ip.weight * Tr.Face->Weight() * fnor(k) * shape(j);
         }
   }
}

void PressureBoundaryNLFIntegrator::AssembleElementVector(Array<const FiniteElement *> &el,
                                                          ElementTransformation &Tr,
                                                          Array<Vector *> &elfun, 
                                                          Array<Vector *> &elvec)
{
   mfem_error("PressureBoundaryNLFIntegrator::AssembleElementVector"
              " is not overloaded!");
}


PressureBoundaryNLFIntegrator::~PressureBoundaryNLFIntegrator()
{ }

void TractionBoundaryNLFIntegrator::AssembleRHSElementVector(Array<const FiniteElement *> &el,
                                                             FaceElementTransformations &Tr,
                                                             Array<Vector *> &elfun, 
                                                             Array<Vector *> &elvec)
{
   int dim = el[0]->GetDim();
   int dof_u = el[0]->GetDof();
   int dof_p = el[1]->GetDof();

   shape.SetSize (dof_u);
   nor.SetSize (dim);
   fnor.SetSize (dim);
   elvec[0]->SetSize (dim*dof_u);
   elvec[1]->SetSize (dof_p);

   DSh_u.SetSize(dof_u, dim);
   DS_u.SetSize(dof_u, dim);
   J0i.SetSize(dim);
   J.SetSize(dim);
   Jinv.SetSize(dim);

   PMatI_u.UseExternalData(elfun[0]->GetData(), dof_u, dim);

   int intorder = 2*el[0]->GetOrder() + 3; 
   const IntegrationRule &ir = IntRules.Get(Tr.FaceGeom, intorder);

   *(elvec[0]) = 0.0;
   *(elvec[1]) = 0.0;

   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      IntegrationPoint eip;
      Tr.Loc1.Transform(ip, eip);

      Tr.Face->SetIntPoint(&ip);

      Tr.Elem1->SetIntPoint(&eip);
      CalcInverse(Tr.Elem1->Jacobian(), J0i);

      el[0]->CalcDShape(eip, DSh_u);
      Mult(DSh_u, J0i, DS_u);
      MultAtB(PMatI_u, DS_u, J);

      CalcInverse(J, Jinv);

      function.Eval(nor, *Tr.Face, ip);
      Jinv.MultTranspose(nor, fnor);

      el[0]->CalcShape (eip, shape);
      fnor *= J.Det(); 
      for (int j = 0; j < dof_u; j++)
         for (int k = 0; k < dim; k++)
         {
            (*elvec[0])(dof_u*k+j) += -1.0 * ip.weight * Tr.Face->Weight() * fnor(k) * shape(j);
         }
   }
}

void TractionBoundaryNLFIntegrator::AssembleElementVector(Array<const FiniteElement *> &el,
                                                          ElementTransformation &Tr,
                                                          Array<Vector *> &elfun, 
                                                          Array<Vector *> &elvec)

{
   mfem_error("TractionBoundaryNLFIntegrator::AssembleElementVector"
              " is not overloaded!");
}


TractionBoundaryNLFIntegrator::~TractionBoundaryNLFIntegrator()
{ }

}
