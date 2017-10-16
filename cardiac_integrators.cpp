#include "mfem.hpp"
#include "cardiac_integrators.hpp"

namespace mfem
{

void CardiacModel::GenerateTransform(const Vector &fiber, DenseMatrix &Q, DenseMatrix &QT) const
{
   Q.SetSize(3);
   QT.SetSize(3);
   
  // Determine orthonormal fiber coordinate system
   Vector fib = fiber;
   fib /= fib.Norml2();

   Vector orth1(3);
   Vector orth2(3);
   
   Vector test(3);
   test = 0.0;
   test(2) = -1.0;
   
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
      Q(i,0) = fib(i);
      Q(i,1) = orth1(i);
      Q(i,2) = orth2(i);
   }

   QT.Transpose(Q);
}

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

   GenerateTransform(fiber, orth, orth_transpose);

   // Transform E into fiber coordinates
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

void CardiacModel::AssembleH(const DenseMatrix &J, const double pres, const Vector &fiber, const DenseMatrix &DS_u, const Vector Sh_p, const double weight, Array2D<DenseMatrix *> &elmats) const
{

   int dof_u = DS_u.Height();
   int dim = DS_u.Width();
   int dof_p = Sh_p.Size();
   
   B.SetSize(3);
   C.SetSize(3);
   E.SetSize(3);
   JT.SetSize(3);

   B(0,0) = b_ff;
   B(0,1) = b_fs;
   B(0,2) = b_fn;
   B(1,0) = b_fs;
   B(1,1) = b_ss;
   B(1,2) = b_ns;
   B(2,0) = b_fn;
   B(2,1) = b_ns;
   B(2,2) = b_nn;

   CalcInverseTranspose(J, FinvT);
   JT = J;
   JT.Transpose();
   
   GenerateTransform(fiber, orth, orth_transpose);
   
   Mult(JT, J, C);
   
   Add(0.5, C, -0.5, I, E);

   // Transform E into fiber coordinates
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
   double dJ = J.Det();

   // u,u block
   for (int i_u = 0; i_u < dof_u; i_u++) {
   for (int i_dim = 0; i_dim < dim; i_dim++) {
      for (int j_u = 0; j_u < dof_u; j_u++) {
      for (int j_dim = 0; j_dim < dim; j_dim++) {

         // m = j_dim;
         // k = i_dim;
         for (int n=0; n<dim; n++) {
         for (int l=0; l<dim; l++) {

            for (int i=0; i<dim; i++) {
               if (i_dim==j_dim && i==n) {
                  (*elmats(0,0))(i_u + i_dim*dof_u, j_u + j_dim*dof_u) += fact * B(i,l) * E(i,l) * DS_u(i_u,l) * DS_u(j_u,n) * weight;
               }

               for (int a=0; a<dim; a++) {
                  (*elmats(0,0))(i_u + i_dim*dof_u, j_u + j_dim*dof_u) += fact * J(i_dim,i) * B(n,a) * E(n,a) * J(j_dim,a) * B(i,l) * E(i,l) * DS_u(i_u,l) * DS_u(j_u,n) * weight;                     
                  (*elmats(0,0))(i_u + i_dim*dof_u, j_u + j_dim*dof_u) += fact * J(i_dim,i) * B(a,n) * E(a,n) * J(j_dim,a) * B(i,l) * E(i,l) * DS_u(i_u,l) * DS_u(j_u,n) * weight;                     
               }
                  
               if (i==n) {
                  (*elmats(0,0))(i_u + i_dim*dof_u, j_u + j_dim*dof_u) += fact * J(i_dim,i) * B(i,l) * 0.5 * J(j_dim,l) * DS_u(i_u,l) * DS_u(j_u,n) * weight;                     
               }
               if (l==n) {
                  (*elmats(0,0))(i_u + i_dim*dof_u, j_u + j_dim*dof_u) += fact * J(i_dim,i) * B(i,l) * 0.5 * J(j_dim,i) * DS_u(i_u,l) * DS_u(j_u,n) * weight;                     
               }
            
                  
            }
               
            (*elmats(0,0))(i_u + i_dim*dof_u, j_u + j_dim*dof_u) -= dJ * pres * FinvT(i_dim,l) * FinvT(j_dim,n) * DS_u(i_u,l) * DS_u(j_u,n) * weight;
               
            // a = n;
            // b = m;
            (*elmats(0,0))(i_u + i_dim*dof_u, j_u + j_dim*dof_u) += dJ * pres * FinvT(i_dim,n) * FinvT(j_dim,l) * DS_u(i_u,l) * DS_u(j_u,n) * weight;
            
         }
         }
      }
      }
   }
   }
   // u,p and p,u blocks
   for (int i_p = 0; i_p < dof_p; i_p++) {
      for (int j_u = 0; j_u < dof_u; j_u++) {
         for (int dim_u = 0; dim_u < dim; dim_u++) {
            for (int l=0; l<dim; l++) {
               (*elmats(1,0))(i_p, j_u + dof_u * dim_u) += dJ * FinvT(dim_u,l) * DS_u(j_u,l) * Sh_p(i_p) * weight; 
               (*elmats(0,1))(j_u + dof_u * dim_u, i_p) -= dJ * FinvT(dim_u,l) * DS_u(j_u,l) * Sh_p(i_p) * weight;
            }               
         }
      }
   }   
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

void CardiacNLFIntegrator::AssembleElementGrad(Array<const FiniteElement*> &el,
                                               ElementTransformation &Tr,
                                               Array<Vector *> &elfun, 
                                               Array2D<DenseMatrix *> &elmats)
{
   int dof_u = el[0]->GetDof();
   int dof_p = el[1]->GetDof();

   int dim = el[0]->GetDim();

   elmats(0,0)->SetSize(dof_u*dim, dof_u*dim);
   elmats(0,1)->SetSize(dof_u*dim, dof_p);
   elmats(1,0)->SetSize(dof_p, dof_u*dim);
   elmats(1,1)->SetSize(dof_p, dof_p);

   *elmats(0,0) = 0.0;
   *elmats(0,1) = 0.0;
   *elmats(1,0) = 0.0;
   *elmats(1,1) = 0.0;

   DSh_u.SetSize(dof_u, dim);
   DS_u.SetSize(dof_u, dim);
   J0i.SetSize(dim);
   F.SetSize(dim);
   FinvT.SetSize(dim);
   Finv.SetSize(dim);
   P.SetSize(dim);
   PMatI_u.UseExternalData(elfun[0]->GetData(), dof_u, dim);
   Sh_p.SetSize(dof_p);
   B.SetSize(dim);

   int intorder = 2*el[0]->GetOrder() + 3; // <---
   const IntegrationRule &ir = IntRules.Get(el[0]->GetGeomType(), intorder);

   Vector fiber(3);

   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      Tr.SetIntPoint(&ip);
      CalcInverse(Tr.Jacobian(), J0i);

      el[0]->CalcDShape(ip, DSh_u);
      Mult(DSh_u, J0i, DS_u);
      MultAtB(PMatI_u, DS_u, F);

      el[1]->CalcShape(ip, Sh_p);
      double pres = Sh_p * *elfun[1];

      Q->Eval(fiber, Tr, ip);

      model->AssembleH(F, pres, fiber, DS_u, Sh_p, ip.weight * Tr.Weight(), elmats);

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
      fnor *= function.Eval(*Tr.Face, ip) * J.Det(); 
      for (int j = 0; j < dof_u; j++)
         for (int k = 0; k < dim; k++)
         {
            (*elvec[0])(dof_u*k+j) += ip.weight * Tr.Face->Weight() * fnor(k) * shape(j);
         }
   }
}

void PressureBoundaryNLFIntegrator::AssembleRHSElementGrad(Array<const FiniteElement*> &el,
                                                           FaceElementTransformations &Tr,
                                                           Array<Vector *> &elfun, 
                                                           Array2D<DenseMatrix *> &elmats)
{
   int dof_u = el[0]->GetDof();
   int dof_p = el[1]->GetDof();

   int dim = el[0]->GetDim();

   elmats(0,0)->SetSize(dof_u*dim, dof_u*dim);
   elmats(0,1)->SetSize(dof_u*dim, dof_p);
   elmats(1,0)->SetSize(dof_p, dof_u*dim);
   elmats(1,1)->SetSize(dof_p, dof_p);

   *elmats(0,0) = 0.0;
   *elmats(0,1) = 0.0;
   *elmats(1,0) = 0.0;
   *elmats(1,1) = 0.0;

   shape.SetSize (dof_u);
   nor.SetSize (dim);
   fnor.SetSize (dim);
   
   DSh_u.SetSize(dof_u, dim);
   DS_u.SetSize(dof_u, dim);
   Sh_u.SetSize(dof_u);
   Sh_p.SetSize(dof_p);
   J0i.SetSize(dim);
   J.SetSize(dim);
   Jinv.SetSize(dim);
   JinvT.SetSize(dim);

   PMatI_u.UseExternalData(elfun[0]->GetData(), dof_u, dim);

   int intorder = 2*el[0]->GetOrder() + 3; 
   const IntegrationRule &ir = IntRules.Get(Tr.FaceGeom, intorder);

   double pres, dJ, norm; 

   for (int i = 0; i < ir.GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir.IntPoint(i);
      IntegrationPoint eip;
      Tr.Loc1.Transform(ip, eip);

      Tr.Face->SetIntPoint(&ip);
      CalcOrtho(Tr.Face->Jacobian(), nor);

      //Normalize vector
      norm = nor.Norml2();
      nor /= norm;

      Tr.Elem1->SetIntPoint(&eip);
      CalcInverse(Tr.Elem1->Jacobian(), J0i);

      el[0]->CalcDShape(eip, DSh_u);
      el[0]->CalcShape(eip, Sh_u);
      Mult(DSh_u, J0i, DS_u);
      MultAtB(PMatI_u, DS_u, J);

      CalcInverse(J, Jinv);
      CalcInverseTranspose(J, JinvT);

      dJ = J.Det();

      pres = function.Eval(*Tr.Face, ip);
            
      // u,u block
      for (int i_u = 0; i_u < dof_u; i_u++) {
      for (int i_dim = 0; i_dim < dim; i_dim++) {
         for (int j_u = 0; j_u < dof_u; j_u++) {
         for (int j_dim = 0; j_dim < dim; j_dim++) {

            for (int n=0; n<dim; n++) {
            for (int l=0; l<dim; l++) {
               (*elmats(0,0))(i_u + i_dim*dof_u, j_u + j_dim*dof_u) += dJ * pres * JinvT(i_dim,l) * JinvT(j_dim,n) * nor(n) * Sh_u(i_u) * DS_u(j_u,n) * ip.weight * Tr.Face->Weight();        
               (*elmats(0,0))(i_u + i_dim*dof_u, j_u + j_dim*dof_u) -= dJ * pres * JinvT(i_dim,n) * JinvT(j_dim,l) * nor(n) * Sh_u(i_u) * DS_u(j_u,n) * ip.weight * Tr.Face->Weight();
               
            }
            }
            
            
         }
         }
      }
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
            (*elvec[0])(dof_u*k+j) -= ip.weight * Tr.Face->Weight() * fnor(k) * shape(j);
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
