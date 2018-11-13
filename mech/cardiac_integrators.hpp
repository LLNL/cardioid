#ifndef CARDIAC_INTEG
#define CARDIAC_INTEG

#include "mfem.hpp"
#include "cardiac_coefficients.hpp"

namespace mfem
{

//  The transversely isotropic cardiac hyperelasticity model
class CardiacModel 
{
protected:
   double C_1, b_ff, b_ss, b_nn, b_fs, b_fn, b_ns;

   mutable DenseMatrix C, E, I, JT, PK2, FinvT, B, orth, orth_transpose, dq_dE, dq_dF;
   mutable DenseMatrix Ftilde, Ntilde, N;
   mutable Array4D<double> dP_dF, dFtilde_dF, dE_dFtilde, dNtilde_dF;
   
public:
   CardiacModel(double _C_1, double _b_ff, double _b_ss, double _b_nn,
                double _b_fs, double _b_fn, double _b_ns)
      : C_1(_C_1), b_ff(_b_ff), b_ss(_b_ss), b_nn(_b_nn),
        b_fs(_b_fs), b_fn(_b_fn), b_ns(_b_ns) { }

   virtual void EvalP(const DenseMatrix &J, const double pres, const Vector &fiber, DenseMatrix &P) const;

   virtual void AssembleH(const DenseMatrix &J, const double pres, const Vector &fiber, const DenseMatrix &DS, const Vector Sh_p, const double weight, const Array2D<DenseMatrix *>&elmats) const;

   virtual void GenerateTransform(const Vector &fiber, DenseMatrix &Q, DenseMatrix &QT) const;

   virtual ~CardiacModel();

};


/// Pressure boundary integrator 
class PressureBoundaryNLFIntegrator : public BlockNonlinearFormIntegrator
{
private:
   Coefficient &function;
   VectorCoefficient &vol_function;
   mutable DenseMatrix DSh_u, DS_u, J0i, J, Jinv, JinvT, PMatI_u;
   mutable Vector shape, nor, fnor, Sh_p, Sh_u;
   
public:
   PressureBoundaryNLFIntegrator(Coefficient &f, VectorCoefficient &vf) : function(f), vol_function(vf) { }

   virtual void AssembleElementVector(const Array<const FiniteElement *> &el,
                                      ElementTransformation &Tr,
                                      const Array<const Vector *> &elfun, 
                                      const Array<Vector *> &elvec);

   virtual void AssembleFaceVector(const Array<const FiniteElement *> &el1,
                                   const Array<const FiniteElement *> &el2,
                                   FaceElementTransformations &Tr,
                                   const Array<const Vector *> &elfun, 
                                   const Array<Vector *> &elvec);

   virtual void AssembleFaceGrad(const Array<const FiniteElement*> &el1,
                                 const Array<const FiniteElement *> &el2,
                                 FaceElementTransformations &Tr,
                                 const Array<const Vector *> &elfun, 
                                 const Array2D<DenseMatrix *> &elmats);

   virtual double GetElementVolume(const Array<const FiniteElement *>&el,
                                   FaceElementTransformations &Tr,
                                   const Array<const Vector *>&elfun);
   
   virtual void AssembleVolumeGradient(const Array<const FiniteElement *> &el1,
                                       FaceElementTransformations &Tr,
                                       const Array<const Vector *> &elfun,
                                       const Array<Vector *> &elvect);
   
   
   virtual ~PressureBoundaryNLFIntegrator();
};


/// Active tension integrator 
class ActiveTensionNLFIntegrator : public BlockNonlinearFormIntegrator
{
private:
   QuadratureFunctionCoefficient &function;
   VectorCoefficient &Q;
   mutable DenseMatrix DSh_u, DS_u, J0i, J, P, PMatI_u, PMatO_u, PMatO_p, Qshap, AStress, AStress_current;
   mutable Vector Sh_p, shape;
   
public:
   ActiveTensionNLFIntegrator(QuadratureFunctionCoefficient &fq, VectorCoefficient &fv) : function(fq), Q(fv) { }

   virtual void AssembleElementVector(const Array<const FiniteElement *> &el,
                                      ElementTransformation &Tr,
                                      const Array<const Vector *> &elfun, 
                                      const Array<Vector *> &elvec);

   virtual void AssembleElementGrad(const Array<const FiniteElement*> &el,
                                    ElementTransformation &Tr,
                                    const Array<const Vector *> &elfun, 
                                    const Array2D<DenseMatrix *> &elmats);
   
   virtual ~ActiveTensionNLFIntegrator();
};

/// Incompressible nonlinear cardiac integrator
class CardiacNLFIntegrator : public BlockNonlinearFormIntegrator
{
private:
   CardiacModel *model;
   mutable DenseMatrix DSh_u, DS_u, J0i, J, P, PMatI_u, PMatO_u, PMatO_p, F, Finv, FinvT, B;
   mutable Vector Sh_p, Qvec;
   VectorCoefficient *Q;

public:
   CardiacNLFIntegrator(CardiacModel *m, VectorCoefficient &fib) : model(m), Q(&fib) { }

   // Assemble the element residual
   virtual void AssembleElementVector(const Array<const FiniteElement *> &el,
                                      ElementTransformation &Tr,
                                      const Array<const Vector *> &elfun, 
                                      const Array<Vector *> &elvec);

   // Assemble the element gradient
   virtual void AssembleElementGrad(const Array<const FiniteElement*> &el,
                                    ElementTransformation &Tr,
                                    const Array<const Vector *> &elfun, 
                                    const Array2D<DenseMatrix *> &elmats);

   virtual ~CardiacNLFIntegrator();
};


}

#endif
