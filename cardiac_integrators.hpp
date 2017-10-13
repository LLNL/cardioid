#ifndef CARDIAC_INTEG
#define CARDIAC_INTEG

#include "mfem.hpp"

namespace mfem
{

//  The transversely isotropic cardiac hyperelasticity model
class CardiacModel 
{
protected:
   double C_1, b_ff, b_ss, b_nn, b_fs, b_fn, b_ns;

   mutable DenseMatrix C, E, I, JT, PK2, FinvT;

public:
   CardiacModel(double _C_1, double _b_ff, double _b_ss, double _b_nn,
                double _b_fs, double _b_fn, double _b_ns)
      : C_1(_C_1), b_ff(_b_ff), b_ss(_b_ss), b_nn(_b_nn),
        b_fs(_b_fs), b_fn(_b_fn), b_ns(_b_ns) { }

   virtual void EvalP(const DenseMatrix &J, const double pres, const Vector &fiber, DenseMatrix &P) const;

   virtual double EvalW(const DenseMatrix &J) const;

   virtual void EvalP(const DenseMatrix &J, DenseMatrix &P) const;

   virtual void AssembleH(const DenseMatrix &J, const DenseMatrix &DS,
                          const double weight, DenseMatrix &A) const;


};


/// Body force integrator 
class BodyForceNLFIntegrator : public BlockNonlinearFormIntegrator
{
private:
   VectorCoefficient &function;
   Vector shape, Qvec;
   
public:
   BodyForceNLFIntegrator(VectorCoefficient &f) : function(f) { }


   virtual void AssembleElementVector(Array<const FiniteElement *> &el,
                                      ElementTransformation &Tr,
                                      Array<Vector *> &elfun, 
                                      Array<Vector *> &elvec);

   virtual ~BodyForceNLFIntegrator();
};

/// Pressure boundary integrator 
class PressureBoundaryNLFIntegrator : public BlockNonlinearFormIntegrator
{
private:
   Coefficient &function;
   DenseMatrix DSh_u, DS_u, J0i, J, Jinv, PMatI_u;
   Vector shape, nor, fnor;
   
public:
   PressureBoundaryNLFIntegrator(Coefficient &f) : function(f) { }

   virtual void AssembleElementVector(Array<const FiniteElement *> &el,
                                      ElementTransformation &Tr,
                                      Array<Vector *> &elfun, 
                                      Array<Vector *> &elvec);

   virtual void AssembleRHSElementVector(Array<const FiniteElement *> &el,
                                         FaceElementTransformations &Tr,
                                         Array<Vector *> &elfun, 
                                         Array<Vector *> &elvec);

   virtual ~PressureBoundaryNLFIntegrator();
};

/// Traction boundary integrator 
class TractionBoundaryNLFIntegrator : public BlockNonlinearFormIntegrator
{
private:
   VectorCoefficient &function;
   DenseMatrix DSh_u, DS_u, J0i, J, Jinv, PMatI_u;
   Vector shape, nor, fnor;
   
public:
   TractionBoundaryNLFIntegrator(VectorCoefficient &f) : function(f) { }

   virtual void AssembleElementVector(Array<const FiniteElement *> &el,
                                      ElementTransformation &Tr,
                                      Array<Vector *> &elfun, 
                                      Array<Vector *> &elvec);

   virtual void AssembleRHSElementVector(Array<const FiniteElement *> &el,
                                         FaceElementTransformations &Tr,
                                         Array<Vector *> &elfun, 
                                         Array<Vector *> &elvec);

   virtual ~TractionBoundaryNLFIntegrator();
};

/// Active tension integrator 
class ActiveTensionNLFIntegrator : public BlockNonlinearFormIntegrator
{
private:
   MatrixCoefficient &function;
   DenseMatrix DSh_u, DS_u, J0i, J, P, PMatI_u, PMatO_u, PMatO_p, Qshap, AStress, AStress_current;
   Vector Sh_p, shape;
   
public:
   ActiveTensionNLFIntegrator(MatrixCoefficient &f) : function(f) { }

   virtual void AssembleElementVector(Array<const FiniteElement *> &el,
                                      ElementTransformation &Tr,
                                      Array<Vector *> &elfun, 
                                      Array<Vector *> &elvec);

   virtual ~ActiveTensionNLFIntegrator();
};

/// Incompressible nonlinear cardiac integrator
class CardiacNLFIntegrator : public BlockNonlinearFormIntegrator
{
private:
   CardiacModel *model;
   DenseMatrix DSh_u, DS_u, J0i, J, P, PMatI_u, PMatO_u, PMatO_p;
   Vector Sh_p, Qvec;
   VectorCoefficient *Q;

public:
   CardiacNLFIntegrator(CardiacModel *m, VectorCoefficient &fib) : model(m), Q(&fib) { }

   virtual void AssembleElementVector(Array<const FiniteElement *> &el,
                                      ElementTransformation &Tr,
                                      Array<Vector *> &elfun, 
                                      Array<Vector *> &elvec);

   virtual ~CardiacNLFIntegrator();
};


}

#endif
