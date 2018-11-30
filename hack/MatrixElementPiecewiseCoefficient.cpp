#include "MatrixElementPiecewiseCoefficient.hpp"
#include <cassert>


mfem::DenseMatrix quat2rot(const mfem::Vector& q)
{
   MFEM_ASSERT(q.Size()==4, "quat2rot: Dimension of quaternion should be 4");
   //MFEM_ASSERT(vecisnorm(q), "quat2rot: quaternion is not normalized");
   mfem::DenseMatrix Q(3);

   double w=q(0);
   double x=q(1);
   double y=q(2);
   double z=q(3);

   double x2=x*x;
   double y2=y*y;
   double z2=z*z;
   double xy=x*y;
   double xz=x*z;
   double yz=y*z;
   double wx=w*x;
   double wy=w*y;
   double wz=w*z;

   Q(0,0)=1-2*y2-2*z2;
   Q(1,0)=2*xy-2*wz;
   Q(2,0)=2*xz+2*wy;

   Q(0,1)=2*xy+2*wz;
   Q(1,1)=1-2*x2-2*z2;
   Q(2,1)=2*yz-2*wx;

   Q(0,2)=2*xz-2*wy;
   Q(1,2)=2*yz+2*wx;
   Q(2,2)=1-2*x2-2*y2;

   return Q;
}

void MatrixElementPiecewiseCoefficient::Eval(
   mfem::DenseMatrix &K,
   mfem::ElementTransformation& T,
   const mfem::IntegrationPoint &ip)
{
   std::unordered_map<int,mfem::Vector>::iterator iter = heartConductivities_.find(T.Attribute);
   if (iter != heartConductivities_.end()) {
      mfem::Vector direction(3);
      if (1) {
         p_gf_->GetVectorValue(T.ElementNo, ip, direction);
      } else {
         direction = 0.0;
      }

      mfem::Vector quat(4);
      double w2 = 1;
      for (int ii=0; ii<3; ii++) {
         quat(ii+1) = direction(ii);
         w2 -= direction(ii)*direction(ii);
      }
      quat(0) = sqrt(w2);
         
      mfem::DenseMatrix VVV = quat2rot(quat);
      MultADAt(VVV,iter->second,K);
   }
   else {
      std::unordered_map<int,double>::iterator iter = bathConductivities_.find(T.Attribute);
      assert(iter != bathConductivities_.end());
      K=0.0;
      for (int ii=0; ii<3; ii++) {
         K(ii,ii) = iter->second;
      }
   }
}

