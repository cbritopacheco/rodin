#include "Sum.h"

#include "UnaryMinus.h"

namespace Rodin::Variational::FormLanguage
{
   double
   UnaryMinus<ScalarCoefficientBase>
   ::getValue(mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
   const
   {
      return -1.0 * getOperand().getValue(trans, ip);
   }

   void
   UnaryMinus<LinearFormIntegratorBase>
   ::getElementVector(const mfem::FiniteElement& fe,
         mfem::ElementTransformation& trans, mfem::Vector& vec)
   {
      getOperand().getElementVector(fe, trans, vec);
      vec *= -1.0;
   }

   void
   UnaryMinus<BilinearFormIntegratorBase>
   ::getElementMatrix(
         const mfem::FiniteElement& trial, const mfem::FiniteElement& test,
         mfem::ElementTransformation& trans, mfem::DenseMatrix& mat)
   {
      getOperand().getElementMatrix(trial, test, trans, mat);
      mat *= -1.0;
   }

   Sum<ScalarCoefficientBase, ScalarCoefficientBase>
   operator-(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
   {
      return Sum<ScalarCoefficientBase, ScalarCoefficientBase>(lhs, UnaryMinus(rhs));
   }

   UnaryMinus<ScalarCoefficientBase>
   operator-(const ScalarCoefficientBase& op)
   {
      return UnaryMinus(op);
   }

   UnaryMinus<LinearFormIntegratorBase>
   operator-(const LinearFormIntegratorBase& op)
   {
      return UnaryMinus(op);
   }

   UnaryMinus<LinearFormIntegratorSum>
   operator-(const LinearFormIntegratorSum& op)
   {
      return UnaryMinus(op);
   }

   UnaryMinus<BilinearFormIntegratorBase> operator-(const BilinearFormIntegratorBase& op)
   {
      return UnaryMinus(op);
   }

   UnaryMinus<BilinearFormIntegratorSum> operator-(const BilinearFormIntegratorSum& op)
   {
      return UnaryMinus(op);
   }
}
