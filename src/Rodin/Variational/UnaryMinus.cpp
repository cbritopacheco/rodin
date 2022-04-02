/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Sum.h"

#include "UnaryMinus.h"

namespace Rodin::Variational
{
   double
   UnaryMinus<ScalarCoefficientBase>
   ::getValue(mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
   const
   {
      return -1.0 * getOperand().getValue(trans, ip);
   }

   UnaryMinus<ScalarCoefficientBase>
   operator-(const ScalarCoefficientBase& op)
   {
      return UnaryMinus(op);
   }

   void
   UnaryMinus<VectorCoefficientBase>
   ::getValue(
         mfem::Vector& value,
         mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
   const
   {
      getOperand().getValue(value, trans, ip);
      value.Neg();
   }

   UnaryMinus<VectorCoefficientBase>
   operator-(const VectorCoefficientBase& op)
   {
      return UnaryMinus(op);
   }

   void
   UnaryMinus<LinearFormIntegratorBase>
   ::getElementVector(const mfem::FiniteElement& fe,
         mfem::ElementTransformation& trans, mfem::Vector& vec)
   const
   {
      getOperand().getElementVector(fe, trans, vec);
      vec *= -1.0;
   }

   UnaryMinus<LinearFormIntegratorBase>
   operator-(const LinearFormIntegratorBase& op)
   {
      return UnaryMinus(op);
   }

   void
   UnaryMinus<BilinearFormIntegratorBase>
   ::getElementMatrix(
         const mfem::FiniteElement& trial, const mfem::FiniteElement& test,
         mfem::ElementTransformation& trans, mfem::DenseMatrix& mat) const
   {
      getOperand().getElementMatrix(trial, test, trans, mat);
      mat *= -1.0;
   }

   UnaryMinus<BilinearFormIntegratorBase> operator-(const BilinearFormIntegratorBase& op)
   {
      return UnaryMinus(op);
   }

   Sum<ScalarCoefficientBase, ScalarCoefficientBase>
   operator-(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
   {
      return Sum<ScalarCoefficientBase, ScalarCoefficientBase>(lhs, UnaryMinus(rhs));
   }

   UnaryMinus<FormLanguage::LinearFormIntegratorSum>
   operator-(const FormLanguage::LinearFormIntegratorSum& op)
   {
      return UnaryMinus(op);
   }

   UnaryMinus<FormLanguage::BilinearFormIntegratorSum> operator-(
         const FormLanguage::BilinearFormIntegratorSum& op)
   {
      return UnaryMinus(op);
   }
}
