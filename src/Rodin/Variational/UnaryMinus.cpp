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
   UnaryMinus<ScalarFunctionBase>
   ::getValue(mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
   const
   {
      return -1.0 * getOperand().getValue(trans, ip);
   }

   UnaryMinus<ScalarFunctionBase>
   operator-(const ScalarFunctionBase& op)
   {
      return UnaryMinus(op);
   }

   void
   UnaryMinus<VectorFunctionBase>
   ::getValue(
         mfem::Vector& value,
         mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip)
   const
   {
      getOperand().getValue(value, trans, ip);
      value.Neg();
   }

   UnaryMinus<VectorFunctionBase>
   operator-(const VectorFunctionBase& op)
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

   Sum<ScalarFunctionBase, ScalarFunctionBase>
   operator-(const ScalarFunctionBase& lhs, const ScalarFunctionBase& rhs)
   {
      return Sum<ScalarFunctionBase, ScalarFunctionBase>(lhs, UnaryMinus(rhs));
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
