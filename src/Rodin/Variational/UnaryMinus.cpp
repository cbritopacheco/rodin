/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "UnaryMinus.h"
#include "RangeShape.h"

namespace Rodin::Variational
{
   // ---- FunctionBase ------------------------------------------------------
   UnaryMinus<FunctionBase>::UnaryMinus(const FunctionBase& op)
      : m_op(op.copy())
   {}

   UnaryMinus<FunctionBase>::UnaryMinus(const UnaryMinus& other)
      :  FunctionBase(other),
         m_op(other.m_op->copy())
   {}

   UnaryMinus<FunctionBase>::UnaryMinus(UnaryMinus&& other)
      : FunctionBase(std::move(other)),
        m_op(std::move(other.m_op))
   {}

   void UnaryMinus<FunctionBase>::getValue(
         mfem::DenseMatrix& value,
         mfem::ElementTransformation& trans,
         const mfem::IntegrationPoint& ip) const
   {
      m_op->getValue(value, trans, ip);
      value.Neg();
   }

   RangeShape UnaryMinus<FunctionBase>::getRangeShape() const
   {
      return m_op->getRangeShape();
   }

   UnaryMinus<FunctionBase> operator-(const FunctionBase& op)
   {
      return UnaryMinus(op);
   }

   // ---- LinearFormIntegratorBase ------------------------------------------
   UnaryMinus<LinearFormIntegratorBase>::UnaryMinus(const LinearFormIntegratorBase& op)
      : m_op(op.copy())
   {}

   UnaryMinus<LinearFormIntegratorBase>::UnaryMinus(const UnaryMinus& other)
      :  LinearFormIntegratorBase(other),
         m_op(other.m_op->copy())
   {}

   UnaryMinus<LinearFormIntegratorBase>::UnaryMinus(UnaryMinus&& other)
      : LinearFormIntegratorBase(std::move(other)),
        m_op(std::move(other.m_op))
   {}

   const std::set<int>&
   UnaryMinus<LinearFormIntegratorBase>::getAttributes() const
   {
      return m_op->getAttributes();
   }

   IntegratorRegion
   UnaryMinus<LinearFormIntegratorBase>::getIntegratorRegion() const
   {
      return m_op->getIntegratorRegion();
   }

   const ShapeFunctionBase<TestSpace>&
   UnaryMinus<LinearFormIntegratorBase>::getTestFunction() const
   {
      return m_op->getTestFunction();
   }

   void
   UnaryMinus<LinearFormIntegratorBase>
   ::getElementVector(const mfem::FiniteElement& fe,
         mfem::ElementTransformation& trans, mfem::Vector& vec)
   const
   {
      m_op->getElementVector(fe, trans, vec);
      vec *= -1.0;
   }

   UnaryMinus<LinearFormIntegratorBase>
   operator-(const LinearFormIntegratorBase& op)
   {
      return UnaryMinus(op);
   }

   // ---- BilinearFormIntegratorBase ----------------------------------------
   UnaryMinus<BilinearFormIntegratorBase>::UnaryMinus(const BilinearFormIntegratorBase& op)
      : m_op(op.copy())
   {}

   UnaryMinus<BilinearFormIntegratorBase>::UnaryMinus(const UnaryMinus& other)
      :  BilinearFormIntegratorBase(other),
         m_op(other.m_op->copy())
   {}

   UnaryMinus<BilinearFormIntegratorBase>::UnaryMinus(UnaryMinus&& other)
      : BilinearFormIntegratorBase(std::move(other)),
        m_op(std::move(other.m_op))
   {}

   const ShapeFunctionBase<TrialSpace>&
   UnaryMinus<BilinearFormIntegratorBase>::getTrialFunction() const
   {
      return m_op->getTrialFunction();
   }

   const ShapeFunctionBase<TestSpace>&
   UnaryMinus<BilinearFormIntegratorBase>::getTestFunction() const
   {
      return m_op->getTestFunction();
   }

   const std::set<int>&
   UnaryMinus<BilinearFormIntegratorBase>::getAttributes() const
   {
      return m_op->getAttributes();
   }

   IntegratorRegion
   UnaryMinus<BilinearFormIntegratorBase>::getIntegratorRegion() const
   {
      return m_op->getIntegratorRegion();
   }

   void
   UnaryMinus<BilinearFormIntegratorBase>
   ::getElementMatrix(
         const mfem::FiniteElement& trial, const mfem::FiniteElement& test,
         mfem::ElementTransformation& trans, mfem::DenseMatrix& mat) const
   {
      m_op->getElementMatrix(trial, test, trans, mat);
      mat.Neg();
   }

   UnaryMinus<BilinearFormIntegratorBase> operator-(const BilinearFormIntegratorBase& op)
   {
      return UnaryMinus(op);
   }

   UnaryMinus<FormLanguage::LinearFormIntegratorSum>
   operator-(const FormLanguage::LinearFormIntegratorSum& op)
   {
      return UnaryMinus<FormLanguage::LinearFormIntegratorSum>(op);
   }

   UnaryMinus<FormLanguage::BilinearFormIntegratorSum> operator-(
         const FormLanguage::BilinearFormIntegratorSum& op)
   {
      return UnaryMinus<FormLanguage::BilinearFormIntegratorSum>(op);
   }
}
