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

   bool
   UnaryMinus<LinearFormIntegratorBase>::isSupported(Linear::Assembly::Type t)
   const
   {
      return m_op->isSupported(t);
   }

   void
   UnaryMinus<LinearFormIntegratorBase>::getElementVector(const Linear::Assembly::Device& as)
   const
   {
      m_op->getElementVector(as);
      as.vec *= -1.0;
   }

   void
   UnaryMinus<LinearFormIntegratorBase>::getElementVector(const Linear::Assembly::Common& as)
   const
   {
      m_op->getElementVector(as);
      as.vec *= -1.0;
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
   ::getElementMatrix(const Bilinear::Assembly::Common& as) const
   {
      m_op->getElementMatrix(as);
      as.mat.Neg();
   }

   UnaryMinus<BilinearFormIntegratorBase> operator-(const BilinearFormIntegratorBase& op)
   {
      return UnaryMinus(op);
   }

   UnaryMinus<LinearFormIntegratorSum>
   operator-(const LinearFormIntegratorSum& op)
   {
      return UnaryMinus<LinearFormIntegratorSum>(op);
   }

   UnaryMinus<BilinearFormIntegratorSum> operator-(const BilinearFormIntegratorSum& op)
   {
      return UnaryMinus<BilinearFormIntegratorSum>(op);
   }
}
