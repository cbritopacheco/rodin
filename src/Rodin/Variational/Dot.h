/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DOT_H
#define RODIN_VARIATIONAL_DOT_H

#include "Rodin/FormLanguage/Base.h"

#include "ForwardDecls.h"
#include "Function.h"
#include "ShapeFunction.h"
#include "ScalarFunction.h"

#include "Exceptions.h"

namespace Rodin::Variational
{
   template <>
   class Dot<FunctionBase, FunctionBase> : public ScalarFunctionBase
   {
      public:
         Dot(const FunctionBase& a, const FunctionBase& b);

         Dot(const Dot& other);

         Dot(Dot&& other);

         Dot& traceOf(const std::set<int>& attrs) override;

         virtual double getValue(
            mfem::ElementTransformation& trans,
            const mfem::IntegrationPoint& ip) const override;

         virtual Dot* copy() const noexcept override
         {
            return new Dot(*this);
         }
      private:
         std::unique_ptr<FunctionBase> m_a, m_b;
   };
   Dot(const FunctionBase&, const FunctionBase&) -> Dot<FunctionBase, FunctionBase>;

   /**
    * @brief Dot product between a FunctionBase and a ShapeFunctionBase.
    *
    * @f[
    *    \Lambda : A(u)
    * @f]
    * with @f$ A(u) \in \mathbb{R}^{p \times q} @f$,
    * @f$ \Lambda \in \mathbb{R}^{p \times q} @f$.
    */
   template <ShapeFunctionSpaceType Space>
   class Dot<FunctionBase, ShapeFunctionBase<Space>>
      : public ShapeFunctionBase<Space>
   {
      public:
         Dot(const FunctionBase& lhs, const ShapeFunctionBase<Space>& rhs)
            : m_f(lhs.copy()), m_u(rhs.copy())
         {
            if (lhs.getRangeShape() != rhs.getRangeShape())
               RangeShapeMismatchException(lhs.getRangeShape(), rhs.getRangeShape()).raise();
         }

         Dot(const ShapeFunctionBase<Space>& lhs, const FunctionBase& rhs)
            : Dot(rhs, lhs)
         {}

         Dot(const Dot& other)
            :  ShapeFunctionBase<Space>(other),
               m_f(other.m_f->copy()), m_u(other.m_u->copy())
         {}

         Dot(Dot&& other)
            :  ShapeFunctionBase<Space>(std::move(other)),
               m_f(std::move(other.m_lhs)), m_u(std::move(other.m_rhs))
         {}

         const FunctionBase& getFunction() const
         {
            return *m_f;
         }

         const ShapeFunctionBase<Space>& getShapeFunction() const
         {
            return *m_u;
         }

         const ShapeFunctionBase<Space>& getLeaf() const override
         {
            return getShapeFunction().getLeaf();
         }

         int getRows() const override
         {
            return 1;
         }

         int getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            return getShapeFunction().getDOFs(fe, trans);
         }

         int getColumns() const override
         {
            return 1;
         }

         virtual void getOperator(
               DenseBasisOperator& op,
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip,
               ShapeComputator& compute) const override
         {
            mfem::DenseMatrix v;
            getFunction().getValue(v, trans, ip);
            DenseBasisOperator tmp;
            getShapeFunction().getOperator(tmp, fe, trans, ip, compute);
            int opDofs = getDOFs(fe, trans);
            op.setSize(1, 1, opDofs);
            for (int i = 0; i < opDofs; i++)
               op(0, 0, i) = v * tmp(i);
         }

         FiniteElementSpaceBase& getFiniteElementSpace() override
         {
            return m_u->getFiniteElementSpace();
         }

         const FiniteElementSpaceBase& getFiniteElementSpace() const override
         {
            return getShapeFunction().getFiniteElementSpace();
         }

         virtual Dot* copy() const noexcept override
         {
            return new Dot(*this);
         }
      private:
         std::unique_ptr<FunctionBase> m_f;
         std::unique_ptr<ShapeFunctionBase<Space>> m_u;
   };
   template <ShapeFunctionSpaceType Space>
   Dot(const FunctionBase&, const ShapeFunctionBase<Space>&)
      -> Dot<FunctionBase, ShapeFunctionBase<Space>>;
   template <ShapeFunctionSpaceType Space>
   Dot(const ShapeFunctionBase<Space>&, const FunctionBase&)
      -> Dot<FunctionBase, ShapeFunctionBase<Space>>;

   template <>
   class Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>
      : public FormLanguage::Base
   {
      public:
         Dot(const ShapeFunctionBase<TrialSpace>& lhs, const ShapeFunctionBase<TestSpace>& rhs)
            : m_trial(lhs.copy()), m_test(rhs.copy())
         {
            if (lhs.getRangeShape() != rhs.getRangeShape())
               RangeShapeMismatchException(lhs.getRangeShape(), rhs.getRangeShape()).raise();
         }

         Dot(const Dot& other)
            :  Base(other),
               m_trial(other.m_trial->copy()), m_test(other.m_test->copy())
         {}

         Dot(Dot&& other)
            :  Base(std::move(other)),
               m_trial(std::move(other.m_trial)), m_test(std::move(other.m_test))
         {}

         ShapeFunctionBase<TrialSpace>& getLHS()
         {
            assert(m_trial);
            return *m_trial;
         }

         ShapeFunctionBase<TestSpace>& getRHS()
         {
            assert(m_test);
            return *m_test;
         }

         const ShapeFunctionBase<TrialSpace>& getLHS() const
         {
            assert(m_trial);
            return *m_trial;
         }

         const ShapeFunctionBase<TestSpace>& getRHS() const
         {
            assert(m_test);
            return *m_test;
         }

         virtual mfem::DenseMatrix getElementMatrix(
               const mfem::FiniteElement& trialElement,
               const mfem::FiniteElement& testElement,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip,
               ShapeComputator& compute) const
         {
            assert(m_trial->getRangeShape() == m_test->getRangeShape());
            DenseBasisOperator trialOp, testOp;
            m_trial->getOperator(trialOp, trialElement, trans, ip, compute);
            m_test->getOperator(testOp, testElement, trans, ip, compute);
            mfem::DenseMatrix result(testOp.getDOFs(), trialOp.getDOFs());
            for (int i = 0; i < testOp.getDOFs(); i++)
               for (int j = 0; j < trialOp.getDOFs(); j++)
                  result(i, j) = testOp(i) * trialOp(j);
            return result;
         }

         virtual Dot* copy() const noexcept override
         {
            return new Dot(*this);
         }

      private:
         std::unique_ptr<ShapeFunctionBase<TrialSpace>> m_trial;
         std::unique_ptr<ShapeFunctionBase<TestSpace>>  m_test;
   };
   Dot(const ShapeFunctionBase<TrialSpace>&, const ShapeFunctionBase<TestSpace>&)
      -> Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>;

   /* ||-- OPTIMIZATIONS -----------------------------------------------------
    * Dot<FunctionBase, ShapeFunctionBase<Space>>
    * ---------------------------------------------------------------------->>
    */

   /**
    * @internal
    *
    * @f[
    *    f \cdot u
    * @f]
    * where $f$ is a function (scalar or vector valued).
    */
   template <class FES, ShapeFunctionSpaceType Space>
   class Dot<FunctionBase, ShapeFunction<FES, Space>>
      : public Dot<FunctionBase, ShapeFunctionBase<Space>>
   {
      public:
         constexpr
         Dot(const FunctionBase& f, const ShapeFunction<FES, Space>& u)
            : Dot<FunctionBase, ShapeFunctionBase<Space>>(f, u)
         {
            if (f.getRangeType() != RangeType::Scalar && f.getRangeType() != RangeType::Vector)
               UnexpectedRangeTypeException({RangeType::Scalar, RangeType::Vector}, f.getRangeType()).raise();
         }

         constexpr
         Dot(const Dot& other)
            : Dot<FunctionBase, ShapeFunctionBase<Space>>(other)
         {}

         constexpr
         Dot(Dot&& other)
            : Dot<FunctionBase, ShapeFunctionBase<Space>>(other)
         {}

         Dot* copy() const noexcept override
         {
            return new Dot(*this);
         }
   };
   template <class FES, ShapeFunctionSpaceType Space>
   Dot(const FunctionBase&, const ShapeFunction<FES, Space>&)
      -> Dot<FunctionBase, ShapeFunction<FES, Space>>;

   /* <<-- OPTIMIZATIONS -----------------------------------------------------
    * Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>
    * ----------------------------------------------------------------------||
    */

   /* ||-- OPTIMIZATIONS -----------------------------------------------------
    * Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>
    * ---------------------------------------------------------------------->>
    */

   /**
    * @internal
    *
    * Represents the following expression:
    * @f[
    *    \nabla u \cdot \nabla v
    * @f]
    */
   template <class FES>
   class Dot<Grad<ShapeFunction<FES, TrialSpace>>, Grad<ShapeFunction<FES, TestSpace>>>
      : public Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>
   {
      public:
         constexpr
         Dot(const Grad<ShapeFunction<FES, TrialSpace>>& lhs, const Grad<ShapeFunction<FES, TestSpace>>& rhs)
            : Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>(lhs, rhs)
         {
            if (lhs.getRangeShape() != rhs.getRangeShape())
               RangeShapeMismatchException(lhs.getRangeShape(), rhs.getRangeShape()).raise();
         }

         constexpr
         Dot(const Dot& other)
            : Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>(other)
         {}

         constexpr
         Dot(Dot&& other)
            : Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>(std::move(other))
         {}
   };
   template <class FES>
   Dot(const Grad<ShapeFunction<FES, TrialSpace>>&, const Grad<ShapeFunction<FES, TestSpace>>&)
      -> Dot<Grad<ShapeFunction<FES, TrialSpace>>, Grad<ShapeFunction<FES, TestSpace>>>;

   /* <<-- OPTIMIZATIONS -----------------------------------------------------
    * Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>
    * ----------------------------------------------------------------------||
    */
}

#endif
