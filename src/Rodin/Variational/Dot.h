/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DOT_H
#define RODIN_VARIATIONAL_DOT_H

#include "FormLanguage/Base.h"

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

         double getValue(
            mfem::ElementTransformation& trans,
            const mfem::IntegrationPoint& ip) const override;

         Dot* copy() const noexcept override
         {
            return new Dot(*this);
         }
      private:
         std::unique_ptr<FunctionBase> m_a, m_b;
   };
   Dot(const FunctionBase&, const FunctionBase&) -> Dot<FunctionBase, FunctionBase>;

   /**
    * @brief Dot product between FunctionBase and ShapeFunctionBase.
    */
   template <ShapeFunctionSpaceType Space>
   class Dot<FunctionBase, ShapeFunctionBase<Space>>
      : public ShapeFunctionBase<Space>
   {
      public:
         Dot(const ShapeFunctionBase<Space>& lhs, const FunctionBase& rhs)
            : Dot(rhs, lhs)
         {
            if (lhs.getRangeShape() != rhs.getRangeShape())
               RangeShapeMismatchException(lhs.getRangeShape(), rhs.getRangeShape()).raise();
         }

         Dot(const FunctionBase& lhs, const ShapeFunctionBase<Space>& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Dot(const Dot& other)
            :  ShapeFunctionBase<Space>(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Dot(Dot&& other)
            :  ShapeFunctionBase<Space>(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         const FunctionBase& getLHS() const
         {
            return *m_lhs;
         }

         const ShapeFunctionBase<Space>& getRHS() const
         {
            return *m_rhs;
         }

         const ShapeFunctionBase<Space>& getLeaf() const override
         {
            return getRHS().getLeaf();
         }

         int getRows() const override
         {
            return 1;
         }

         int getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            return getRHS().getDOFs(fe, trans);
         }

         int getColumns() const override
         {
            return 1;
         }

         std::unique_ptr<Internal::Rank3Operator> getOperator(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans) const override
         {
            mfem::DenseMatrix v;
            getLHS().getValue(v, trans, trans.GetIntPoint());
            return getRHS().getOperator(fe, trans)->MatrixDot(v);
         }

         FiniteElementSpaceBase& getFiniteElementSpace() override
         {
            return m_rhs->getFiniteElementSpace();
         }

         const FiniteElementSpaceBase& getFiniteElementSpace() const override
         {
            return getRHS().getFiniteElementSpace();
         }

         Dot* copy() const noexcept override
         {
            return new Dot(*this);
         }
      private:
         std::unique_ptr<FunctionBase> m_lhs;
         std::unique_ptr<ShapeFunctionBase<Space>> m_rhs;
   };
   template <ShapeFunctionSpaceType Space>
   Dot(const FunctionBase&, const ShapeFunctionBase<Space>&)
      -> Dot<FunctionBase, ShapeFunctionBase<Space>>;
   template <ShapeFunctionSpaceType Space>
   Dot(const ShapeFunctionBase<Space>&, const FunctionBase&)
      -> Dot<FunctionBase, ShapeFunctionBase<Space>>;

   template <>
   class Dot<
      ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>
      : public FormLanguage::Base
   {
      public:
         Dot(const ShapeFunctionBase<TrialSpace>& lhs, const ShapeFunctionBase<TestSpace>& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {
            if (lhs.getRangeShape() != rhs.getRangeShape())
               RangeShapeMismatchException(lhs.getRangeShape(), rhs.getRangeShape()).raise();
         }

         Dot(const Dot& other)
            :  Base(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Dot(Dot&& other)
            :  Base(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         ShapeFunctionBase<TrialSpace>& getLHS()
         {
            assert(m_lhs);
            return *m_lhs;
         }

         ShapeFunctionBase<TestSpace>& getRHS()
         {
            assert(m_rhs);
            return *m_rhs;
         }

         const ShapeFunctionBase<TrialSpace>& getLHS() const
         {
            assert(m_lhs);
            return *m_lhs;
         }

         const ShapeFunctionBase<TestSpace>& getRHS() const
         {
            assert(m_rhs);
            return *m_rhs;
         }

         mfem::DenseMatrix getElementMatrix(
               const mfem::FiniteElement& trialElement, const mfem::FiniteElement& testElement,
               mfem::ElementTransformation& trans) const
         {
            auto& trial = getLHS();
            auto& test  = getRHS();
            assert(trial.getRangeShape() == test.getRangeShape());
            return test.getOperator(testElement, trans)->OperatorDot(
                     *trial.getOperator(trialElement, trans));
         }

         Dot* copy() const noexcept override
         {
            return new Dot(*this);
         }

      private:
         std::unique_ptr<ShapeFunctionBase<TrialSpace>> m_lhs;
         std::unique_ptr<ShapeFunctionBase<TestSpace>>  m_rhs;
   };
   Dot(const ShapeFunctionBase<TrialSpace>&, const ShapeFunctionBase<TestSpace>&)
      -> Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>;
}

#endif
