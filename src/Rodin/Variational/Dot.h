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

         std::unique_ptr<BasisOperator> getOperator(
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

   /**
    * @f$ u @f$ with @f$ n @f$ degrees of freedom
    * @f$ v @f$ with @f$ m @f$ degrees of freedom
    *
    * Stiffness matrix of @f$ m \times n @f$.
    */
   template <>
   class Dot<BasisOperator, BasisOperator>
   {
      public:
         Dot(std::unique_ptr<BasisOperator> lhs, std::unique_ptr<BasisOperator> rhs);

         Dot(const Dot& other);

         Dot(Dot&& other);

         operator mfem::DenseMatrix&()
         {
            return m_result;
         }

         operator const mfem::DenseMatrix&() const
         {
            return m_result;
         }

      private:
         mfem::DenseMatrix m_result;
   };
   Dot(std::unique_ptr<BasisOperator>, std::unique_ptr<BasisOperator>)
      -> Dot<BasisOperator, BasisOperator>;

   template <>
   class Dot<
      ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>
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

         mfem::DenseMatrix getElementMatrix(
               const mfem::FiniteElement& trialElement, const mfem::FiniteElement& testElement,
               mfem::ElementTransformation& trans) const
         {
            assert(m_trial->getRangeShape() == m_test->getRangeShape());
            return Dot<BasisOperator, BasisOperator>(
                  m_trial->getOperator(trialElement, trans),
                  m_test->getOperator(testElement, trans));
         }

         Dot* copy() const noexcept override
         {
            return new Dot(*this);
         }

      private:
         std::unique_ptr<ShapeFunctionBase<TrialSpace>> m_trial;
         std::unique_ptr<ShapeFunctionBase<TestSpace>>  m_test;
   };
   Dot(const ShapeFunctionBase<TrialSpace>&, const ShapeFunctionBase<TestSpace>&)
      -> Dot<ShapeFunctionBase<TrialSpace>, ShapeFunctionBase<TestSpace>>;
}

#endif
