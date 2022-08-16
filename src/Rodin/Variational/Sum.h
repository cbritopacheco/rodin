/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_SUM_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_SUM_H

#include "Rodin/FormLanguage/Base.h"

#include "ForwardDecls.h"
#include "Function.h"
#include "ShapeFunction.h"

#include "Exceptions.h"

namespace Rodin::Variational
{
   template <>
   class Sum<FunctionBase, FunctionBase> : public FunctionBase
   {
      public:
         Sum(const FunctionBase& lhs, const FunctionBase& rhs);

         Sum(const Sum& other);

         Sum(Sum&& other);

         Sum& operator+=(const FunctionBase& lhs);

         RangeShape getRangeShape() const override;

         Sum& traceOf(const std::set<int>& attrs) override;

         void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) const override;

         Sum* copy() const noexcept override
         {
            return new Sum(*this);
         }

      private:
         std::unique_ptr<FunctionBase> m_lhs;
         std::unique_ptr<FunctionBase> m_rhs;
   };
   Sum(const FunctionBase&, const FunctionBase&) -> Sum<FunctionBase, FunctionBase>;

   Sum<FunctionBase, FunctionBase>
   operator+(const FunctionBase& lhs, const FunctionBase& rhs);

   template <class T>
   std::enable_if_t<std::is_arithmetic_v<T>, Sum<FunctionBase, FunctionBase>>
   operator+(const FunctionBase& lhs, T v)
   {
      return Sum(lhs, ScalarFunction(v));
   }

   template <class T>
   std::enable_if_t<std::is_arithmetic_v<T>, Sum<FunctionBase, FunctionBase>>
   operator+(T v, const FunctionBase& rhs)
   {
      return Sum(ScalarFunction(v), rhs);
   }

   template <ShapeFunctionSpaceType Space>
   class Sum<ShapeFunctionBase<Space>, ShapeFunctionBase<Space>>
      : public ShapeFunctionBase<Space>
   {
      public:
         Sum(const ShapeFunctionBase<Space>& lhs, const ShapeFunctionBase<Space>& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {
            assert(lhs.getLeaf().getUUID() == rhs.getLeaf().getUUID());
            if (lhs.getRangeShape() != rhs.getRangeShape())
               RangeShapeMismatchException(lhs.getRangeShape(), rhs.getRangeShape()).raise();
         }

         Sum(const Sum& other)
            : m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Sum(Sum&& other)
            : m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         ShapeFunctionBase<Space>& getLHS()
         {
            return *m_lhs;
         }

         ShapeFunctionBase<Space>& getRHS()
         {
            return *m_rhs;
         }

         const ShapeFunctionBase<Space>& getLHS() const
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
            return getLHS().getRows();
         }

         int getColumns() const override
         {
            return getLHS().getColumns();
         }

         int getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            assert(getLHS().getDOFs(fe, trans) == getRHS().getDOFs(fe, trans));
            return getLHS().getDOFs(fe, trans);
         }

         void getOperator(
               DenseBasisOperator& op,
               const mfem::FiniteElement& fe,
               ShapeComputator& comp) const override
         {
            getLHS().getOperator(op, fe, comp);
            DenseBasisOperator tmp;
            getRHS().getOperator(tmp, fe, comp);
            op += tmp;
         }

         FiniteElementSpaceBase& getFiniteElementSpace() override
         {
            return getLHS().getFiniteElementSpace();
         }

         const FiniteElementSpaceBase& getFiniteElementSpace() const override
         {
            return getLHS().getFiniteElementSpace();
         }

         Sum* copy() const noexcept override
         {
            return new Sum(*this);
         }
      private:
         std::unique_ptr<ShapeFunctionBase<Space>> m_lhs;
         std::unique_ptr<ShapeFunctionBase<Space>> m_rhs;
   };
   template <ShapeFunctionSpaceType Space>
   Sum(const ShapeFunctionBase<Space>&, const ShapeFunctionBase<Space>&)
      -> Sum<ShapeFunctionBase<Space>, ShapeFunctionBase<Space>>;

   template <ShapeFunctionSpaceType Space>
   Sum<ShapeFunctionBase<Space>, ShapeFunctionBase<Space>>
   operator+(const ShapeFunctionBase<Space>& lhs, const ShapeFunctionBase<Space>& rhs)
   {
      return Sum(lhs, rhs);
   }
}

#endif
