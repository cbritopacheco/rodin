/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_SUM_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_SUM_H

#include <memory>
#include <type_traits>

#include "FormLanguage/Base.h"

#include "ForwardDecls.h"
#include "ShapeFunction.h"
#include "ScalarFunction.h"
#include "VectorFunction.h"
#include "MatrixFunction.h"

namespace Rodin::Variational
{
   template <>
   class Sum<ScalarFunctionBase, ScalarFunctionBase>
      : public ScalarFunctionBase
   {
      public:
         Sum(const ScalarFunctionBase& lhs, const ScalarFunctionBase& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Sum(const Sum& other)
            :  ScalarFunctionBase(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Sum(Sum&& other)
            :  ScalarFunctionBase(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         ScalarFunctionBase& getLHS()
         {
            return *m_lhs;
         }

         ScalarFunctionBase& getRHS()
         {
            return *m_rhs;
         }

         const ScalarFunctionBase& getLHS() const
         {
            return *m_lhs;
         }

         const ScalarFunctionBase& getRHS() const
         {
            return *m_rhs;
         }

         double getValue(
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override;

         Sum* copy() const noexcept override
         {
            return new Sum(*this);
         }

         Sum& operator+=(const ScalarFunctionBase& lhs)
         {
            auto sum = new Sum(lhs, *m_rhs);
            m_rhs.reset(sum);
            return *this;
         }

      private:
         std::unique_ptr<ScalarFunctionBase> m_lhs;
         std::unique_ptr<ScalarFunctionBase> m_rhs;
   };
   Sum(const ScalarFunctionBase&, const ScalarFunctionBase&)
      -> Sum<ScalarFunctionBase, ScalarFunctionBase>;

   Sum<ScalarFunctionBase, ScalarFunctionBase>
   operator+(const ScalarFunctionBase& lhs, const ScalarFunctionBase& rhs);

   template <>
   class Sum<VectorFunctionBase, VectorFunctionBase>
      : public VectorFunctionBase
   {
      public:
         Sum(const VectorFunctionBase& lhs, const VectorFunctionBase& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {
            assert(lhs.getDimension() == rhs.getDimension());
         }

         Sum(const Sum& other)
            :  VectorFunctionBase(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Sum(Sum&& other)
            :  VectorFunctionBase(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         VectorFunctionBase& getLHS()
         {
            return *m_lhs;
         }

         VectorFunctionBase& getRHS()
         {
            return *m_rhs;
         }

         const VectorFunctionBase& getLHS() const
         {
            return *m_lhs;
         }

         const VectorFunctionBase& getRHS() const
         {
            return *m_rhs;
         }

         int getDimension() const override
         {
            return m_lhs->getDimension();
         }

         void getValue(mfem::Vector& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override
         {
            getLHS().getValue(value, trans, ip);
            mfem::Vector tmp;
            getRHS().getValue(tmp, trans, ip);
            value += tmp;
         }

         Sum& operator+=(const VectorFunctionBase& lhs)
         {
            auto sum = new Sum(lhs, *m_rhs);
            m_rhs.reset(sum);
            return *this;
         }

         Sum* copy() const noexcept override
         {
            return new Sum(*this);
         }

      private:
         std::unique_ptr<VectorFunctionBase> m_lhs;
         std::unique_ptr<VectorFunctionBase> m_rhs;
   };
   Sum(const VectorFunctionBase&, const VectorFunctionBase&)
      -> Sum<VectorFunctionBase, VectorFunctionBase>;

   Sum<VectorFunctionBase, VectorFunctionBase>
   operator+(const VectorFunctionBase& lhs, const VectorFunctionBase& rhs);

   /**
    * @brief %Sum of two MatrixFunctionBase instances.
    */
   template <>
   class Sum<MatrixFunctionBase, MatrixFunctionBase>
      : public MatrixFunctionBase
   {
      public:
         Sum(const MatrixFunctionBase& lhs, const MatrixFunctionBase& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Sum(const Sum& other)
            :  MatrixFunctionBase(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Sum(Sum&& other)
            :  MatrixFunctionBase(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         MatrixFunctionBase& getLHS()
         {
            return *m_lhs;
         }

         MatrixFunctionBase& getRHS()
         {
            return *m_rhs;
         }

         const MatrixFunctionBase& getLHS() const
         {
            return *m_lhs;
         }

         const MatrixFunctionBase& getRHS() const
         {
            return *m_rhs;
         }

         int getRows() const override;

         int getColumns() const override;

         void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override;

         Sum& operator+=(const MatrixFunctionBase& lhs)
         {
            auto sum = new Sum(lhs, *m_rhs);
            m_rhs.reset(sum);
            return *this;
         }

         Sum* copy() const noexcept override
         {
            return new Sum(*this);
         }

      private:
         std::unique_ptr<MatrixFunctionBase> m_lhs;
         std::unique_ptr<MatrixFunctionBase> m_rhs;
   };
   Sum(const MatrixFunctionBase&, const MatrixFunctionBase&)
      -> Sum<MatrixFunctionBase, MatrixFunctionBase>;
   Sum<MatrixFunctionBase, MatrixFunctionBase>
   operator+(const MatrixFunctionBase& lhs, const MatrixFunctionBase& rhs);

   template <ShapeFunctionSpaceType Space>
   class Sum<ShapeFunctionBase<Space>, ShapeFunctionBase<Space>>
      : public ShapeFunctionBase<Space>
   {
      public:
         Sum(const ShapeFunctionBase<Space>& lhs, const ShapeFunctionBase<Space>& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {
            assert(lhs.getLeaf().getUUID() == rhs.getLeaf().getUUID());
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

         ShapeFunctionBase<Space>& getLeaf() override
         {
            return getRHS().getLeaf();
         }

         const ShapeFunctionBase<Space>& getLeaf() const override
         {
            return getRHS().getLeaf();
         }

         int getRows(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            assert(getLHS().getRows(fe, trans) == getRHS().getRows(fe, trans));
            return getLHS().getRows(fe, trans);
         }

         int getColumns(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            assert(getLHS().getColumns(fe, trans) == getRHS().getColumns(fe, trans));
            return getLHS().getColumns(fe, trans);
         }

         int getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            assert(getLHS().getDOFs(fe, trans) == getRHS().getDOFs(fe, trans));
            return getLHS().getDOFs(fe, trans);
         }

         std::unique_ptr<Rank3Operator> getOperator(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans) const override
         {
            return getLHS().getOperator(fe, trans)->OperatorSum(*getRHS().getOperator(fe, trans));
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

   template <class T>
   std::enable_if_t<std::is_arithmetic_v<T>, Sum<ScalarFunctionBase, ScalarFunctionBase>>
   operator+(const ScalarFunctionBase& lhs, T v)
   {
      return Sum(lhs, ScalarFunction<T>(v));
   }
}

#endif
