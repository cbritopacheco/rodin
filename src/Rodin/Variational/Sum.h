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
#include "MatrixCoefficient.h"

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

      private:
         std::unique_ptr<ScalarFunctionBase> m_lhs;
         std::unique_ptr<ScalarFunctionBase> m_rhs;
   };
   Sum(const ScalarFunctionBase&, const ScalarFunctionBase&)
      -> Sum<ScalarFunctionBase, ScalarFunctionBase>;

   Sum<ScalarFunctionBase, ScalarFunctionBase>
   operator+(const ScalarFunctionBase& lhs, const ScalarFunctionBase& rhs);

   /**
    * @brief %Sum of two MatrixCoefficientBase instances.
    */
   template <>
   class Sum<MatrixCoefficientBase, MatrixCoefficientBase>
      : public MatrixCoefficientBase
   {
      public:
         Sum(const MatrixCoefficientBase& lhs, const MatrixCoefficientBase& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Sum(const Sum& other)
            :  MatrixCoefficientBase(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Sum(Sum&& other)
            :  MatrixCoefficientBase(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         MatrixCoefficientBase& getLHS()
         {
            return *m_lhs;
         }

         MatrixCoefficientBase& getRHS()
         {
            return *m_rhs;
         }

         const MatrixCoefficientBase& getLHS() const
         {
            return *m_lhs;
         }

         const MatrixCoefficientBase& getRHS() const
         {
            return *m_rhs;
         }

         int getRows() const override;

         int getColumns() const override;

         void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override;

         Sum* copy() const noexcept override
         {
            return new Sum(*this);
         }

      private:
         std::unique_ptr<MatrixCoefficientBase> m_lhs;
         std::unique_ptr<MatrixCoefficientBase> m_rhs;
   };
   Sum(const MatrixCoefficientBase&, const MatrixCoefficientBase&)
      -> Sum<MatrixCoefficientBase, MatrixCoefficientBase>;
   Sum<MatrixCoefficientBase, MatrixCoefficientBase>
   operator+(const MatrixCoefficientBase& lhs, const MatrixCoefficientBase& rhs);

   template <ShapeFunctionSpaceType Space>
   class Sum<ShapeFunctionBase<Space>, ShapeFunctionBase<Space>>
      : public ShapeFunctionBase<Space>
   {
      public:
         Sum(const ShapeFunctionBase<Space>& lhs, const ShapeFunctionBase<Space>& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {
            assert(lhs.getRoot().getUUID() == rhs.getRoot().getUUID());
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

         ShapeFunctionBase<Space>& getRoot() override
         {
            return getRHS().getRoot();
         }

         const ShapeFunctionBase<Space>& getRoot() const override
         {
            return getRHS().getRoot();
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
}

#endif
