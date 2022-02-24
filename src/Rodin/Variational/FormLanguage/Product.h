/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FORMLANGUAGE_PRODUCT_H
#define RODIN_VARIATIONAL_FORMLANGUAGE_PRODUCT_H

#include <memory>
#include <type_traits>

#include "Rodin/Alert.h"
#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Variational/ScalarCoefficient.h"
#include "Rodin/Variational/VectorCoefficient.h"
#include "Rodin/Variational/MatrixCoefficient.h"
#include "Rodin/Variational/TestFunction.h"
#include "Rodin/Variational/TrialFunction.h"

#include "ForwardDecls.h"
#include "Base.h"

namespace Rodin::Variational::FormLanguage
{
   /**
    * @brief Product between instances of Lhs and Rhs
    * @tparam Lhs Left-hand side operand type
    * @tparam Rhs Right-hand side operand type
    */
   template <class Lhs, class Rhs>
   class Product : public Base
   {
      static_assert(std::is_base_of_v<Base, Lhs>,
            "Lhs must be derived from FormLanguage::Base");
      static_assert(std::is_base_of_v<Base, Rhs>,
            "Rhs must be derived from FormLanguage::Base");

      public:
         Product(const Lhs& lhs, const Rhs& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Product(const Product& other)
            : m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Product(Product&&) = default;

         virtual ~Product() = default;

         Lhs& getLHS()
         {
            assert(m_lhs);
            return *m_lhs;
         }

         Rhs& getRHS()
         {
            assert(m_rhs);
            return *m_rhs;
         }

         const Lhs& getLHS() const
         {
            assert(m_lhs);
            return *m_lhs;
         }

         const Rhs& getRHS() const
         {
            assert(m_lhs);
            return *m_rhs;
         }

         Product* copy() const noexcept override
         {
            return new Product(*this);
         }

      private:
         std::unique_ptr<Lhs> m_lhs;
         std::unique_ptr<Rhs> m_rhs;
   };

   /**
    * @brief Product between two instances of ScalarCoefficientBase
    */
   template <>
   class Product<ScalarCoefficientBase, ScalarCoefficientBase>
      : public ScalarCoefficientBase
   {
      public:
         Product(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Product(const Product& other)
            :  ScalarCoefficientBase(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Product(Product&& other)
            :  ScalarCoefficientBase(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         ScalarCoefficientBase& getLHS()
         {
            return *m_lhs;
         }

         ScalarCoefficientBase& getRHS()
         {
            return *m_rhs;
         }

         const ScalarCoefficientBase& getLHS() const
         {
            return *m_lhs;
         }

         const ScalarCoefficientBase& getRHS() const
         {
            return *m_rhs;
         }

         double getValue(
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override
         {
            return getLHS().getValue(trans, ip) * getRHS().getValue(trans, ip);
         }

         Product* copy() const noexcept override
         {
            return new Product(*this);
         }
      private:
         std::unique_ptr<ScalarCoefficientBase> m_lhs;
         std::unique_ptr<ScalarCoefficientBase> m_rhs;
   };

   /**
    * @brief Product between instances of ScalarCoefficientBase and MatrixCoefficientBase
    */
   template <>
   class Product<ScalarCoefficientBase, MatrixCoefficientBase>
      : public MatrixCoefficientBase
   {
      public:
         Product(const ScalarCoefficientBase& lhs, const MatrixCoefficientBase& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Product(const Product& other)
            :  MatrixCoefficientBase(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Product(Product&& other)
            :  MatrixCoefficientBase(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         ScalarCoefficientBase& getLHS()
         {
            return *m_lhs;
         }

         MatrixCoefficientBase& getRHS()
         {
            return *m_rhs;
         }

         const ScalarCoefficientBase& getLHS() const
         {
            return *m_lhs;
         }

         const MatrixCoefficientBase& getRHS() const
         {
            return *m_rhs;
         }

         int getRows() const override
         {
            return getRHS().getRows();
         }

         int getColumns() const override
         {
            return getRHS().getColumns();
         }

         void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) override
         {
            getRHS().getValue(value, trans, ip);
            value *= getLHS().getValue(trans, ip);
         }

         Product* copy() const noexcept override
         {
            return new Product(*this);
         }
      private:
         std::unique_ptr<ScalarCoefficientBase> m_lhs;
         std::unique_ptr<MatrixCoefficientBase> m_rhs;
   };

   template <ShapeFunctionSpaceType Space>
   class Product<ScalarCoefficientBase, ShapeFunctionBase<Space>>
      : public ShapeFunctionBase<Space>
   {
      public:
         Product(const ScalarCoefficientBase& lhs, const ShapeFunctionBase<Space>& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Product(const Product& other)
            :  ShapeFunctionBase<Space>(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Product(Product&& other)
            :  ShapeFunctionBase<Space>(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         ScalarCoefficientBase& getLHS()
         {
            return *m_lhs;
         }

         ShapeFunctionBase<Space>& getRHS()
         {
            return *m_rhs;
         }

         const ScalarCoefficientBase& getLHS() const
         {
            return *m_lhs;
         }

         const ShapeFunctionBase<Space>& getRHS() const
         {
            return *m_rhs;
         }

         size_t getRows(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            return getRHS().getRows(fe, trans);
         }

         size_t getColumns(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            return getRHS().getColumns(fe, trans);
         }

         void getOperator(
              const mfem::FiniteElement& fe,
              mfem::ElementTransformation& trans,
              mfem::DenseMatrix& op) const override
         {
            getRHS().getOperator(fe, trans, op);
            const mfem::IntegrationPoint& ip = trans.GetIntPoint();
            op *= getLHS().getValue(trans, ip);
         }

         FiniteElementSpaceBase& getFiniteElementSpace() override
         {
            return getRHS().getFiniteElementSpace();
         }

         const FiniteElementSpaceBase& getFiniteElementSpace() const override
         {
            return getRHS().getFiniteElementSpace();
         }

         Product* copy() const noexcept override
         {
            return new Product(*this);
         }
      private:
         std::unique_ptr<ScalarCoefficientBase> m_lhs;
         std::unique_ptr<ShapeFunctionBase<Space>> m_rhs;
   };

   template <ShapeFunctionSpaceType Space>
   class Product<VectorCoefficientBase, ShapeFunctionBase<Space>>
      : public ShapeFunctionBase<Space>
   {
      public:
         Product(const VectorCoefficientBase& lhs, const ShapeFunctionBase<Space>& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Product(const Product& other)
            :  ShapeFunctionBase<Space>(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Product(Product&& other)
            :  ShapeFunctionBase<Space>(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         VectorCoefficientBase& getLHS()
         {
            return *m_lhs;
         }

         ShapeFunctionBase<Space>& getRHS()
         {
            return *m_rhs;
         }

         const VectorCoefficientBase& getLHS() const
         {
            return *m_lhs;
         }

         const ShapeFunctionBase<Space>& getRHS() const
         {
            return *m_rhs;
         }

         size_t getRows(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            assert(false);
         }

         size_t getColumns(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            assert(false);
         }

         void getOperator(
              const mfem::FiniteElement& fe,
              mfem::ElementTransformation& trans,
              mfem::DenseMatrix& op) const override
         {
            assert(false);
         }

         FiniteElementSpaceBase& getFiniteElementSpace() override
         {
            return getRHS().getFiniteElementSpace();
         }

         const FiniteElementSpaceBase& getFiniteElementSpace() const override
         {
            return getRHS().getFiniteElementSpace();
         }

         Product* copy() const noexcept override
         {
            return new Product(*this);
         }
      private:
         std::unique_ptr<VectorCoefficientBase> m_lhs;
         std::unique_ptr<ShapeFunctionBase<Space>> m_rhs;
   };

   template <ShapeFunctionSpaceType Space>
   Product<ScalarCoefficientBase, ShapeFunctionBase<Space>>
   operator*(const ScalarCoefficientBase& lhs, const ShapeFunctionBase<Space>& rhs)
   {
      return Product(lhs, rhs);
   }

   template <ShapeFunctionSpaceType Space>
   Product<VectorCoefficientBase, ShapeFunctionBase<Space>>
   operator*(const VectorCoefficientBase& lhs, const ShapeFunctionBase<Space>& rhs)
   {
      return Product(lhs, rhs);
   }

   template <class FEC>
   Product<ScalarCoefficientBase, GridFunction<FEC>>
   operator*(const ScalarCoefficientBase& lhs, const GridFunction<FEC>& rhs)
   {
      return Product(lhs, ScalarCoefficient(rhs));
   }

   Product<ScalarCoefficientBase, ScalarCoefficientBase>
   operator*(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs);

   Product<ScalarCoefficientBase, MatrixCoefficientBase>
   operator*(const ScalarCoefficientBase& lhs, const MatrixCoefficientBase& rhs);

   Product<ScalarCoefficientBase, MatrixCoefficientBase>
   operator*(const MatrixCoefficientBase& lhs, const ScalarCoefficientBase& rhs);

   Product<ShapeFunctionBase<Trial>, ShapeFunctionBase<Test>>
   operator*(const ShapeFunctionBase<Trial>& lhs, const ShapeFunctionBase<Test>& rhs);
}

#endif
