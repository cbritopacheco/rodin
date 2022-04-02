/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_MULT_H
#define RODIN_VARIATIONAL_MULT_H

#include <memory>
#include <type_traits>

#include "Rodin/Alert.h"
#include "FormLanguage/Base.h"

#include "ForwardDecls.h"

#include "GridFunction.h"
#include "ScalarFunction.h"
#include "VectorFunction.h"
#include "MatrixCoefficient.h"
#include "TestFunction.h"
#include "TrialFunction.h"

namespace Rodin::Variational
{
   /**
    * @brief Multiplication of two ScalarFunctionBase instances.
    */
   template <>
   class Mult<ScalarFunctionBase, ScalarFunctionBase>
      : public ScalarFunctionBase
   {
      public:
         Mult(const ScalarFunctionBase& lhs, const ScalarFunctionBase& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Mult(const Mult& other)
            :  ScalarFunctionBase(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Mult(Mult&& other)
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
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override
         {
            return getLHS().getValue(trans, ip) * getRHS().getValue(trans, ip);
         }

         Mult* copy() const noexcept override
         {
            return new Mult(*this);
         }
      private:
         std::unique_ptr<ScalarFunctionBase> m_lhs;
         std::unique_ptr<ScalarFunctionBase> m_rhs;
   };
   Mult(const ScalarFunctionBase&, const ScalarFunctionBase&)
      -> Mult<ScalarFunctionBase, ScalarFunctionBase>;

   Mult<ScalarFunctionBase, ScalarFunctionBase>
   operator*(const ScalarFunctionBase& lhs, const ScalarFunctionBase& rhs);

   template <class T>
   std::enable_if_t<std::is_arithmetic_v<T>, Mult<ScalarFunctionBase, ScalarFunctionBase>>
   operator*(T lhs, const ScalarFunctionBase& rhs)
   {
      return Mult(ScalarFunction(lhs), rhs);
   }

   template <class T>
   std::enable_if_t<std::is_arithmetic_v<T>, Mult<ScalarFunctionBase, ScalarFunctionBase>>
   operator*(const ScalarFunctionBase& rhs, T lhs)
   {
      return Mult(rhs, ScalarFunction(lhs));
   }

   /**
    * @brief Multiplication of ScalarFunctionBase and VectorFunctionBase.
    */
   template <>
   class Mult<ScalarFunctionBase, VectorFunctionBase>
      : public VectorFunctionBase
   {
      public:
         Mult(const ScalarFunctionBase& lhs, const VectorFunctionBase& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Mult(const Mult& other)
            :  VectorFunctionBase(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Mult(Mult&& other)
            :  VectorFunctionBase(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         ScalarFunctionBase& getLHS()
         {
            return *m_lhs;
         }

         VectorFunctionBase& getRHS()
         {
            return *m_rhs;
         }

         const ScalarFunctionBase& getLHS() const
         {
            return *m_lhs;
         }

         const VectorFunctionBase& getRHS() const
         {
            return *m_rhs;
         }

         int getDimension() const override
         {
            return getRHS().getDimension();
         }

         void getValue(
               mfem::Vector& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override
         {
            getRHS().getValue(value, trans, ip);
            value *= getLHS().getValue(trans, ip);
         }

         Mult* copy() const noexcept override
         {
            return new Mult(*this);
         }
      private:
         std::unique_ptr<ScalarFunctionBase> m_lhs;
         std::unique_ptr<VectorFunctionBase> m_rhs;
   };
   Mult(const ScalarFunctionBase&, const VectorFunctionBase&)
      -> Mult<ScalarFunctionBase, VectorFunctionBase>;

   Mult<ScalarFunctionBase, VectorFunctionBase>
   operator*(const ScalarFunctionBase& lhs, const VectorFunctionBase& rhs);

   Mult<ScalarFunctionBase, VectorFunctionBase>
   operator*(const VectorFunctionBase& lhs, const ScalarFunctionBase& rhs);

   template <class T>
   std::enable_if_t<std::is_arithmetic_v<T>, Mult<ScalarFunctionBase, VectorFunctionBase>>
   operator*(T lhs, const VectorFunctionBase& rhs)
   {
      return Mult(ScalarFunction(lhs), rhs);
   }

   template <class T>
   std::enable_if_t<std::is_arithmetic_v<T>, Mult<ScalarFunctionBase, VectorFunctionBase>>
   operator*(const VectorFunctionBase& lhs, T rhs)
   {
      return Mult(ScalarFunction(rhs), lhs);
   }

   /**
    * @brief Multiplication of ScalarFunctionBase and MatrixCoefficientBase.
    */
   template <>
   class Mult<ScalarFunctionBase, MatrixCoefficientBase>
      : public MatrixCoefficientBase
   {
      public:
         Mult(const ScalarFunctionBase& lhs, const MatrixCoefficientBase& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Mult(const Mult& other)
            :  MatrixCoefficientBase(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Mult(Mult&& other)
            :  MatrixCoefficientBase(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         ScalarFunctionBase& getLHS()
         {
            return *m_lhs;
         }

         MatrixCoefficientBase& getRHS()
         {
            return *m_rhs;
         }

         const ScalarFunctionBase& getLHS() const
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
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override
         {
            getRHS().getValue(value, trans, ip);
            value *= getLHS().getValue(trans, ip);
         }

         Mult* copy() const noexcept override
         {
            return new Mult(*this);
         }
      private:
         std::unique_ptr<ScalarFunctionBase> m_lhs;
         std::unique_ptr<MatrixCoefficientBase> m_rhs;
   };
   Mult(const ScalarFunctionBase&, const MatrixCoefficientBase&)
      -> Mult<ScalarFunctionBase, MatrixCoefficientBase>;

   Mult<ScalarFunctionBase, MatrixCoefficientBase>
   operator*(const ScalarFunctionBase& lhs, const MatrixCoefficientBase& rhs);

   Mult<ScalarFunctionBase, MatrixCoefficientBase>
   operator*(const MatrixCoefficientBase& lhs, const ScalarFunctionBase& rhs);

   template <class T>
   std::enable_if_t<std::is_arithmetic_v<T>, Mult<ScalarFunctionBase, MatrixCoefficientBase>>
   operator*(T lhs, const MatrixCoefficientBase& rhs)
   {
      return Mult(ScalarFunction(lhs), rhs);
   }

   template <class T>
   std::enable_if_t<std::is_arithmetic_v<T>, Mult<ScalarFunctionBase, MatrixCoefficientBase>>
   operator*(const MatrixCoefficientBase& lhs, T rhs)
   {
      return Mult(ScalarFunction(rhs), lhs);
   }

   template <ShapeFunctionSpaceType Space>
   class Mult<ScalarFunctionBase, ShapeFunctionBase<Space>>
      : public ShapeFunctionBase<Space>
   {
      public:
         Mult(const ScalarFunctionBase& lhs, const ShapeFunctionBase<Space>& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Mult(const Mult& other)
            :  ShapeFunctionBase<Space>(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Mult(Mult&& other)
            :  ShapeFunctionBase<Space>(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         ScalarFunctionBase& getLHS()
         {
            return *m_lhs;
         }

         ShapeFunctionBase<Space>& getRHS()
         {
            return *m_rhs;
         }

         const ScalarFunctionBase& getLHS() const
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
            return getRHS().getRows(fe, trans);
         }

         int getColumns(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            return getRHS().getColumns(fe, trans);
         }

         int getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            return getRHS().getDOFs(fe, trans);
         }

         std::unique_ptr<Rank3Operator> getOperator(
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans) const override
         {
            auto result = getRHS().getOperator(fe, trans);
            (*result) *= getLHS().getValue(trans, trans.GetIntPoint());
            return result;
         }

         FiniteElementSpaceBase& getFiniteElementSpace() override
         {
            return getRHS().getFiniteElementSpace();
         }

         const FiniteElementSpaceBase& getFiniteElementSpace() const override
         {
            return getRHS().getFiniteElementSpace();
         }

         Mult* copy() const noexcept override
         {
            return new Mult(*this);
         }
      private:
         std::unique_ptr<ScalarFunctionBase> m_lhs;
         std::unique_ptr<ShapeFunctionBase<Space>> m_rhs;
   };
   template <ShapeFunctionSpaceType Space>
   Mult(const ScalarFunctionBase&, const ShapeFunctionBase<Space>&)
      -> Mult<ScalarFunctionBase, ShapeFunctionBase<Space>>;

   template <ShapeFunctionSpaceType Space>
   Mult<ScalarFunctionBase, ShapeFunctionBase<Space>>
   operator*(const ScalarFunctionBase& lhs, const ShapeFunctionBase<Space>& rhs)
   {
      return Mult(lhs, rhs);
   }

   template <ShapeFunctionSpaceType Space>
   Mult<ScalarFunctionBase, ShapeFunctionBase<Space>>
   operator*(const ShapeFunctionBase<Space>& lhs, const ScalarFunctionBase& rhs)
   {
      return Mult(rhs, lhs);
   }

   template <class T, ShapeFunctionSpaceType Space>
   std::enable_if_t<std::is_arithmetic_v<T>, Mult<ScalarFunctionBase, ShapeFunctionBase<Space>>>
   operator*(T lhs, const ShapeFunctionBase<Space>& rhs)
   {
      return Mult(ScalarFunction(lhs), rhs);
   }

   template <class T, ShapeFunctionSpaceType Space>
   std::enable_if_t<std::is_arithmetic_v<T>, Mult<ScalarFunctionBase, ShapeFunctionBase<Space>>>
   operator*(const ShapeFunctionBase<Space>& lhs, T rhs)
   {
      return Mult(lhs, ScalarFunction(rhs));
   }
}

#endif
