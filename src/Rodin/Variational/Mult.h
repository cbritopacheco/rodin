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
#include "Rodin/Variational/GridFunction.h"
#include "Rodin/Variational/ScalarCoefficient.h"
#include "Rodin/Variational/VectorCoefficient.h"
#include "Rodin/Variational/MatrixCoefficient.h"
#include "Rodin/Variational/TestFunction.h"
#include "Rodin/Variational/TrialFunction.h"
#include "FormLanguage/Base.h"
#include "ForwardDecls.h"

namespace Rodin::Variational
{
   /**
    * @brief Multiplication of two ScalarCoefficientBase instances.
    */
   template <>
   class Mult<ScalarCoefficientBase, ScalarCoefficientBase>
      : public ScalarCoefficientBase
   {
      public:
         Mult(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Mult(const Mult& other)
            :  ScalarCoefficientBase(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Mult(Mult&& other)
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

         Mult* copy() const noexcept override
         {
            return new Mult(*this);
         }
      private:
         std::unique_ptr<ScalarCoefficientBase> m_lhs;
         std::unique_ptr<ScalarCoefficientBase> m_rhs;
   };
   Mult(const ScalarCoefficientBase&, const ScalarCoefficientBase&)
      -> Mult<ScalarCoefficientBase, ScalarCoefficientBase>;

   Mult<ScalarCoefficientBase, ScalarCoefficientBase>
   operator*(const ScalarCoefficientBase& lhs, const ScalarCoefficientBase& rhs);

   template <class T>
   std::enable_if_t<std::is_arithmetic_v<T>, Mult<ScalarCoefficientBase, ScalarCoefficientBase>>
   operator*(T lhs, const ScalarCoefficientBase& rhs)
   {
      return Mult(ScalarCoefficient(lhs), rhs);
   }

   template <class T>
   std::enable_if_t<std::is_arithmetic_v<T>, Mult<ScalarCoefficientBase, ScalarCoefficientBase>>
   operator*(const ScalarCoefficientBase& rhs, T lhs)
   {
      return Mult(rhs, ScalarCoefficient(lhs));
   }

   /**
    * @brief Multiplication of ScalarCoefficientBase and VectorCoefficientBase.
    */
   template <>
   class Mult<ScalarCoefficientBase, VectorCoefficientBase>
      : public VectorCoefficientBase
   {
      public:
         Mult(const ScalarCoefficientBase& lhs, const VectorCoefficientBase& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Mult(const Mult& other)
            :  VectorCoefficientBase(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Mult(Mult&& other)
            :  VectorCoefficientBase(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         ScalarCoefficientBase& getLHS()
         {
            return *m_lhs;
         }

         VectorCoefficientBase& getRHS()
         {
            return *m_rhs;
         }

         const ScalarCoefficientBase& getLHS() const
         {
            return *m_lhs;
         }

         const VectorCoefficientBase& getRHS() const
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
         std::unique_ptr<ScalarCoefficientBase> m_lhs;
         std::unique_ptr<VectorCoefficientBase> m_rhs;
   };
   Mult(const ScalarCoefficientBase&, const VectorCoefficientBase&)
      -> Mult<ScalarCoefficientBase, VectorCoefficientBase>;

   Mult<ScalarCoefficientBase, VectorCoefficientBase>
   operator*(const ScalarCoefficientBase& lhs, const VectorCoefficientBase& rhs);

   Mult<ScalarCoefficientBase, VectorCoefficientBase>
   operator*(const VectorCoefficientBase& lhs, const ScalarCoefficientBase& rhs);

   template <class T>
   std::enable_if_t<std::is_arithmetic_v<T>, Mult<ScalarCoefficientBase, VectorCoefficientBase>>
   operator*(T lhs, const VectorCoefficientBase& rhs)
   {
      return Mult(ScalarCoefficient(lhs), rhs);
   }

   template <class T>
   std::enable_if_t<std::is_arithmetic_v<T>, Mult<ScalarCoefficientBase, VectorCoefficientBase>>
   operator*(const VectorCoefficientBase& lhs, T rhs)
   {
      return Mult(ScalarCoefficient(rhs), lhs);
   }

   /**
    * @brief Multiplication of ScalarCoefficientBase and MatrixCoefficientBase.
    */
   template <>
   class Mult<ScalarCoefficientBase, MatrixCoefficientBase>
      : public MatrixCoefficientBase
   {
      public:
         Mult(const ScalarCoefficientBase& lhs, const MatrixCoefficientBase& rhs)
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
         std::unique_ptr<ScalarCoefficientBase> m_lhs;
         std::unique_ptr<MatrixCoefficientBase> m_rhs;
   };
   Mult(const ScalarCoefficientBase&, const MatrixCoefficientBase&)
      -> Mult<ScalarCoefficientBase, MatrixCoefficientBase>;

   Mult<ScalarCoefficientBase, MatrixCoefficientBase>
   operator*(const ScalarCoefficientBase& lhs, const MatrixCoefficientBase& rhs);

   Mult<ScalarCoefficientBase, MatrixCoefficientBase>
   operator*(const MatrixCoefficientBase& lhs, const ScalarCoefficientBase& rhs);

   template <class T>
   std::enable_if_t<std::is_arithmetic_v<T>, Mult<ScalarCoefficientBase, MatrixCoefficientBase>>
   operator*(T lhs, const MatrixCoefficientBase& rhs)
   {
      return Mult(ScalarCoefficient(lhs), rhs);
   }

   template <class T>
   std::enable_if_t<std::is_arithmetic_v<T>, Mult<ScalarCoefficientBase, MatrixCoefficientBase>>
   operator*(const MatrixCoefficientBase& lhs, T rhs)
   {
      return Mult(ScalarCoefficient(rhs), lhs);
   }

   template <ShapeFunctionSpaceType Space>
   class Mult<ScalarCoefficientBase, ShapeFunctionBase<Space>>
      : public ShapeFunctionBase<Space>
   {
      public:
         Mult(const ScalarCoefficientBase& lhs, const ShapeFunctionBase<Space>& rhs)
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
         std::unique_ptr<ScalarCoefficientBase> m_lhs;
         std::unique_ptr<ShapeFunctionBase<Space>> m_rhs;
   };
   template <ShapeFunctionSpaceType Space>
   Mult(const ScalarCoefficientBase&, const ShapeFunctionBase<Space>&)
      -> Mult<ScalarCoefficientBase, ShapeFunctionBase<Space>>;

   template <ShapeFunctionSpaceType Space>
   Mult<ScalarCoefficientBase, ShapeFunctionBase<Space>>
   operator*(const ScalarCoefficientBase& lhs, const ShapeFunctionBase<Space>& rhs)
   {
      return Mult(lhs, rhs);
   }

   template <ShapeFunctionSpaceType Space>
   Mult<ScalarCoefficientBase, ShapeFunctionBase<Space>>
   operator*(const ShapeFunctionBase<Space>& lhs, const ScalarCoefficientBase& rhs)
   {
      return Mult(rhs, lhs);
   }

   template <class T, ShapeFunctionSpaceType Space>
   std::enable_if_t<std::is_arithmetic_v<T>, Mult<ScalarCoefficientBase, ShapeFunctionBase<Space>>>
   operator*(T lhs, const ShapeFunctionBase<Space>& rhs)
   {
      return Mult(ScalarCoefficient(lhs), rhs);
   }

   template <class T, ShapeFunctionSpaceType Space>
   std::enable_if_t<std::is_arithmetic_v<T>, Mult<ScalarCoefficientBase, ShapeFunctionBase<Space>>>
   operator*(const ShapeFunctionBase<Space>& lhs, T rhs)
   {
      return Mult(lhs, ScalarCoefficient(rhs));
   }
}

#endif
