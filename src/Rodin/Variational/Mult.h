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
#include "Rodin/FormLanguage/Base.h"

#include "ForwardDecls.h"
#include "Function.h"
#include "ScalarFunction.h"
#include "ShapeFunction.h"

#include "Exceptions.h"

namespace Rodin::Variational
{
   template <class Lhs, class Rhs>
   class Mult
   {};

   /**
    * @brief Multiplication of two FunctionBase instances.
    */
   template <>
   class Mult<FunctionBase, FunctionBase> : public FunctionBase
   {
      public:
         Mult(const FunctionBase& lhs, const FunctionBase& rhs);

         Mult(const Mult& other);

         Mult(Mult&& other);

         RangeShape getRangeShape() const override;

         Mult& traceOf(const std::set<int>& attrs) override;

         void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) const override;

         Mult* copy() const noexcept override
         {
            return new Mult(*this);
         }
      private:
         std::unique_ptr<FunctionBase> m_lhs;
         std::unique_ptr<FunctionBase> m_rhs;
   };
   Mult(const FunctionBase&, const FunctionBase&) -> Mult<FunctionBase, FunctionBase>;

   Mult<FunctionBase, FunctionBase>
   operator*(const FunctionBase& lhs, const FunctionBase& rhs);

   template <class T>
   std::enable_if_t<std::is_arithmetic_v<T>, Mult<FunctionBase, FunctionBase>>
   operator*(T lhs, const FunctionBase& rhs)
   {
      return Mult(ScalarFunction(lhs), rhs);
   }

   template <class T>
   std::enable_if_t<std::is_arithmetic_v<T>, Mult<FunctionBase, FunctionBase>>
   operator*(const FunctionBase& lhs, T rhs)
   {
      return Mult(lhs, ScalarFunction(rhs));
   }

   /**
    * @brief Left Multiplication by a FunctionBase of a ShapeFunctionBase
    *
    * Represents the following expression:
    * @f[
    *    f A(u)
    * @f]
    * where @f$ f @f$ is the function (scalar, vector or matrix valued), and
    * $A(u)$ is the shape operator (scalar, vector, or matrix valued).
    */
   template <ShapeFunctionSpaceType Space>
   class Mult<FunctionBase, ShapeFunctionBase<Space>>
      : public ShapeFunctionBase<Space>
   {
      public:
         constexpr
         Mult(const FunctionBase& lhs, const ShapeFunctionBase<Space>& rhs)
            : m_f(lhs.copy()), m_u(rhs.copy())
         {}

         constexpr
         Mult(const Mult& other)
            :  ShapeFunctionBase<Space>(other),
               m_f(other.m_f->copy()), m_u(other.m_u->copy())
         {}

         constexpr
         Mult(Mult&& other)
            :  ShapeFunctionBase<Space>(std::move(other)),
               m_f(std::move(other.m_f)), m_u(std::move(other.m_u))
         {}

         const ShapeFunctionBase<Space>& getLeaf() const override
         {
            return m_u->getLeaf();
         }

         int getRows() const override
         {
            return m_u->getRangeShape().height();
         }

         int getColumns() const override
         {
            return m_u->getRangeShape().width();
         }

         int getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            return m_u->getDOFs(fe, trans);
         }

         void getOperator(
               DenseBasisOperator& op,
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip,
               ShapeComputator& compute) const override
         {
            switch (m_f->getRangeType())
            {
               case RangeType::Scalar:
               {
                  m_u->getOperator(op, fe, trans, ip, compute);
                  mfem::DenseMatrix v;
                  m_f->getValue(v, trans, ip);
                  op *= v(0, 0);
                  break;
               }
               default:
                  assert(false); // Unimplemented
            }
         }

         FiniteElementSpaceBase& getFiniteElementSpace() override
         {
            return m_u->getFiniteElementSpace();
         }

         const FiniteElementSpaceBase& getFiniteElementSpace() const override
         {
            return m_u->getFiniteElementSpace();
         }

         virtual Mult* copy() const noexcept override
         {
            return new Mult(*this);
         }
      private:
         std::unique_ptr<FunctionBase> m_f;
         std::unique_ptr<ShapeFunctionBase<Space>> m_u;
   };
   template <ShapeFunctionSpaceType Space>
   Mult(const FunctionBase&, const ShapeFunctionBase<Space>&)
      -> Mult<FunctionBase, ShapeFunctionBase<Space>>;

   template <ShapeFunctionSpaceType Space>
   Mult<FunctionBase, ShapeFunctionBase<Space>>
   operator*(const FunctionBase& lhs, const ShapeFunctionBase<Space>& rhs)
   {
      return Mult(lhs, rhs);
   }

   template <class T, ShapeFunctionSpaceType Space>
   std::enable_if_t<std::is_arithmetic_v<T>, Mult<FunctionBase, ShapeFunctionBase<Space>>>
   operator*(T v, const ShapeFunctionBase<Space>& rhs)
   {
      return Mult(ScalarFunction(v), rhs);
   }

   /* ||--- OPTIMIZATIONS ----------------------------------------------------
    * Mult<FunctionBase, ShapeFunctionBase<Space>>
    * ---------------------------------------------------------------------->>
    */

   /**
    * @internal
    * @brief Left Multiplication of a ShapeFunction by a FunctionBase
    *
    * Represents the following expression:
    * @f[
    *    f u
    * @f]
    * where @f$ f @f$ is a function (scalar or or matrix valued).
    */
   template <ShapeFunctionSpaceType Space, class FES>
   class Mult<FunctionBase, ShapeFunction<FES, Space>>
      : public Mult<FunctionBase, ShapeFunctionBase<Space>>
   {
      public:
         constexpr
         Mult(const FunctionBase& lhs, const ShapeFunction<FES, Space>& rhs)
            : Mult<FunctionBase, ShapeFunctionBase<Space>>(lhs, rhs)
         {
            if (lhs.getRangeType() != RangeType::Scalar && lhs.getRangeType() != RangeType::Matrix)
               UnexpectedRangeTypeException({RangeType::Scalar, RangeType::Matrix}, lhs.getRangeType()).raise();
         }

         constexpr
         Mult(const Mult& other)
            : Mult<FunctionBase, ShapeFunctionBase<Space>>(other)
         {}

         constexpr
         Mult(Mult&& other)
            : Mult<FunctionBase, ShapeFunctionBase<Space>>(std::move(other))
         {}

         Mult* copy() const noexcept override
         {
            return new Mult(*this);
         }
   };
   template <ShapeFunctionSpaceType Space, class FES>
   Mult(const FunctionBase&, const ShapeFunction<FES, Space>& rhs)
      -> Mult<FunctionBase, ShapeFunction<FES, Space>>;

   template <ShapeFunctionSpaceType Space, class FES>
   Mult<FunctionBase, ShapeFunctionBase<Space>>
   operator*(const FunctionBase& lhs, const ShapeFunction<FES, Space>& rhs)
   {
      return Mult(lhs, rhs);
   }

   /**
    * @internal
    * @brief Left Multiplication of the gradient of a ShapeFunction by a FunctionBase
    *
    * Represents the following expression:
    * @f[
    *    f \nabla u
    * @f]
    * where @f$ f @f$ is a function (scalar or matrix valued).
    */
   template <ShapeFunctionSpaceType Space, class FES>
   class Mult<FunctionBase, Grad<ShapeFunction<FES, Space>>>
      : public Mult<FunctionBase, ShapeFunctionBase<Space>>
   {
      public:
         constexpr
         Mult(const FunctionBase& lhs, const Grad<ShapeFunction<FES, Space>>& rhs)
            : Mult<FunctionBase, ShapeFunctionBase<Space>>(lhs, rhs)
         {
            if (lhs.getRangeType() != RangeType::Scalar && lhs.getRangeType() != RangeType::Matrix)
               UnexpectedRangeTypeException({RangeType::Scalar, RangeType::Matrix}, lhs.getRangeType()).raise();
         }

         constexpr
         Mult(const Mult& other)
            : Mult<FunctionBase, ShapeFunctionBase<Space>>(other)
         {}

         constexpr
         Mult(Mult&& other)
            : Mult<FunctionBase, ShapeFunctionBase<Space>>(std::move(other))
         {}

         Mult* copy() const noexcept override
         {
            return new Mult(*this);
         }
   };
   template <ShapeFunctionSpaceType Space, class FES>
   Mult(const FunctionBase&, const Grad<ShapeFunction<FES, Space>>& rhs)
      -> Mult<FunctionBase, Grad<ShapeFunction<FES, Space>>>;

   template <ShapeFunctionSpaceType Space, class FES>
   Mult<FunctionBase, ShapeFunctionBase<Space>>
   operator*(const FunctionBase& lhs, const Grad<ShapeFunction<FES, Space>>& rhs)
   {
      return Mult(lhs, rhs);
   }

   /* <<-- OPTIMIZATIONS -----------------------------------------------------
    * Mult<FunctionBase, ShapeFunctionBase<Space>>
    * ----------------------------------------------------------------------||
    */
}

#endif
