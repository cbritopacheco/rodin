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

   template <ShapeFunctionSpaceType Space>
   class Mult<FunctionBase, ShapeFunctionBase<Space>>
      : public ShapeFunctionBase<Space>
   {
      public:
         Mult(const FunctionBase& lhs, const ShapeFunctionBase<Space>& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {
            if (lhs.getRangeType() != RangeType::Scalar)
               UnexpectedRangeTypeException(RangeType::Scalar, lhs.getRangeType()).raise();
         }

         Mult(const Mult& other)
            :  ShapeFunctionBase<Space>(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Mult(Mult&& other)
            :  ShapeFunctionBase<Space>(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         const ShapeFunctionBase<Space>& getLeaf() const override
         {
            return m_rhs->getLeaf();
         }

         int getRows() const override
         {
            return m_rhs->getRangeShape().height();
         }

         int getColumns() const override
         {
            return m_rhs->getRangeShape().width();
         }

         int getDOFs(
               const mfem::FiniteElement& fe,
               const mfem::ElementTransformation& trans) const override
         {
            return m_rhs->getDOFs(fe, trans);
         }

         void getOperator(
               DenseBasisOperator& op,
               const mfem::FiniteElement& fe,
               mfem::ElementTransformation& trans) const override
         {
            switch (m_lhs->getRangeType())
            {
               case RangeType::Scalar:
               {
                  m_rhs->getOperator(op, fe, trans);
                  mfem::DenseMatrix v;
                  m_lhs->getValue(v, trans, trans.GetIntPoint());
                  op *= v(0, 0);
                  break;
               }
               default:
                  assert(false); // Unimplemented
            }
         }

         FiniteElementSpaceBase& getFiniteElementSpace() override
         {
            return m_rhs->getFiniteElementSpace();
         }

         const FiniteElementSpaceBase& getFiniteElementSpace() const override
         {
            return m_rhs->getFiniteElementSpace();
         }

         Mult* copy() const noexcept override
         {
            return new Mult(*this);
         }
      private:
         std::unique_ptr<FunctionBase> m_lhs;
         std::unique_ptr<ShapeFunctionBase<Space>> m_rhs;
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

   template <ShapeFunctionSpaceType Space>
   Mult<FunctionBase, ShapeFunctionBase<Space>>
   operator*(const ShapeFunctionBase<Space>& lhs, const FunctionBase& rhs)
   {
      return Mult(rhs, lhs);
   }

   template <class T, ShapeFunctionSpaceType Space>
   std::enable_if_t<std::is_arithmetic_v<T>, Mult<FunctionBase, ShapeFunctionBase<Space>>>
   operator*(const ShapeFunctionBase<Space>& lhs, T v)
   {
      return Mult(ScalarFunction(v), lhs);
   }

   template <class T, ShapeFunctionSpaceType Space>
   std::enable_if_t<std::is_arithmetic_v<T>, Mult<FunctionBase, ShapeFunctionBase<Space>>>
   operator*(T v, const ShapeFunctionBase<Space>& rhs)
   {
      return Mult(ScalarFunction(v), rhs);
   }
}

#endif
