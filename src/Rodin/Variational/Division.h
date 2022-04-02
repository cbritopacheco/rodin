/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_DIVISION_H
#define RODIN_VARIATIONAL_DIVISION_H

#include "ForwardDecls.h"
#include "ScalarFunction.h"
#include "VectorCoefficient.h"

namespace Rodin::Variational
{
   /**
    * @brief Division of VectorCoefficientBase by ScalarFunctionBase.
    */
   template <>
   class Division<VectorCoefficientBase, ScalarFunctionBase>
      : public VectorCoefficientBase
   {
      public:
         Division(const VectorCoefficientBase& lhs, const ScalarFunctionBase& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Division(const Division& other)
            :  VectorCoefficientBase(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Division(Division&& other)
            :  VectorCoefficientBase(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         VectorCoefficientBase& getLHS()
         {
            return *m_lhs;
         }

         ScalarFunctionBase& getRHS()
         {
            return *m_rhs;
         }

         const VectorCoefficientBase& getLHS() const
         {
            return *m_lhs;
         }

         const ScalarFunctionBase& getRHS() const
         {
            return *m_rhs;
         }

         int getDimension() const override
         {
            return getLHS().getDimension();
         }

         void getValue(
               mfem::Vector& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const override
         {
            getLHS().getValue(value, trans, ip);
            value /= getRHS().getValue(trans, ip);
         }

         Division* copy() const noexcept override
         {
            return new Division(*this);
         }
      private:
         std::unique_ptr<VectorCoefficientBase> m_lhs;
         std::unique_ptr<ScalarFunctionBase> m_rhs;
   };
   Division(const VectorCoefficientBase&, const ScalarFunctionBase&)
      -> Division<VectorCoefficientBase, ScalarFunctionBase>;

   Division<VectorCoefficientBase, ScalarFunctionBase>
   operator/(const VectorCoefficientBase& lhs, const ScalarFunctionBase& rhs);

   template <class T>
   std::enable_if_t<std::is_arithmetic_v<T>,
      Division<VectorCoefficientBase, ScalarFunctionBase>>
   operator*(const VectorCoefficientBase& lhs, T rhs)
   {
      return Division(lhs, ScalarFunction(rhs));
   }
}
#endif
