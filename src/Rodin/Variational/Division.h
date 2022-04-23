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
#include "VectorFunction.h"

namespace Rodin::Variational
{
   /**
    * @brief Division of VectorFunctionBase by ScalarFunctionBase.
    */
   template <>
   class Division<ScalarFunctionBase, ScalarFunctionBase>
      : public ScalarFunctionBase
   {
      public:
         Division(const ScalarFunctionBase& lhs, const ScalarFunctionBase& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Division(const Division& other)
            :  ScalarFunctionBase(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Division(Division&& other)
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
            return getLHS().getValue(trans, ip) / getRHS().getValue(trans, ip);
         }

         Division* copy() const noexcept override
         {
            return new Division(*this);
         }
      private:
         std::unique_ptr<ScalarFunctionBase> m_lhs;
         std::unique_ptr<ScalarFunctionBase> m_rhs;
   };
   Division(const ScalarFunctionBase&, const ScalarFunctionBase&)
      -> Division<ScalarFunctionBase, ScalarFunctionBase>;

   Division<ScalarFunctionBase, ScalarFunctionBase>
   operator/(const ScalarFunctionBase& lhs, const ScalarFunctionBase& rhs);

   template <class T>
   std::enable_if_t<std::is_arithmetic_v<T>,
      Division<ScalarFunctionBase, ScalarFunctionBase>>
   operator/(const ScalarFunctionBase& lhs, T rhs)
   {
      return Division(lhs, ScalarFunction(rhs));
   }

   /**
    * @brief Division of VectorFunctionBase by ScalarFunctionBase.
    */
   template <>
   class Division<VectorFunctionBase, ScalarFunctionBase>
      : public VectorFunctionBase
   {
      public:
         Division(const VectorFunctionBase& lhs, const ScalarFunctionBase& rhs)
            : m_lhs(lhs.copy()), m_rhs(rhs.copy())
         {}

         Division(const Division& other)
            :  VectorFunctionBase(other),
               m_lhs(other.m_lhs->copy()), m_rhs(other.m_rhs->copy())
         {}

         Division(Division&& other)
            :  VectorFunctionBase(std::move(other)),
               m_lhs(std::move(other.m_lhs)), m_rhs(std::move(other.m_rhs))
         {}

         VectorFunctionBase& getLHS()
         {
            return *m_lhs;
         }

         ScalarFunctionBase& getRHS()
         {
            return *m_rhs;
         }

         const VectorFunctionBase& getLHS() const
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
         std::unique_ptr<VectorFunctionBase> m_lhs;
         std::unique_ptr<ScalarFunctionBase> m_rhs;
   };
   Division(const VectorFunctionBase&, const ScalarFunctionBase&)
      -> Division<VectorFunctionBase, ScalarFunctionBase>;

   Division<VectorFunctionBase, ScalarFunctionBase>
   operator/(const VectorFunctionBase& lhs, const ScalarFunctionBase& rhs);

   template <class T>
   std::enable_if_t<std::is_arithmetic_v<T>,
      Division<VectorFunctionBase, ScalarFunctionBase>>
   operator/(const VectorFunctionBase& lhs, T rhs)
   {
      return Division(lhs, ScalarFunction(rhs));
   }
}
#endif
