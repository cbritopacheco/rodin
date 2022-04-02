#ifndef RODIN_VARIATIONAL_POW_H
#define RODIN_VARIATIONAL_POW_H

#include <cmath>

#include "ScalarFunction.h"

namespace Rodin::Variational
{
   template <class T>
   class Pow : public ScalarFunctionBase
   {
      static_assert(std::is_arithmetic_v<T>, "T must be an arithmetic type");
      public:
         Pow(const ScalarFunctionBase& s, T p)
            : m_s(s.copy()),
              m_p(p)
         {}

         Pow(const Pow& other)
            : ScalarFunctionBase(other),
              m_s(other.m_s->copy()),
              m_p(other.m_p)
         {}

         Pow(Pow&& other)
            : ScalarFunctionBase(std::move(other)),
              m_s(std::move(other.m_s)),
              m_p(other.m_p)
         {}

         double getValue(
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) const override
         {
            return std::pow(m_s->getValue(trans, ip), m_p);
         }

         Pow* copy() const noexcept override
         {
            return new Pow(*this);
         }
      private:
         std::unique_ptr<ScalarFunctionBase> m_s;
         T m_p;
   };
}

#endif
