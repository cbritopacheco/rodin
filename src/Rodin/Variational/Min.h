/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_MIN_H
#define RODIN_VARIATIONAL_MIN_H

#include <cmath>

#include "ForwardDecls.h"
#include "Function.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
   /**
    * @brief Represent the min function.
    */
   class Min : public ScalarFunctionBase
   {
      public:
         Min(const FunctionBase& a, double b);

         Min(double a, const FunctionBase& b);

         Min(const FunctionBase& a, const FunctionBase& b);

         Min(const Min& other);

         Min(Min&& other);

         Min& traceOf(const std::set<int>& attrs) override;

         double getValue(
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) const override;

         Min* copy() const noexcept override;

      private:
         std::unique_ptr<FunctionBase> m_a;
         std::unique_ptr<FunctionBase> m_b;
   };
}

#endif
