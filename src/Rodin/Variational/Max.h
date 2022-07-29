/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_MAX_H
#define RODIN_VARIATIONAL_MAX_H

#include <cmath>

#include "ForwardDecls.h"
#include "Function.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
   /**
    * @brief Represent the max function.
    */
   class Max : public ScalarFunctionBase
   {
      public:
         Max(const FunctionBase& a, double b);

         Max(double a, const FunctionBase& b);

         /**
          * @bref Constructs the power object
          * @param[in] s Base value
          * @param[in] p Power
          */
         Max(const FunctionBase& a, const FunctionBase& b);

         Max(const Max& other);

         Max(Max&& other);

         Max& traceOf(const std::set<int>& attrs) override;

         double getValue(
               mfem::ElementTransformation& trans,
               const mfem::IntegrationPoint& ip) const override;

         Max* copy() const noexcept override;

      private:
         std::unique_ptr<FunctionBase> m_a;
         std::unique_ptr<FunctionBase> m_b;
   };
}

#endif
