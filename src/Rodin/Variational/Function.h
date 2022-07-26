/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FUNCTION_H
#define RODIN_VARIATIONAL_FUNCTION_H

#include <tuple>
#include <mfem.hpp>

namespace Rodin::Variational
{
   enum class RangeType
   {
      Scalar,
      Vector,
      Matrix
   };

   class FunctionBase
   {
      public:
         virtual RangeType getRangeType() const = 0;

         virtual std::tuple<int, int> getRangeShape() const = 0;

         virtual void getValue(
               mfem::DenseMatrix& value,
               mfem::ElementTransformation& trans, const mfem::IntegrationPoint& ip) const = 0;
   };
}

#endif
