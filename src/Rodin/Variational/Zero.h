/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_ZERO_H
#define RODIN_VARIATIONAL_ZERO_H

#include <cmath>

#include <Rodin/Math/Common.h>

#include "ForwardDecls.h"

#include "RangeShape.h"
#include "Function.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup ZeroSpecializations Zero Template Specializations
   * @brief Template specializations of the Zero class.
   * @see Zero
   */

  template <>
  class Zero<Scalar> final
    : public ScalarFunctionBase<Zero<Scalar>>
  {
    public:
      /// Parent class
      using Parent = ScalarFunctionBase<Zero<Scalar>>;

      /**
       * @brief Constructs the Zeroer object
       * @param[in] s Base value
       * @param[in] p Zeroer
       */
      Zero() {}

      Zero(const Zero& other)
        : Parent(std::move(other))
      {}

      Zero(Zero&& other)
        : Parent(std::move(other))
      {}

      inline
      constexpr
      Zero& traceOf(Geometry::Attribute attrs)
      {
        return *this;
      }

      inline
      constexpr
      auto getValue(const Geometry::Point&) const
      {
        return Scalar(0);
      }

      inline Zero* copy() const noexcept override
      {
        return new Zero(*this);
      }
  };

  /**
   * @brief CTAD for Zero.
   */
  Zero() -> Zero<Scalar>;
}

#endif

