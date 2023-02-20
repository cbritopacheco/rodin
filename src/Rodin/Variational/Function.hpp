/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_FUNCTION_HPP
#define RODIN_VARIATIONAL_FUNCTION_HPP

#include "Function.h"

#include "Traits.h"

namespace Rodin::Variational
{
  template <class Derived>
  inline
  constexpr
  RangeType FunctionBase<Derived>::getRangeType() const
  {
    using R = typename FormLanguage::Traits<FunctionBase<Derived>>::RangeType;
    if constexpr (std::is_same_v<R, Boolean>)
    {
      return RangeType::Boolean;
    }
    else if constexpr (std::is_same_v<R, Scalar>)
    {
      return RangeType::Scalar;
    }
    else
    {
      assert(false);
      return RangeType::Scalar;
    }
  }
}

#endif

