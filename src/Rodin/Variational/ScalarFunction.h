/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_SCALARFUNCTION_H
#define RODIN_VARIATIONAL_SCALARFUNCTION_H

#include "ForwardDecls.h"

#include "Function.h"

namespace Rodin::FormLanguage
{
  template <class Scalar, class Derived>
  struct Traits<Variational::ScalarFunctionBase<Scalar, Derived>>
  {
    using ScalarType = Scalar;
    using DerivedType = Derived;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup RealFunctionSpecializations RealFunction Template Specializations
   * @brief Template specializations of the RealFunction class.
   * @see RealFunction
   */

  template <class Scalar, class Derived>
  class ScalarFunctionBase
    : public FunctionBase<ScalarFunctionBase<Scalar, Derived>>
  {
    public:
      using ScalarType = Scalar;

      using Parent = FunctionBase<ScalarFunctionBase<ScalarType, Derived>>;

      using Parent::traceOf;

      using Parent::operator();


      ScalarFunctionBase() = default;

      ScalarFunctionBase(const ScalarFunctionBase& other)
        : Parent(other)
      {}

      ScalarFunctionBase(ScalarFunctionBase&& other)
        : Parent(std::move(other))
      {}

      virtual ~ScalarFunctionBase() = default;

      constexpr
      decltype(auto) getValue(const Geometry::Point& p) const
      {
        return static_cast<const Derived&>(*this).getValue(p);
      }

      virtual ScalarFunctionBase* copy() const noexcept override = 0;
  };
}

#endif

