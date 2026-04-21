/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Input.h
 * @brief Extensible input injection mechanism for fluid constitutive evaluation.
 */
#ifndef RODIN_FLUID_LOCAL_INPUT_H
#define RODIN_FLUID_LOCAL_INPUT_H

#include <functional>

#include "FlowPoint.h"

namespace Rodin::Fluid
{
  /**
   * @brief CRTP base class for FlowPoint input providers.
   *
   * Derived classes implement:
   * @code
   * void populate(FlowPoint& fp) const;
   * @endcode
   * to inject auxiliary data at each quadrature point.
   */
  template <class Derived>
  class Input
  {
    public:
      void populate(FlowPoint& fp) const
      {
        static_cast<const Derived&>(*this).populate(fp);
      }

    protected:
      Input() = default;
      Input(const Input&) = default;
      Input(Input&&) = default;
      Input& operator=(const Input&) = default;
      Input& operator=(Input&&) = default;
  };

  /// @brief Type-erased callable for FlowPoint input injection.
  using InputFunction = std::function<void(FlowPoint&)>;
}

#endif
