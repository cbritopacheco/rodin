/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Input.h
 * @brief Extensible input injection mechanism for constitutive evaluation.
 *
 * An Input populates a ConstitutivePoint with auxiliary data at each
 * quadrature point during integration.  This is the injection site for
 * material-specific fields (fiber directions, activation parameters,
 * region-wise material properties, etc.) without hard-coding them into the
 * integrator or the constitutive law.
 *
 * ## Usage with CRTP
 *
 * @code
 * struct MyInput : Solid::Input<MyInput>
 * {
 *   void populate(Solid::ConstitutivePoint& cp) const
 *   {
 *     // e.g., interpolate fiber direction from a GridFunction
 *     Math::SpatialVector<Real> fiber = ...;
 *     cp.set<Solid::Tags::FiberDirection>(fiber);
 *   }
 * };
 *
 * Solid::InternalForce force(law, v);
 * force.setInput(MyInput{});
 * @endcode
 *
 * ## Usage with lambda / std::function
 *
 * @code
 * force.setInput([&](Solid::ConstitutivePoint& cp) {
 *   cp.set<Solid::Tags::Activation>(activationValue);
 * });
 * @endcode
 */
#ifndef RODIN_SOLID_LOCAL_INPUT_H
#define RODIN_SOLID_LOCAL_INPUT_H

#include <functional>

#include "ConstitutivePoint.h"

namespace Rodin::Solid
{
  /**
   * @brief CRTP base class for inputs.
   *
   * Derived classes implement a `populate(ConstitutivePoint&)` method that
   * injects auxiliary data into the ConstitutivePoint at each quadrature
   * point during assembly.
   *
   * @tparam Derived The concrete input type (CRTP)
   */
  template <class Derived>
  class Input
  {
    public:
      /**
       * @brief Populates a ConstitutivePoint with auxiliary data.
       *
       * Called by integrators at each quadrature point after the kinematic
       * state and geometric context have been set.
       *
       * @param[in,out] cp The constitutive point to populate
       */
      void populate(ConstitutivePoint& cp) const
      {
        static_cast<const Derived&>(*this).populate(cp);
      }

    protected:
      Input() = default;
      Input(const Input&) = default;
      Input(Input&&) = default;
      Input& operator=(const Input&) = default;
      Input& operator=(Input&&) = default;
  };

  /// @brief Type-erased callable for input injection into ConstitutivePoint.
  using InputFunction = std::function<void(ConstitutivePoint&)>;
}

#endif
