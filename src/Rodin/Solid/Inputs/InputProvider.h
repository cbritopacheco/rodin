/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file InputProvider.h
 * @brief Extensible input injection mechanism for constitutive evaluation.
 *
 * An InputProvider populates a ConstitutivePoint with auxiliary data at each
 * quadrature point during integration.  This is the injection site for
 * material-specific fields (fiber directions, activation parameters,
 * region-wise material properties, etc.) without hard-coding them into the
 * integrator or the constitutive law.
 *
 * ## Usage with CRTP
 *
 * @code
 * struct MyProvider : Solid::InputProvider<MyProvider>
 * {
 *   void populate(Solid::ConstitutivePoint& cp) const
 *   {
 *     // e.g., interpolate fiber direction from a GridFunction
 *     struct FiberTag {};
 *     Math::SpatialVector<Real> fiber = ...;
 *     cp.setAuxiliaryData<FiberTag>(fiber);
 *   }
 * };
 *
 * Solid::InternalForce force(law, v);
 * force.setInputProvider(MyProvider{});
 * @endcode
 *
 * ## Usage with lambda / std::function
 *
 * @code
 * force.setInputProvider([&](Solid::ConstitutivePoint& cp) {
 *   struct ActivationTag {};
 *   cp.setAuxiliaryData<ActivationTag>(activationValue);
 * });
 * @endcode
 */
#ifndef RODIN_SOLID_INPUTS_INPUTPROVIDER_H
#define RODIN_SOLID_INPUTS_INPUTPROVIDER_H

#include <functional>

#include "ConstitutivePoint.h"

namespace Rodin::Solid
{
  /**
   * @brief CRTP base class for input providers.
   *
   * Derived classes implement a `populate(ConstitutivePoint&)` method that
   * injects auxiliary data into the ConstitutivePoint at each quadrature
   * point during assembly.
   *
   * @tparam Derived The concrete input provider type (CRTP)
   */
  template <class Derived>
  class InputProvider
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
      InputProvider() = default;
      InputProvider(const InputProvider&) = default;
      InputProvider(InputProvider&&) = default;
      InputProvider& operator=(const InputProvider&) = default;
      InputProvider& operator=(InputProvider&&) = default;
  };

  /// @brief Type-erased callable for input injection into ConstitutivePoint.
  using InputProviderFunction = std::function<void(ConstitutivePoint&)>;
}

#endif
