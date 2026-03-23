/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file HyperElasticLaw.h
 * @brief Base class for hyperelastic constitutive laws.
 *
 * Provides the CRTP base class for all hyperelastic constitutive laws.
 * Each derived law must provide:
 * - A `Cache` struct for precomputed data
 * - `setCache(Cache&, const KinematicState&)` to populate the cache
 * - `getFirstPiolaKirchhoffStress(SpatialMatrix&, const Cache&, const KinematicState&)`
 *   to compute @f$ \mathbf{P} = \partial W / \partial \mathbf{F} @f$
 * - `getMaterialTangent(SpatialMatrix&, const Cache&, const KinematicState&, const SpatialMatrix&)`
 *   to compute @f$ D\mathbf{P}[\delta\mathbf{F}] @f$
 */
#ifndef RODIN_SOLID_CONSTITUTIVE_HYPERELASTICLAW_H
#define RODIN_SOLID_CONSTITUTIVE_HYPERELASTICLAW_H

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"

#include "Rodin/Solid/Kinematics/KinematicState.h"

namespace Rodin::Solid
{
  /**
   * @brief CRTP base class for hyperelastic constitutive laws.
   *
   * A hyperelastic material is characterized by a stored energy density
   * function @f$ W(\mathbf{F}) @f$ from which all stress measures derive.
   *
   * Each derived class implements:
   * - **Cache setup**: Precomputes law-specific data from the kinematic state.
   * - **First Piola-Kirchhoff stress**: @f$ \mathbf{P} = \frac{\partial W}{\partial \mathbf{F}} @f$
   * - **Material tangent**: Directional derivative
   *   @f$ D\mathbf{P}[\delta\mathbf{F}] = \frac{\partial^2 W}{\partial \mathbf{F}^2} : \delta\mathbf{F} @f$
   *
   * @tparam Derived The derived constitutive law type (CRTP)
   */
  template <class Derived>
  class HyperElasticLaw
  {
    public:
      /**
       * @brief Sets the cache from the kinematic state.
       * @param[out] cache Cache to populate
       * @param[in] state Current kinematic state
       */
      template <class Cache>
      void setCache(Cache& cache, const KinematicState& state) const
      {
        static_cast<const Derived&>(*this).setCache(cache, state);
      }

      /**
       * @brief Computes the first Piola-Kirchhoff stress tensor.
       *
       * @f[
       *   \mathbf{P} = \frac{\partial W}{\partial \mathbf{F}}
       * @f]
       *
       * @param[out] P Output stress tensor
       * @param[in] cache Precomputed cache
       * @param[in] state Current kinematic state
       */
      template <class Cache>
      void getFirstPiolaKirchhoffStress(
          Math::Matrix<Real>& P,
          const Cache& cache,
          const KinematicState& state) const
      {
        static_cast<const Derived&>(*this).getFirstPiolaKirchhoffStress(P, cache, state);
      }

      /**
       * @brief Computes the material tangent on a perturbation.
       *
       * Evaluates the directional derivative of the first Piola-Kirchhoff
       * stress:
       * @f[
       *   D\mathbf{P}[\delta\mathbf{F}]
       *     = \frac{\partial^2 W}{\partial \mathbf{F}^2} : \delta\mathbf{F}
       * @f]
       *
       * @param[out] dP Output tangent stress increment
       * @param[in] cache Precomputed cache
       * @param[in] state Current kinematic state
       * @param[in] dF Perturbation of the deformation gradient
       */
      template <class Cache>
      void getMaterialTangent(
          Math::Matrix<Real>& dP,
          const Cache& cache,
          const KinematicState& state,
          const Math::Matrix<Real>& dF) const
      {
        static_cast<const Derived&>(*this).getMaterialTangent(dP, cache, state, dF);
      }

      /**
       * @brief Computes the strain energy density.
       *
       * @f[
       *   W(\mathbf{F})
       * @f]
       *
       * @param[in] cache Precomputed cache
       * @param[in] state Current kinematic state
       * @returns Strain energy density value
       */
      template <class Cache>
      Real getStrainEnergyDensity(
          const Cache& cache,
          const KinematicState& state) const
      {
        return static_cast<const Derived&>(*this).getStrainEnergyDensity(cache, state);
      }

    protected:
      HyperElasticLaw() = default;
      HyperElasticLaw(const HyperElasticLaw&) = default;
      HyperElasticLaw(HyperElasticLaw&&) = default;
      HyperElasticLaw& operator=(const HyperElasticLaw&) = default;
      HyperElasticLaw& operator=(HyperElasticLaw&&) = default;
  };
}

#endif
