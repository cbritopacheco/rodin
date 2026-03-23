/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file InternalForce.h
 * @brief Internal force vector computation for hyperelastic formulations.
 *
 * Assembles the nonlinear residual vector (internal force) contribution:
 * @f[
 *   R_{\text{int}}(\mathbf{v}) = \int_\Omega \mathbf{P}(\mathbf{u}) : \nabla \mathbf{v} \, dX
 * @f]
 * where @f$ \mathbf{P} @f$ is the first Piola-Kirchhoff stress and
 * @f$ \mathbf{v} @f$ is the test function.
 *
 * This class evaluates the stress at each quadrature point using the
 * supplied constitutive law and kinematic state, then provides the values
 * needed for the linear form assembly in Rodin's variational framework.
 */
#ifndef RODIN_SOLID_INTEGRATORS_INTERNALFORCE_H
#define RODIN_SOLID_INTEGRATORS_INTERNALFORCE_H

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"

#include "Rodin/Solid/Kinematics/KinematicState.h"
#include "Rodin/Solid/Constitutive/HyperElasticLaw.h"

namespace Rodin::Solid
{
  /**
   * @brief Evaluates the internal force density @f$ \mathbf{P} : \nabla \mathbf{v} @f$
   * at a given quadrature point.
   *
   * This is a helper class that provides the stress tensor at a given
   * kinematic state for use in linear form assembly. It does not assemble
   * by itself; instead, the user constructs the residual linear form using
   * Rodin's existing @c Integral facilities with the stress tensor as a
   * matrix-valued coefficient.
   *
   * ## Usage
   * @code
   * NeoHookean law(lambda, mu);
   * NeoHookean::Cache cache;
   * KinematicState state(d);
   * state.setDisplacementGradient(gradU).update();
   * law.initializeCache(cache, state);
   *
   * InternalForce<NeoHookean> force(law);
   * Math::Matrix<Real> P;
   * force.evaluate(P, cache, state);
   * // P now holds the first Piola-Kirchhoff stress for residual assembly
   * @endcode
   *
   * @tparam LawDerived The hyperelastic law type
   */
  template <class LawDerived>
  class InternalForce
  {
    public:
      using LawType = LawDerived;
      using CacheType = typename LawType::Cache;

      /**
       * @brief Constructs the internal force evaluator.
       * @param law Reference to the constitutive law
       */
      InternalForce(const LawType& law)
        : m_law(law)
      {}

      InternalForce(const InternalForce&) = default;
      InternalForce(InternalForce&&) = default;

      /**
       * @brief Evaluates the first Piola-Kirchhoff stress.
       * @param[out] P Output stress tensor
       * @param[in] cache Precomputed law cache
       * @param[in] state Current kinematic state
       */
      void evaluate(
          Math::Matrix<Real>& P,
          const CacheType& cache,
          const KinematicState& state) const
      {
        m_law.firstPiolaKirchhoffStress(P, cache, state);
      }

      /// @brief Gets the constitutive law.
      const LawType& getLaw() const { return m_law; }

    private:
      const LawType& m_law;
  };
}

#endif
