/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file MaterialTangent.h
 * @brief Material tangent operator for Newton-Raphson linearization.
 *
 * Evaluates the bilinear form arising from the linearization of the
 * internal virtual work:
 * @f[
 *   a(\delta\mathbf{u}, \mathbf{v})
 *     = \int_\Omega D\mathbf{P}[\nabla \delta\mathbf{u}]
 *       : \nabla \mathbf{v} \, dX
 * @f]
 * where @f$ D\mathbf{P}[\cdot] @f$ denotes the directional derivative
 * of the first Piola-Kirchhoff stress.
 *
 * This class evaluates the material tangent action at a given
 * kinematic state and perturbation, providing the values needed for
 * the tangent stiffness bilinear form assembly.
 */
#ifndef RODIN_SOLID_INTEGRATORS_MATERIALTANGENT_H
#define RODIN_SOLID_INTEGRATORS_MATERIALTANGENT_H

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"

#include "Rodin/Solid/Kinematics/KinematicState.h"
#include "Rodin/Solid/Constitutive/HyperElasticLaw.h"

namespace Rodin::Solid
{
  /**
   * @brief Evaluates the material tangent action @f$ D\mathbf{P}[\delta\mathbf{F}] @f$.
   *
   * This is a helper class for assembling the tangent stiffness matrix
   * in a Newton-Raphson iteration. The user provides the perturbation
   * @f$ \delta\mathbf{F} @f$ (gradient of trial function increment) and
   * obtains the linearized stress increment.
   *
   * ## Usage
   * @code
   * NeoHookean law(lambda, mu);
   * NeoHookean::Cache cache;
   * KinematicState state(d);
   * state.setDisplacementGradient(gradU).update();
   * law.initializeCache(cache, state);
   *
   * MaterialTangent<NeoHookean> tangent(law);
   * Math::Matrix<Real> dP;
   * tangent.evaluate(dP, cache, state, dF);
   * @endcode
   *
   * @tparam LawDerived The hyperelastic law type
   */
  template <class LawDerived>
  class MaterialTangent
  {
    public:
      using LawType = LawDerived;
      using CacheType = typename LawType::Cache;

      /**
       * @brief Constructs the material tangent evaluator.
       * @param law Reference to the constitutive law
       */
      MaterialTangent(const LawType& law)
        : m_law(law)
      {}

      MaterialTangent(const MaterialTangent&) = default;
      MaterialTangent(MaterialTangent&&) = default;

      /**
       * @brief Evaluates the material tangent action.
       * @param[out] dP Output linearized stress increment
       * @param[in] cache Precomputed law cache
       * @param[in] state Current kinematic state
       * @param[in] dF Perturbation of the deformation gradient
       */
      void evaluate(
          Math::Matrix<Real>& dP,
          const CacheType& cache,
          const KinematicState& state,
          const Math::Matrix<Real>& dF) const
      {
        m_law.materialTangentAction(dP, cache, state, dF);
      }

      /// @brief Gets the constitutive law.
      const LawType& getLaw() const { return m_law; }

    private:
      const LawType& m_law;
  };
}

#endif
