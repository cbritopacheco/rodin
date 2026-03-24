/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Fields/FirstPiolaKirchhoffStress.h
 * @brief First Piola-Kirchhoff stress tensor computation.
 *
 * Computes the first Piola-Kirchhoff stress @f$ \mathbf{P} @f$ from a
 * hyperelastic constitutive law:
 * @f[
 *   \mathbf{P} = \frac{\partial W}{\partial \mathbf{F}}
 * @f]
 */
#ifndef RODIN_SOLID_FIELDS_FIRSTPIOLAKIRCHHOFFSTRESS_H
#define RODIN_SOLID_FIELDS_FIRSTPIOLAKIRCHHOFFSTRESS_H

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialMatrix.h"

#include "Rodin/Solid/Kinematics/KinematicState.h"
#include "Rodin/Solid/Local/ConstitutivePoint.h"
#include "Rodin/Solid/Constitutive/HyperElasticLaw.h"

namespace Rodin::Solid
{
  /**
   * @brief Computes the first Piola-Kirchhoff stress from a constitutive law.
   *
   * ## Usage
   * @code
   * NeoHookean law(lambda, mu);
   * NeoHookean::Cache cache;
   * KinematicState state(d);
   * state.setDisplacementGradient(gradU);
   * law.setCache(cache, state);
   *
   * FirstPiolaKirchhoffStress<NeoHookean> pk1(law);
    * Math::SpatialMatrix<Real> P;
   * pk1.getFirstPiolaKirchhoffStress(P, cache, state);
   * @endcode
   *
   * @tparam LawDerived The hyperelastic law type
   */
  template <class LawDerived>
  class FirstPiolaKirchhoffStress
  {
    public:
      using LawType = LawDerived;
      using CacheType = typename LawType::Cache;

      /**
       * @brief Constructs the stress evaluator.
       * @param law Reference to the constitutive law
       */
      FirstPiolaKirchhoffStress(const LawType& law)
        : m_law(law)
      {}

      FirstPiolaKirchhoffStress(const FirstPiolaKirchhoffStress&) = default;
      FirstPiolaKirchhoffStress(FirstPiolaKirchhoffStress&&) = default;

      /**
       * @brief Computes @f$ \mathbf{P} = \partial W / \partial \mathbf{F} @f$.
       * @param[out] P Output stress tensor
       * @param[in] cache Precomputed law cache
       * @param[in] cp Constitutive point
       */
      void getFirstPiolaKirchhoffStress(
          Math::SpatialMatrix<Real>& P,
          const CacheType& cache,
          const ConstitutivePoint& cp) const
      {
        m_law.getFirstPiolaKirchhoffStress(P, cache, cp);
      }

      /// @brief Gets the constitutive law.
      const LawType& getLaw() const { return m_law; }

    private:
      const LawType& m_law;
  };
}

#endif
