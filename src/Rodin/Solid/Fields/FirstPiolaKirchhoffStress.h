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
   * # Usage
   * @code
   * NeoHookean law(lambda, mu);
   * KinematicState state(d);
   * state.setDisplacementGradient(gradU);
   * ConstitutivePoint point(state);
   * law.setConstitutivePoint(point);
   *
   * FirstPiolaKirchhoffStress<NeoHookean> pk1(law);
   * Math::SpatialMatrix<Real> P;
   * pk1.getFirstPiolaKirchhoffStress(P);
   * @endcode
   *
   * @tparam LawDerived The hyperelastic law type
   */
  template <class LawDerived>
  class FirstPiolaKirchhoffStress
  {
    public:
      /// @brief Constitutive law type.
      using LawType = LawDerived;

      /**
       * @brief Constructs the stress evaluator.
       * @param law Reference to a prepared constitutive law
       */
      FirstPiolaKirchhoffStress(const LawType& law)
        : m_law(law)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      FirstPiolaKirchhoffStress(const FirstPiolaKirchhoffStress& other) = default;

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      FirstPiolaKirchhoffStress(FirstPiolaKirchhoffStress&& other) = default;

      /**
       * @brief Computes @f$ \mathbf{P} = \partial W / \partial \mathbf{F} @f$.
       * @param[out] P Output stress tensor
       */
      void getFirstPiolaKirchhoffStress(Math::SpatialMatrix<Real>& P) const
      {
        m_law.getFirstPiolaKirchhoffStress(P);
      }

      /**
       * @brief Gets the constitutive law.
       * @returns The constitutive law.
       */
      const LawType& getLaw() const
      {
        return m_law;
      }

    private:
      const LawType& m_law;
  };
}

#endif
