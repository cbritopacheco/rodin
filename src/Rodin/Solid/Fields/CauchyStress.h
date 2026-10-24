/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Fields/CauchyStress.h
 * @brief Cauchy (true) stress tensor computation from the first Piola-Kirchhoff stress.
 *
 * Computes the Cauchy stress tensor via the Piola transform:
 * @f[
 *   \boldsymbol{\sigma} = \frac{1}{J} \mathbf{P} \mathbf{F}^T
 * @f]
 * where @f$ \mathbf{P} @f$ is the first Piola-Kirchhoff stress,
 * @f$ \mathbf{F} @f$ is the deformation gradient, and
 * @f$ J = \det(\mathbf{F}) @f$.
 */
#ifndef RODIN_SOLID_FIELDS_CAUCHYSTRESS_H
#define RODIN_SOLID_FIELDS_CAUCHYSTRESS_H

#include "Rodin/Types.h"
#include "Rodin/Math/Matrix.h"
#include "Rodin/Math/SpatialMatrix.h"

#include "Rodin/Solid/Kinematics/KinematicState.h"
#include "Rodin/Solid/Local/ConstitutivePoint.h"
#include "Rodin/Solid/Constitutive/HyperElasticLaw.h"

namespace Rodin::Solid
{
  /**
   * @brief Computes the Cauchy stress from a constitutive law.
   *
   * ## Usage
   * @code
   * NeoHookean law(lambda, mu);
   * KinematicState state(d);
   * state.setDisplacementGradient(gradU);
   * ConstitutivePoint point(state);
   * law.setConstitutivePoint(point);
   *
   * CauchyStress<NeoHookean> cauchy(law);
   * Math::SpatialMatrix<Real> sigma;
   * cauchy.getCauchyStress(sigma);
   * @endcode
   *
   * @tparam LawDerived The hyperelastic law type
   */
  template <class LawDerived>
  class CauchyStress
  {
    public:
      /// @brief Constitutive law type.
      using LawType = LawDerived;

      /**
       * @brief Constructs the Cauchy stress evaluator.
       * @param law Reference to a prepared constitutive law
       */
      CauchyStress(const LawType& law)
        : m_law(law)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      CauchyStress(const CauchyStress& other) = default;

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      CauchyStress(CauchyStress&& other) = default;

      /**
       * @brief Computes @f$ \boldsymbol{\sigma} = \tfrac{1}{J} \mathbf{P} \mathbf{F}^T @f$.
       * @param[out] sigma Output Cauchy stress tensor
       */
      void getCauchyStress(Math::SpatialMatrix<Real>& sigma) const
      {
        const auto& state = m_law.getKinematicState();
        Math::SpatialMatrix<Real> P;
        m_law.getFirstPiolaKirchhoffStress(P);
        const auto& F = state.getDeformationGradient();
        const Real J = state.getJacobian();
        sigma = (1.0 / J) * P * F.transpose();
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
