/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLID_CONSTITUTIVE_KELVINVOIGT_H
#define RODIN_SOLID_CONSTITUTIVE_KELVINVOIGT_H

#include <cmath>
#include <utility>
#include "Rodin/Alert.h"
#include "Rodin/Solid/Local/ConstitutivePoint.h"

namespace Rodin::Solid
{
  /**
   * @brief Finite-strain Kelvin–Voigt viscosity added to an elastic material.
   *
   * The dissipation potential is @f$\Phi=\nu\dot E:\dot E/2@f$,
   * with @f$\dot E=\operatorname{sym}(F^T G)@f$ and @f$G=\nabla_Xv@f$.
   * Thus @f$P=P_e+\nu F\dot E@f$. The tangent includes both the change in
   * @f$F@f$ and the time-integrator relation @f$\delta G=c_v\delta F@f$.
   *
   * ## Architecture
   *
   * The wrapped material owns elastic kinematics and its stored energy.
   * Binding captures the velocity gradient and its displacement derivative
   * through constitutive tags. Stress and tangent share this rate state.
   * Stored energy excludes dissipation; getDissipationRate() reports its power.
   *
   * @tparam Law Elastic material providing the standard stress/tangent interface.
   */
  template <class Law>
  class KelvinVoigt
  {
    public:
      /**
       * @brief Constructs a viscous material with an elastic branch.
       * @param law Elastic material, owned by value.
       * @param viscosity Nonnegative coefficient in stress times time units.
       */
      KelvinVoigt(Law law, Real viscosity)
        : m_law(std::move(law)),
          m_viscosity(viscosity),
          m_velocityCoefficient(0),
          m_bound(false)
      {
        if (!std::isfinite(viscosity) || viscosity < 0)
          Alert::Exception() << "Invalid solid viscosity" << Alert::Raise;
      }

      /**
       * @brief Captures elastic and rate inputs without borrowing their storage.
       * @param point Constitutive point with velocity tags when viscosity is positive.
       * @returns This material.
       */
      KelvinVoigt& setConstitutivePoint(const ConstitutivePoint& point)
      {
        m_bound = false;
        const auto dimension = point.getKinematicState().getDimension();
        Math::SpatialMatrix<Real> gradient(dimension, dimension);
        gradient.setZero();
        Real coefficient = 0;
        if (m_viscosity > 0)
        {
          if (!point.has<Tags::VelocityGradient>() ||
            !point.has<Tags::VelocityGradientDerivative>())
          {
            Alert::Exception() << "Solid viscosity requires velocity-gradient inputs"
                               << Alert::Raise;
          }
          gradient = point.get<Tags::VelocityGradient>();
          coefficient = point.get<Tags::VelocityGradientDerivative>();
          if (gradient.rows() != dimension || gradient.cols() != dimension ||
            !gradient.getData().topLeftCorner(dimension, dimension).allFinite() ||
            !std::isfinite(coefficient) || coefficient < 0)
          {
            Alert::Exception() << "Invalid solid velocity-gradient inputs"
                               << Alert::Raise;
          }
        }
        m_law.setConstitutivePoint(point);
        m_velocityGradient = gradient;
        m_velocityCoefficient = coefficient;
        const auto& f = m_law.getKinematicState().getDeformationGradient();
        m_strainRate = Real(0.5) * (f.transpose() * gradient + gradient.transpose() * f);
        m_bound = true;
        return *this;
      }

      /// @brief Returns the owned elastic kinematics.
      const KinematicState& getKinematicState() const
      {
        if (!m_bound)
        {
          Alert::Exception() << "Viscous material has no prepared constitutive state"
                             << Alert::Raise;
        }
        return m_law.getKinematicState();
      }

      /// @brief Returns elastic stored energy, excluding the dissipation potential.
      Real getStrainEnergyDensity() const
      {
        getKinematicState();
        return m_law.getStrainEnergyDensity();
      }

      /// @brief Returns nonnegative viscous power per reference volume.
      Real getDissipationRate() const
      {
        getKinematicState();
        return m_viscosity * m_strainRate.dot(m_strainRate);
      }

      /**
       * @brief Adds @f$\nu F\dot E@f$ to elastic first Piola stress.
       * @param stress Output total first Piola stress.
       */
      void getFirstPiolaKirchhoffStress(Math::SpatialMatrix<Real>& stress) const
      {
        const auto& f = getKinematicState().getDeformationGradient();
        m_law.getFirstPiolaKirchhoffStress(stress);
        stress = stress + m_viscosity * f * m_strainRate;
      }

      /**
       * @brief Adds the consistent viscous stress derivative to the elastic tangent.
       * @param tangent Output stress increment.
       * @param increment Deformation-gradient increment.
       */
      void getMaterialTangent(Math::SpatialMatrix<Real>& tangent,
        const Math::SpatialMatrix<Real>& increment) const
      {
        const auto& f = getKinematicState().getDeformationGradient();
        m_law.getMaterialTangent(tangent, increment);
        const Math::SpatialMatrix<Real> rateIncrement = Real(0.5) *
          (increment.transpose() * m_velocityGradient +
            m_velocityGradient.transpose() * increment +
            m_velocityCoefficient *
              (f.transpose() * increment + increment.transpose() * f));
        tangent = tangent + m_viscosity * (increment * m_strainRate + f * rateIncrement);
      }

    private:
      Law m_law;
      Real m_viscosity;
      Real m_velocityCoefficient;
      Math::SpatialMatrix<Real> m_velocityGradient, m_strainRate;
      bool m_bound;
  };
}

#endif
