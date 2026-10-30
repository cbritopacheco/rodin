/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_SOLID_CONSTITUTIVE_DEMIRAYHOLZAPFELGASSEROGDEN_H
#define RODIN_SOLID_CONSTITUTIVE_DEMIRAYHOLZAPFELGASSEROGDEN_H

#include <cmath>
#include "HyperElasticLaw.h"
#include "Rodin/Solid/Local/FiberKinematics.h"

namespace Rodin::Solid
{
  /**
   * @brief Three-dimensional Demiray matrix with one HGO-type reinforcement family.
   *
   * @f[
   * W=\frac{a}{2b}(e^{b(\bar I_1-3)}-1)
   *   +\frac{a_f}{2b_f}(e^{b_f(\bar I_4-1)^2}-1)+\mu_v(\log J)^2,
   * \qquad \bar I_i=J^{-2/3}I_i.
   * @f]
   * The name identifies a hybrid of the Demiray matrix and
   * Holzapfel–Gasser–Ogden-type directional reinforcement. The stated
   * energy defines this single-family variant and its volumetric penalty.
   * The exponential terms use their continuous limits when an exponent is zero.
   * The volumetric coefficient gives initial bulk modulus @f$2\mu_v@f$.
   * Additive constants set the reference energy to zero without changing stress.
   * Stress and tangent are analytic derivatives of this same energy.
   */
  class DemirayHolzapfelGasserOgden final
    : public HyperElasticLaw<DemirayHolzapfelGasserOgden>
  {
    public:
      /// @brief Material coefficients in consistent stress units.
      struct Parameters
      {
        Real matrixScale = 20; ///< Isotropic exponential scale @f$a>0@f$.
        Real matrixExponent = 5; ///< Dimensionless exponent @f$b\geq0@f$.
        Real fiberScale = 20; ///< Fiber exponential scale @f$a_f\geq0@f$.
        Real fiberExponent = 10; ///< Dimensionless exponent @f$b_f\geq0@f$.
        Real volumePenalty = 500; ///< Coefficient @f$\mu_v>0@f$ of squared log volume.
      };

      /**
       * @brief Constructs the material from positive, finite coefficients.
       * @param parameters Elastic coefficients; scales have stress units.
       */
      explicit DemirayHolzapfelGasserOgden(const Parameters& parameters)
        : m_parameters(parameters),
          m_logJ(0),
          m_scale(1),
          m_first(3),
          m_fourth(1),
          m_firstDerivative(0),
          m_fiberDerivative(0),
          m_fiberCurvature(0),
          m_volumeStress(0)
      {
        const auto& p = m_parameters;
        if (!std::isfinite(p.matrixScale) || p.matrixScale <= 0 ||
          !std::isfinite(p.matrixExponent) || p.matrixExponent < 0 ||
          !std::isfinite(p.fiberScale) || p.fiberScale < 0 ||
          !std::isfinite(p.fiberExponent) || p.fiberExponent < 0 ||
          !std::isfinite(p.volumePenalty) || p.volumePenalty <= 0)
        {
          Alert::Exception() << "Invalid Demiray–Holzapfel–Gasser–Ogden coefficients"
                             << Alert::Raise;
        }
      }

      /// @brief Returns the elastic coefficients.
      const Parameters& getParameters() const noexcept
      {
        return m_parameters;
      }

      /**
       * @brief Captures kinematics and the reference fiber, preparing invariants.
       * @param point Constitutive point with three-dimensional kinematics.
       * @returns This material.
       */
      DemirayHolzapfelGasserOgden& setConstitutivePoint(const ConstitutivePoint& point)
      {
        const auto& state = point.getKinematicState();
        if (state.getDimension() != 3 || !std::isfinite(state.getJacobian()) ||
          state.getJacobian() <= 0)
        {
          m_state.reset();
          Alert::Exception() << "Demiray–Holzapfel–Gasser–Ogden requires a positive "
                                "three-dimensional deformation"
                             << Alert::Raise;
        }
        const auto fiber = FiberKinematics::direction(point);
        if (m_state && hasSameKinematics(*m_state, state) &&
          hasSameDirection(m_fiber, fiber))
          return *this;
        m_state.reset();
        m_fiber = fiber;
        const auto& f = state.getDeformationGradient();
        m_logJ = state.getLogJacobian();
        m_scale = std::exp(-Real(2) * m_logJ / 3);
        m_first = m_scale * f.dot(f);
        const auto stretched = f * m_fiber;
        m_fourth = m_scale * stretched.dot(stretched);
        const Real offset = m_fourth - 1;
        m_firstDerivative = Real(0.5) * m_parameters.matrixScale *
          std::exp(m_parameters.matrixExponent * (m_first - 3));
        const Real fiberExponential = m_parameters.fiberScale *
          std::exp(m_parameters.fiberExponent * offset * offset);
        m_fiberDerivative = fiberExponential * offset;
        m_fiberCurvature =
          fiberExponential * (1 + 2 * m_parameters.fiberExponent * offset * offset);
        m_volumeStress = 2 * m_parameters.volumePenalty * m_logJ -
          Real(2) / 3 * (m_firstDerivative * m_first + m_fiberDerivative * m_fourth);
        m_state.emplace(state);
        return *this;
      }

      /// @brief Returns the owned, prepared kinematic state.
      const KinematicState& getKinematicState() const
      {
        this->require(m_state.has_value(), __func__);
        return *m_state;
      }

      /// @brief Returns stored elastic energy, zero in the reference state.
      Real getStrainEnergyDensity() const
      {
        getKinematicState();
        const auto& p = m_parameters;
        const Real firstOffset = m_first - 3;
        const Real fiberOffset = (m_fourth - 1) * (m_fourth - 1);
        const Real matrixEnergy = p.matrixExponent == 0
          ? firstOffset
          : std::expm1(p.matrixExponent * firstOffset) / p.matrixExponent;
        const Real fiberEnergy = p.fiberExponent == 0
          ? fiberOffset
          : std::expm1(p.fiberExponent * fiberOffset) / p.fiberExponent;
        return Real(0.5) * (p.matrixScale * matrixEnergy + p.fiberScale * fiberEnergy) +
          p.volumePenalty * m_logJ * m_logJ;
      }

      /**
       * @brief Computes @f$P=\partial_F W@f$.
       * @param stress Output first Piola stress.
       */
      void getFirstPiolaKirchhoffStress(Math::SpatialMatrix<Real>& stress) const
      {
        const auto& state = getKinematicState();
        const auto& f = state.getDeformationGradient();
        const auto fiberTensor = FiberKinematics::dyad(m_fiber);
        stress =
          (2 * m_scale) * (m_firstDerivative * f + m_fiberDerivative * f * fiberTensor) +
          m_volumeStress * state.getDeformationGradientInverseTranspose();
      }

      /**
       * @brief Computes the analytic directional derivative @f$DP[\delta F]@f$.
       * @param tangent Output stress increment.
       * @param increment Deformation-gradient increment.
       */
      void getMaterialTangent(Math::SpatialMatrix<Real>& tangent,
        const Math::SpatialMatrix<Real>& increment) const
      {
        const auto& state = getKinematicState();
        const auto& f = state.getDeformationGradient();
        const auto& inverse = state.getDeformationGradientInverseTranspose();
        const auto fiberTensor = FiberKinematics::dyad(m_fiber);
        const Real volumeIncrement = inverse.dot(increment);
        const Real scaleIncrement = -Real(2) / 3 * m_scale * volumeIncrement;
        const Real firstIncrement =
          2 * m_scale * f.dot(increment) - Real(2) / 3 * m_first * volumeIncrement;
        const Real fourthIncrement =
          2 * m_scale * (f * m_fiber).dot(increment * m_fiber) -
          Real(2) / 3 * m_fourth * volumeIncrement;
        const Real firstDerivativeIncrement =
          m_parameters.matrixExponent * m_firstDerivative * firstIncrement;
        const Real fiberDerivativeIncrement = m_fiberCurvature * fourthIncrement;
        const Real volumeStressIncrement =
          2 * m_parameters.volumePenalty * volumeIncrement -
          Real(2) / 3 *
            (firstDerivativeIncrement * m_first + m_firstDerivative * firstIncrement +
              fiberDerivativeIncrement * m_fourth + m_fiberDerivative * fourthIncrement);
        tangent = 2 * scaleIncrement *
            (m_firstDerivative * f + m_fiberDerivative * f * fiberTensor) +
          2 * m_scale *
            (firstDerivativeIncrement * f + m_firstDerivative * increment +
              fiberDerivativeIncrement * f * fiberTensor +
              m_fiberDerivative * increment * fiberTensor) +
          volumeStressIncrement * inverse -
          m_volumeStress * inverse * increment.transpose() * inverse;
      }

    private:
      Parameters m_parameters;
      Optional<KinematicState> m_state;
      Math::SpatialVector<Real> m_fiber;
      Real m_logJ, m_scale, m_first, m_fourth;
      Real m_firstDerivative, m_fiberDerivative, m_fiberCurvature;
      Real m_volumeStress;
  };
}

#endif
