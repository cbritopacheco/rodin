/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIROBSERVATIONCOEFFICIENT_H
#define RODIN_ADAPTATION_WNGIROBSERVATIONCOEFFICIENT_H

#include "Rodin/Variational/MatrixFunction.h"

#include "DeformationMap.h"
#include "WNGIRParameters.h"
#include "WNGIRResidualState.h"

namespace Rodin::Adaptation::Detail
{
  /// @brief Matrix coefficient of the WNGIR surface observation metric.
  template <class PhiDerived, class GradDerived, class Displacement,
    class LocatorType>
  class WNGIRObservationCoefficient final
    : public Variational::MatrixFunctionBase<Real,
        WNGIRObservationCoefficient<PhiDerived, GradDerived, Displacement,
          LocatorType>>
  {
    public:
      using ScalarType = Real;
      using RangeType = Math::SpatialMatrix<ScalarType>;
      using Parent = Variational::MatrixFunctionBase<ScalarType,
        WNGIRObservationCoefficient<PhiDerived, GradDerived, Displacement,
          LocatorType>>;
      using PhiType = Variational::RealFunctionBase<PhiDerived>;
      using GradType = Variational::VectorFunctionBase<Real, GradDerived>;

      WNGIRObservationCoefficient(const PhiType& phi, const GradType& grad,
        const Displacement& current, const LocatorType& locator,
        const WNGIRParameters& parameters, Real sigma2, std::size_t dimension)
        : m_phi(phi.copy()),
          m_grad(grad.copy()),
          m_deformation(current, locator),
          m_parameters(parameters),
          m_sigma2(sigma2),
          m_dimension(dimension)
      {}

      WNGIRObservationCoefficient(const WNGIRObservationCoefficient& other)
        : Parent(other),
          m_phi(other.m_phi->copy()),
          m_grad(other.m_grad->copy()),
          m_deformation(other.m_deformation),
          m_parameters(other.m_parameters),
          m_sigma2(other.m_sigma2),
          m_dimension(other.m_dimension)
      {}

      RangeType getValue(const Variational::IntegrationPoint& ip) const
      {
        constexpr Real epsG = Real(1e-12);
        const auto& params = m_parameters.get();
        const auto d = static_cast<std::uint8_t>(m_dimension);
        RangeType m = RangeType::Identity(d, d);

        if (params.observationMetric == WNGIRObservationMetric::Isotropic)
        {
          const WNGIRResidualState state(
            *m_phi, *m_grad, m_deformation, ip, m_sigma2, false);
          const Real g2 = state.getGradient().dot(state.getGradient());
          const Real a = params.gammaObs * (g2 + epsG +
            (params.residualStabilizedObservationMetric
              ? (state.getResidual() * state.getResidual()) / m_sigma2
              : Real(0)));
          m.setZero();
          for (std::uint8_t r = 0; r < d; ++r)
            m(r, r) = a;
          return m;
        }

        const bool hybrid =
          params.observationMetric == WNGIRObservationMetric::HybridRankOneIRLS;
        const bool weighted =
          params.observationMetric == WNGIRObservationMetric::RankOneIRLS || hybrid;
        const WNGIRResidualState state(
          *m_phi, *m_grad, m_deformation, ip, m_sigma2, weighted);
        const auto& g = state.getGradient();
        const Real g2 = g.dot(g);
        const Real tau = hybrid ? params.observationTangentialFloor * g2 : Real(0);
        const Real tangential = hybrid ? tau / (g2 + epsG) : Real(0);
        const Real scale = params.gammaObs * state.getWeight();
        for (std::uint8_t r = 0; r < d; ++r)
        {
          for (std::uint8_t c = 0; c < d; ++c)
          {
            m(r, c) = scale * g(r) * g(c) * (Real(1) - tangential);
            if (r == c)
              m(r, c) += scale * tau;
          }
        }
        return m;
      }

      std::size_t getRows() const noexcept { return m_dimension; }

      std::size_t getColumns() const noexcept { return m_dimension; }

      Optional<std::size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return std::nullopt;
      }

      WNGIRObservationCoefficient* copy() const noexcept override
      {
        return new WNGIRObservationCoefficient(*this);
      }

    private:
      std::unique_ptr<PhiType> m_phi;
      std::unique_ptr<GradType> m_grad;
      DeformationMap<Displacement, LocatorType> m_deformation;
      std::reference_wrapper<const WNGIRParameters> m_parameters;
      Real m_sigma2;
      std::size_t m_dimension;
  };

  template <class PhiDerived, class GradDerived, class Displacement,
    class LocatorType>
  WNGIRObservationCoefficient(const Variational::RealFunctionBase<PhiDerived>&,
    const Variational::VectorFunctionBase<Real, GradDerived>&,
    const Displacement&, const LocatorType&, const WNGIRParameters&, Real,
    std::size_t)
    -> WNGIRObservationCoefficient<PhiDerived, GradDerived, Displacement,
         LocatorType>;
}

#endif
