/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIROBSERVATIONCOEFFICIENT_H
#define RODIN_ADAPTATION_WNGIROBSERVATIONCOEFFICIENT_H

#include "Rodin/Variational/MatrixFunction.h"

#include "DeformationMap.h"
#include "WNGIRLoss.h"
#include "WNGIRParameters.h"
#include "WNGIRResidualState.h"

namespace Rodin::Adaptation::Detail
{
  /// @brief Matrix coefficient of the WNGIR surface observation metric.
  template <class PhiDerived, class GradDerived, class Displacement, class LocatorType>
  class WNGIRObservationCoefficient final
    : public Variational::MatrixFunctionBase<Real,
        WNGIRObservationCoefficient<PhiDerived, GradDerived, Displacement, LocatorType>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;
      /// @brief Range (evaluation value) type.
      using RangeType = Math::SpatialMatrix<ScalarType>;
      /// @brief Parent class type.
      using Parent = Variational::MatrixFunctionBase<ScalarType,
        WNGIRObservationCoefficient<PhiDerived, GradDerived, Displacement, LocatorType>>;
      /// @brief Level-set function type.
      using PhiType = Variational::RealFunctionBase<PhiDerived>;
      /// @brief Level-set gradient function type.
      using GradType = Variational::VectorFunctionBase<Real, GradDerived>;

      /// @brief Constructs the WNGIR observation coefficient.
      WNGIRObservationCoefficient(const PhiType& phi, const GradType& grad,
        const Displacement& current, const LocatorType& locator,
        const WNGIRParameters& parameters, Real sigma2, Real normalization,
        std::size_t dimension)
        : m_phi(phi.copy()),
          m_grad(grad.copy()),
          m_deformation(current, locator),
          m_parameters(parameters),
          m_loss(std::sqrt(sigma2)),
          m_normalization(normalization),
          m_dimension(dimension)
      {}

      /// @brief Copy constructor.
      WNGIRObservationCoefficient(const WNGIRObservationCoefficient& other)
        : Parent(other),
          m_phi(other.m_phi->copy()),
          m_grad(other.m_grad->copy()),
          m_deformation(other.m_deformation),
          m_parameters(other.m_parameters),
          m_loss(other.m_loss),
          m_normalization(other.m_normalization),
          m_dimension(other.m_dimension)
      {}

      /// @brief Evaluates the coefficient at a point.
      RangeType getValue(const Variational::IntegrationPoint& ip) const
      {
        constexpr Real epsG = Real(1e-12);
        const auto& params = m_parameters.get();
        const auto d = static_cast<std::uint8_t>(m_dimension);
        RangeType m = RangeType::Identity(d, d);

        const WNGIRResidualState state(*m_phi, *m_grad, m_deformation, ip, m_loss, true);
        const auto& g = state.getGradient();
        const Real g2 = g.dot(g);
        const Real tau = params.tauTan * g2;
        const Real normalizedG2 = m_normalization * g2;
        const Real tangential =
          params.tauTan * normalizedG2 / (normalizedG2 + epsG);
        const Real scale = params.kappaObs * m_normalization * state.getWeight();
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

      /// @brief Number of rows of the matrix value.
      std::size_t getRows() const noexcept
      {
        return m_dimension;
      }

      /// @brief Number of columns of the matrix value.
      std::size_t getColumns() const noexcept
      {
        return m_dimension;
      }

      /// @brief Reports no intrinsic polynomial order.
      Optional<std::size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return std::nullopt;
      }

      /// @brief Clones this object.
      WNGIRObservationCoefficient* copy() const noexcept override
      {
        return new WNGIRObservationCoefficient(*this);
      }

    private:
      std::unique_ptr<PhiType> m_phi;
      std::unique_ptr<GradType> m_grad;
      DeformationMap<Displacement, LocatorType> m_deformation;
      std::reference_wrapper<const WNGIRParameters> m_parameters;
      WNGIRLoss m_loss;
      Real m_normalization;
      std::size_t m_dimension;
  };

  template <class PhiDerived, class GradDerived, class Displacement, class LocatorType>
  WNGIRObservationCoefficient(const Variational::RealFunctionBase<PhiDerived>&,
    const Variational::VectorFunctionBase<Real, GradDerived>&, const Displacement&,
    const LocatorType&, const WNGIRParameters&, Real, Real, std::size_t)
    -> WNGIRObservationCoefficient<PhiDerived, GradDerived, Displacement, LocatorType>;
}

#endif
