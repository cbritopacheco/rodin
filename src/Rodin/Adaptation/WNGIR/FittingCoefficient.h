/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIR_FITTINGCOEFFICIENT_H
#define RODIN_ADAPTATION_WNGIR_FITTINGCOEFFICIENT_H

#include "Rodin/Variational/MatrixFunction.h"

#include "../DeformationMap.h"
#include "Parameters.h"

namespace Rodin::Adaptation
{
  /// @brief Matrix coefficient of the WNGIR surface observation metric.
  /// Hessian of half the squared residual, with D2 phi omitted and fixed normalization.
  template <class GradDerived, class Displacement, class LocatorType>
  class WNGIRObservationCoefficient final
    : public Variational::MatrixFunctionBase<Real,
        WNGIRObservationCoefficient<GradDerived, Displacement, LocatorType>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;
      /// @brief Range (evaluation value) type.
      using RangeType = Math::SpatialMatrix<ScalarType>;
      /// @brief Parent class type.
      using Parent = Variational::MatrixFunctionBase<ScalarType,
        WNGIRObservationCoefficient<GradDerived, Displacement, LocatorType>>;
      /// @brief Level-set gradient function type.
      using GradType = Variational::VectorFunctionBase<Real, GradDerived>;

      /// @brief Constructs the WNGIR observation coefficient.
      WNGIRObservationCoefficient(const GradType& grad, const Displacement& current,
        const LocatorType& locator, const WNGIRParameters& parameters, Real normalization,
        std::size_t dimension)
        : m_grad(grad.copy()),
          m_deformation(current, locator),
          m_parameters(parameters),
          m_normalization(normalization),
          m_dimension(dimension)
      {}

      /// @brief Copy constructor.
      WNGIRObservationCoefficient(const WNGIRObservationCoefficient& other)
        : Parent(other),
          m_grad(other.m_grad->copy()),
          m_deformation(other.m_deformation),
          m_parameters(other.m_parameters),
          m_normalization(other.m_normalization),
          m_dimension(other.m_dimension)
      {}

      /// @brief Evaluates the coefficient at a point.
      RangeType getValue(const Variational::IntegrationPoint& ip) const
      {
        const auto& params = m_parameters.get();
        const auto d = static_cast<std::uint8_t>(m_dimension);
        RangeType m = RangeType::Identity(d, d);
        m.setZero();

        const auto g = m_grad->getValue(m_deformation.getMovedPoint(ip));
        const Real scale = params.kappaF * m_normalization;
        for (std::uint8_t r = 0; r < d; ++r)
        {
          for (std::uint8_t c = 0; c < d; ++c)
            m(r, c) = scale * g(r) * g(c);
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
      std::unique_ptr<GradType> m_grad;
      DeformationMap<Displacement, LocatorType> m_deformation;
      std::reference_wrapper<const WNGIRParameters> m_parameters;
      Real m_normalization;
      std::size_t m_dimension;
  };

  template <class GradDerived, class Displacement, class LocatorType>
  WNGIRObservationCoefficient(const Variational::VectorFunctionBase<Real, GradDerived>&,
    const Displacement&, const LocatorType&, const WNGIRParameters&, Real,
    std::size_t) -> WNGIRObservationCoefficient<GradDerived, Displacement, LocatorType>;
}

#endif
