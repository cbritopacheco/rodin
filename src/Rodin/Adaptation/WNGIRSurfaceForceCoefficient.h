/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIRSURFACEFORCECOEFFICIENT_H
#define RODIN_ADAPTATION_WNGIRSURFACEFORCECOEFFICIENT_H

#include "Rodin/Variational/VectorFunction.h"

#include "DeformationMap.h"
#include "WNGIRResidualState.h"
#include "WNGIRValidationWeights.h"

namespace Rodin::Adaptation::Detail
{
  /// @brief Negative first-variation coefficient of the Welsch interface energy.
  template <class PhiDerived, class GradDerived, class Displacement, class LocatorType>
  class WNGIRSurfaceForceCoefficient final
    : public Variational::VectorFunctionBase<Real,
        WNGIRSurfaceForceCoefficient<PhiDerived, GradDerived, Displacement, LocatorType>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;
      /// @brief Range (evaluation value) type.
      using RangeType = Math::SpatialVector<ScalarType>;
      /// @brief Parent class type.
      using Parent = Variational::VectorFunctionBase<ScalarType,
        WNGIRSurfaceForceCoefficient<PhiDerived, GradDerived, Displacement, LocatorType>>;
      /// @brief Level-set function type.
      using PhiType = Variational::RealFunctionBase<PhiDerived>;
      /// @brief Level-set gradient function type.
      using GradType = Variational::VectorFunctionBase<Real, GradDerived>;

      /// @brief Constructs the w n g i r surface force coefficient.
      WNGIRSurfaceForceCoefficient(const PhiType& phi, const GradType& grad,
        const Displacement& current, const LocatorType& locator, const WNGIRParameters&,
        Real sigma2, std::size_t dimension)
        : m_phi(phi.copy()),
          m_grad(grad.copy()),
          m_deformation(current, locator),
          m_sigma2(sigma2),
          m_dimension(dimension)
      {}

      /// @brief Copy constructor.
      WNGIRSurfaceForceCoefficient(const WNGIRSurfaceForceCoefficient& other)
        : Parent(other),
          m_phi(other.m_phi->copy()),
          m_grad(other.m_grad->copy()),
          m_deformation(other.m_deformation),
          m_sigma2(other.m_sigma2),
          m_dimension(other.m_dimension)
      {}

      /// @brief Evaluates the coefficient at a point.
      RangeType getValue(const Variational::IntegrationPoint& ip) const
      {
        const WNGIRResidualState state(
          *m_phi, *m_grad, m_deformation, ip, m_sigma2, true);
        const Real quadratureCorrection = WNGIRValidationWeights::getCorrection(ip);
        return (-quadratureCorrection * state.getWeight() * state.getResidual()) *
          state.getGradient();
      }

      /// @brief Dimension of the vector value.
      std::size_t getDimension() const noexcept
      {
        return m_dimension;
      }

      /// @brief Reports no intrinsic polynomial order.
      Optional<std::size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return std::nullopt;
      }

      WNGIRSurfaceForceCoefficient* copy() const noexcept override
      {
        return new WNGIRSurfaceForceCoefficient(*this);
      }

    private:
      std::unique_ptr<PhiType> m_phi;
      std::unique_ptr<GradType> m_grad;
      DeformationMap<Displacement, LocatorType> m_deformation;
      Real m_sigma2;
      std::size_t m_dimension;
  };

  template <class PhiDerived, class GradDerived, class Displacement, class LocatorType>
  WNGIRSurfaceForceCoefficient(const Variational::RealFunctionBase<PhiDerived>&,
    const Variational::VectorFunctionBase<Real, GradDerived>&, const Displacement&,
    const LocatorType&, const WNGIRParameters&, Real, std::size_t)
    -> WNGIRSurfaceForceCoefficient<PhiDerived, GradDerived, Displacement, LocatorType>;
}

#endif
