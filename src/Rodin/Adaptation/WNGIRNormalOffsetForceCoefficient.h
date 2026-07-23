/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIRNORMALOFFSETFORCECOEFFICIENT_H
#define RODIN_ADAPTATION_WNGIRNORMALOFFSETFORCECOEFFICIENT_H

#include "Rodin/Variational/VectorFunction.h"

#include "WNGIRNormalOffset.h"

namespace Rodin::Adaptation::Detail
{
  /// @brief Load coefficient of the normal-offset initial guess.
  template <class PhiDerived, class GradDerived>
  class WNGIRNormalOffsetForceCoefficient final
    : public Variational::VectorFunctionBase<Real,
        WNGIRNormalOffsetForceCoefficient<PhiDerived, GradDerived>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;
      /// @brief Range (evaluation value) type.
      using RangeType = Math::SpatialVector<ScalarType>;
      /// @brief Parent class type.
      using Parent = Variational::VectorFunctionBase<ScalarType,
        WNGIRNormalOffsetForceCoefficient<PhiDerived, GradDerived>>;
      /// @brief Level-set function type.
      using PhiType = Variational::RealFunctionBase<PhiDerived>;
      /// @brief Level-set gradient function type.
      using GradType = Variational::VectorFunctionBase<Real, GradDerived>;

      /// @brief Constructs the w n g i r normal offset force coefficient.
      WNGIRNormalOffsetForceCoefficient(const PhiType& phi, const GradType& grad,
        const WNGIRParameters& parameters, Real sigma2, std::size_t dimension)
        : m_phi(phi.copy()),
          m_grad(grad.copy()),
          m_parameters(parameters),
          m_sigma2(sigma2),
          m_dimension(dimension)
      {}

      /// @brief Copy constructor.
      WNGIRNormalOffsetForceCoefficient(const WNGIRNormalOffsetForceCoefficient& other)
        : Parent(other),
          m_phi(other.m_phi->copy()),
          m_grad(other.m_grad->copy()),
          m_parameters(other.m_parameters),
          m_sigma2(other.m_sigma2),
          m_dimension(other.m_dimension)
      {}

      /// @brief Evaluates the coefficient at a point.
      RangeType getValue(const Geometry::Point& p) const
      {
        const WNGIRNormalOffset offset(*m_phi, *m_grad, p, m_parameters.get(), m_sigma2);
        if (offset.isDegenerate())
          return RangeType::Zero(m_dimension);
        return (m_parameters.get().initialGuessGamma * offset.getAttenuation() *
                 offset.getOffset()) *
          offset.getNormal();
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

      WNGIRNormalOffsetForceCoefficient* copy() const noexcept override
      {
        return new WNGIRNormalOffsetForceCoefficient(*this);
      }

    private:
      std::unique_ptr<PhiType> m_phi;
      std::unique_ptr<GradType> m_grad;
      std::reference_wrapper<const WNGIRParameters> m_parameters;
      Real m_sigma2;
      std::size_t m_dimension;
  };

  template <class PhiDerived, class GradDerived>
  WNGIRNormalOffsetForceCoefficient(const Variational::RealFunctionBase<PhiDerived>&,
    const Variational::VectorFunctionBase<Real, GradDerived>&, const WNGIRParameters&,
    Real, std::size_t) -> WNGIRNormalOffsetForceCoefficient<PhiDerived, GradDerived>;
}

#endif
