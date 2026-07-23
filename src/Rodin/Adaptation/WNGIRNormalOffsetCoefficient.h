/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_WNGIRNORMALOFFSETCOEFFICIENT_H
#define RODIN_ADAPTATION_WNGIRNORMALOFFSETCOEFFICIENT_H

#include "Rodin/Variational/MatrixFunction.h"

#include "WNGIRNormalOffset.h"

namespace Rodin::Adaptation::Detail
{
  /// @brief Rank-one metric coefficient of the normal-offset initial guess.
  template <class PhiDerived, class GradDerived>
  class WNGIRNormalOffsetCoefficient final
    : public Variational::MatrixFunctionBase<Real,
        WNGIRNormalOffsetCoefficient<PhiDerived, GradDerived>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;
      /// @brief Range (evaluation value) type.
      using RangeType = Math::SpatialMatrix<ScalarType>;
      /// @brief Parent class type.
      using Parent = Variational::MatrixFunctionBase<ScalarType,
        WNGIRNormalOffsetCoefficient<PhiDerived, GradDerived>>;
      /// @brief Level-set function type.
      using PhiType = Variational::RealFunctionBase<PhiDerived>;
      /// @brief Level-set gradient function type.
      using GradType = Variational::VectorFunctionBase<Real, GradDerived>;

      /// @brief Constructs the w n g i r normal offset coefficient.
      WNGIRNormalOffsetCoefficient(const PhiType& phi, const GradType& grad,
        const WNGIRParameters& parameters, Real sigma2, std::size_t dimension)
        : m_phi(phi.copy()),
          m_grad(grad.copy()),
          m_parameters(parameters),
          m_sigma2(sigma2),
          m_dimension(dimension)
      {}

      /// @brief Copy constructor.
      WNGIRNormalOffsetCoefficient(const WNGIRNormalOffsetCoefficient& other)
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
        const auto d = static_cast<std::uint8_t>(m_dimension);
        RangeType m = RangeType::Identity(d, d);
        m.setZero();
        const WNGIRNormalOffset offset(*m_phi, *m_grad, p, m_parameters.get(), m_sigma2);
        if (offset.isDegenerate())
          return m;
        const Real scale = m_parameters.get().initialGuessGamma * offset.getAttenuation();
        for (std::uint8_t r = 0; r < d; ++r)
          for (std::uint8_t c = 0; c < d; ++c)
            m(r, c) = scale * offset.getNormal()(r) * offset.getNormal()(c);
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
      WNGIRNormalOffsetCoefficient* copy() const noexcept override
      {
        return new WNGIRNormalOffsetCoefficient(*this);
      }

    private:
      std::unique_ptr<PhiType> m_phi;
      std::unique_ptr<GradType> m_grad;
      std::reference_wrapper<const WNGIRParameters> m_parameters;
      Real m_sigma2;
      std::size_t m_dimension;
  };

  template <class PhiDerived, class GradDerived>
  WNGIRNormalOffsetCoefficient(const Variational::RealFunctionBase<PhiDerived>&,
    const Variational::VectorFunctionBase<Real, GradDerived>&, const WNGIRParameters&,
    Real, std::size_t) -> WNGIRNormalOffsetCoefficient<PhiDerived, GradDerived>;
}

#endif
