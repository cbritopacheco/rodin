/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_SWIFT_FITTINGTENSOR_H
#define RODIN_ADAPTATION_SWIFT_FITTINGTENSOR_H

#include <type_traits>

#include "Rodin/Variational/MatrixFunction.h"
#include "Rodin/Variational/VectorFunction.h"

#include "../DeformationMap.h"
#include "Parameters.h"

namespace Rodin::Adaptation::SWIFT
{
  /**
   * @brief Matrix coefficient of the SWIFT surface fitting metric.
   * Hessian of half the squared residual, with @f$D^2\phi@f$ omitted
   * and fixed gradient-scale normalization.
   */
  template <class GradDerived, class Displacement, class LocatorType>
  class FittingTensor final : public Variational::MatrixFunctionBase<Real,
                                FittingTensor<GradDerived, Displacement, LocatorType>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;
      /// @brief Range (evaluation value) type.
      using RangeType = Math::SpatialMatrix<ScalarType>;
      /// @brief Parent class type.
      using Parent = Variational::MatrixFunctionBase<ScalarType,
        FittingTensor<GradDerived, Displacement, LocatorType>>;
      /// @brief Level-set gradient function type.
      using GradType = Variational::VectorFunctionBase<Real, GradDerived>;

      /**
       * @brief Constructs the fitting metric coefficient.
       * @param grad Gradient of the target level set.
       * @param current Current displacement field.
       * @param locator Locator for deformed evaluation points.
       * @param parameters Model coefficients retained by reference.
       * @param normalization Fixed gradient-scale normalization.
       * @param dimension Spatial dimension.
       */
      FittingTensor(const GradType& grad, const Displacement& current,
        const LocatorType& locator, const Parameters& parameters, Real normalization,
        std::size_t dimension)
        : m_grad(grad.copy()),
          m_deformation(current, locator),
          m_parameters(parameters),
          m_normalization(normalization),
          m_dimension(dimension)
      {}

      /**
       * @brief Copy constructor.
       * @param other Coefficient to copy, cloning its target gradient.
       */
      FittingTensor(const FittingTensor& other)
        : Parent(other),
          m_grad(other.m_grad->copy()),
          m_deformation(other.m_deformation),
          m_parameters(other.m_parameters),
          m_normalization(other.m_normalization),
          m_dimension(other.m_dimension)
      {}

      /**
       * @brief Evaluates the coefficient at a point.
       * @param ip Integration point at which the expression is evaluated.
       * @returns Value of the expression at the supplied evaluation point.
       */
      RangeType getValue(const Variational::IntegrationPoint& ip) const
      {
        const auto& params = m_parameters.get();
        const auto d = static_cast<std::uint8_t>(m_dimension);
        RangeType m = RangeType::Identity(d, d);
        m.setZero();

        const auto g = m_grad->getValue(m_deformation.getMovedPoint(ip));
        const Real scale = params.model.fit * m_normalization;
        for (std::uint8_t r = 0; r < d; ++r)
        {
          for (std::uint8_t c = 0; c < d; ++c)
            m(r, c) = scale * g(r) * g(c);
        }
        return m;
      }

      /**
       * @brief Number of rows of the matrix value.
       * @returns The rows.
       */
      std::size_t getRows() const noexcept
      {
        return m_dimension;
      }

      /**
       * @brief Number of columns of the matrix value.
       * @returns The columns.
       */
      std::size_t getColumns() const noexcept
      {
        return m_dimension;
      }

      /**
       * @brief Polynomial order of the composed fitting tensor when known.
       * A globally constant vector is unchanged by deformation. Elementwise
       * orders cannot be forwarded through point location: a moved entity may
       * cross target cells, including cells with different constant gradients.
       * @returns Zero for a globally constant gradient, or an empty optional.
       * @param polytope Reference entity on which the tensor is integrated.
       */
      Optional<std::size_t> getOrder(
        [[maybe_unused]] const Geometry::Polytope& polytope) const noexcept
      {
        if constexpr (std::is_same_v<GradDerived,
                        Variational::VectorFunction<Math::Vector<Real>>>)
          return m_grad->getOrder(polytope);
        return std::nullopt;
      }

      /**
       * @brief Clones this coefficient.
       * @returns Newly allocated copy owned by the caller.
       */
      FittingTensor* copy() const noexcept override
      {
        return new FittingTensor(*this);
      }

    private:
      std::unique_ptr<GradType> m_grad;
      DeformationMap<Displacement, LocatorType> m_deformation;
      std::reference_wrapper<const Parameters> m_parameters;
      Real m_normalization;
      std::size_t m_dimension;
  };
  /**
   * @brief Deduction guide for FittingTensor.
   * @param grad Gradient of the observation field.
   * @param current Current displacement field.
   * @param locator Point locator used to find mesh entities.
   * @param parameters Parameters configuring the operation.
   * @param normalization Normalization factor.
   * @param dimension Spatial dimension.
   */

  template <class GradDerived, class Displacement, class LocatorType>
  FittingTensor(const Variational::VectorFunctionBase<Real, GradDerived>& grad,
    const Displacement& current, const LocatorType& locator, const Parameters& parameters,
    Real normalization,
    std::size_t dimension) -> FittingTensor<GradDerived, Displacement, LocatorType>;
}

#endif
