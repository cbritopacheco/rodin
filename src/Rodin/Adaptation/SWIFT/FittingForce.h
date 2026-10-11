/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#ifndef RODIN_ADAPTATION_SWIFT_FITTINGFORCE_H
#define RODIN_ADAPTATION_SWIFT_FITTINGFORCE_H

#include "Rodin/Variational/VectorFunction.h"

#include "../DeformationMap.h"
#include "Loss.h"
#include "Residual.h"

namespace Rodin::Adaptation::SWIFT
{
  /// @brief Negative first-variation coefficient of the robust interface energy.
  template <class PhiDerived, class GradDerived, class Displacement, class LocatorType>
  class FittingForce final
    : public Variational::VectorFunctionBase<Real,
        FittingForce<PhiDerived, GradDerived, Displacement, LocatorType>>
  {
    public:
      /// @brief Scalar value type.
      using ScalarType = Real;
      /// @brief Range (evaluation value) type.
      using RangeType = Math::SpatialVector<ScalarType>;
      /// @brief Parent class type.
      using Parent = Variational::VectorFunctionBase<ScalarType,
        FittingForce<PhiDerived, GradDerived, Displacement, LocatorType>>;
      /// @brief Level-set function type.
      using PhiType = Variational::RealFunctionBase<PhiDerived>;
      /// @brief Level-set gradient function type.
      using GradType = Variational::VectorFunctionBase<Real, GradDerived>;

      /**
       * @brief Constructs the robust fitting force coefficient.
       * @param phi Target level set.
       * @param grad Gradient of the target level set.
       * @param current Current displacement field.
       * @param locator Locator for deformed evaluation points.
       * @param loss Robust residual loss.
       * @param normalization Fixed gradient-scale normalization.
       * @param dimension Spatial dimension.
       */
      FittingForce(const PhiType& phi, const GradType& grad, const Displacement& current,
        const LocatorType& locator, const Loss& loss, Real normalization,
        std::size_t dimension)
        : m_phi(phi.copy()),
          m_grad(grad.copy()),
          m_deformation(current, locator),
          m_loss(loss),
          m_normalization(normalization),
          m_dimension(dimension)
      {}

      /**
       * @brief Copy constructor.
       * @param other Coefficient to copy, cloning its target expressions.
       */
      FittingForce(const FittingForce& other)
        : Parent(other),
          m_phi(other.m_phi->copy()),
          m_grad(other.m_grad->copy()),
          m_deformation(other.m_deformation),
          m_loss(other.m_loss),
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
        const ResidualState state(*m_phi, *m_grad, m_deformation, ip, m_loss);
        return (-m_normalization * state.getWeight() * state.getResidual()) *
          state.getGradient();
      }

      /**
       * @brief Dimension of the vector value.
       * @returns The dimension.
       */
      std::size_t getDimension() const noexcept
      {
        return m_dimension;
      }

      /**
       * @brief Reports no intrinsic polynomial order.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       * @param polytope Mesh entity; the reported order is independent of this argument.
       */
      Optional<std::size_t> getOrder(
        [[maybe_unused]] const Geometry::Polytope& polytope) const noexcept
      {
        return std::nullopt;
      }

      /**
       * @brief Clones this coefficient.
       * @returns Newly allocated copy owned by the caller.
       */
      FittingForce* copy() const noexcept override
      {
        return new FittingForce(*this);
      }

    private:
      std::unique_ptr<PhiType> m_phi;
      std::unique_ptr<GradType> m_grad;
      DeformationMap<Displacement, LocatorType> m_deformation;
      Loss m_loss;
      Real m_normalization;
      std::size_t m_dimension;
  };
  /**
   * @brief Deduction guide for FittingForce.
   * @param phi Observation field.
   * @param grad Gradient of the observation field.
   * @param current Current displacement field.
   * @param locator Point locator used to find mesh entities.
   * @param loss Robust residual loss.
   * @param normalization Normalization factor.
   * @param dimension Spatial dimension.
   */

  template <class PhiDerived, class GradDerived, class Displacement, class LocatorType>
  FittingForce(const Variational::RealFunctionBase<PhiDerived>& phi,
    const Variational::VectorFunctionBase<Real, GradDerived>& grad,
    const Displacement& current, const LocatorType& locator, const Loss& loss,
    Real normalization, std::size_t dimension)
    -> FittingForce<PhiDerived, GradDerived, Displacement, LocatorType>;
}

#endif
