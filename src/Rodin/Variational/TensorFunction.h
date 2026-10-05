/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file TensorFunction.h
 * @brief Constant or callable spatial-tensor coefficients in variational forms.
 */
#ifndef RODIN_VARIATIONAL_TENSORFUNCTION_H
#define RODIN_VARIATIONAL_TENSORFUNCTION_H
#include "Rodin/Math/SpatialTensor.h"
#include "Function.h"

namespace Rodin::Variational
{
  /**
   * @brief Tensor coefficient with value semantics.
   *
   * A constant @ref Rodin::Math::SpatialTensor or a callable returning one is
   * evaluated through the same @ref FunctionBase interface as other ranges.
   */
  template <class Value>
  class TensorFunction final : public FunctionBase<TensorFunction<Value>>
  {
    public:
      /// @brief CRTP or finite element base class.
      using Parent = FunctionBase<TensorFunction>;
      /// @brief Owns a constant or callable tensor coefficient.
      /// @param value Value to store or assign.
      explicit TensorFunction(Value value)
        : m_value(std::move(value))
      {}
      /// @brief Owns a constant or callable tensor coefficient.
      /// @param other Object to copy from.
      TensorFunction(const TensorFunction& other)
        : Parent(other),
          m_value(other.m_value)
      {}
      /// @brief Owns a constant or callable tensor coefficient.
      /// @param other Object to move from.
      TensorFunction(TensorFunction&& other)
        : Parent(std::move(other)),
          m_value(std::move(other.m_value))
      {}
      /// @brief Evaluates the expression at the supplied physical or integration point.
      /// @param point Point at which the operation is evaluated.
      /// @returns Value of the expression at the supplied evaluation point.
      auto getValue(const Geometry::Point& point) const
      {
        if constexpr (std::is_invocable_v<Value, const Geometry::Point&>)
          return m_value(point);
        else
          return m_value;
      }
      /// @brief Evaluates the expression at the supplied physical or integration point.
      /// @param point Point at which the operation is evaluated.
      /// @returns Value of the expression at the supplied evaluation point.
      auto getValue(const IntegrationPoint& point) const
      {
        if constexpr (std::is_invocable_v<Value, const IntegrationPoint&>)
          return m_value(point);
        else
          return getValue(point.getPoint());
      }
      /// @brief Returns the polynomial order when it is known.
      /// @returns Polynomial order on the entity, or an empty optional when no order is available.
      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        if constexpr (FormLanguage::IsTensorRange<Value>::Value)
          return 0;
        else
          return std::nullopt;
      }
      TensorFunction* copy() const noexcept override
      {
        return new TensorFunction(*this);
      }

    private:
      Value m_value;
  };
  /// @brief Deduces the matrix space or coefficient type from constructor arguments.
  template <class Value>
  TensorFunction(Value) -> TensorFunction<Value>;
}
#endif
