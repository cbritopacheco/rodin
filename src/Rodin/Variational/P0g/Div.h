/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_P0G_DIV_H
#define RODIN_VARIATIONAL_P0G_DIV_H

/**
 * @file
 * @brief Divergence operator specialization for P0g (global constant) vector functions.
 *
 * For a globally constant vector field u in P0g:
 *   div(u) = 0.
 *
 * This holds on cells and on faces/boundary (trace choice irrelevant since it is zero).
 */

#include <type_traits>
#include <utility>

#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Variational/Div.h"
#include "Rodin/Variational/ShapeFunction.h"

#include "Rodin/Variational/P0g/ForwardDecls.h"

namespace Rodin::FormLanguage
{
  /// @brief Type traits for @c Div over a grid function: exposes the finite element
  /// space, the scalar type and the operand type.
  template <class Scalar, class Data, class Mesh>
  struct Traits<Variational::Div<Variational::GridFunction<Variational::P0g<Math::SpatialVector<Scalar>, Mesh>, Data>>>
  {
      /// @brief Finite element space type.
      using FESType = Variational::P0g<Math::SpatialVector<Scalar>, Mesh>;
      /// @brief Scalar value type.
      using ScalarType = Scalar;
      /// @brief Operand type.
      using OperandType = Variational::GridFunction<FESType, Data>;
  };

  /// @brief Type traits for @c Div over a shape function: exposes the finite element
  /// space, the shape function space, the scalar type and the operand type.
  template <class NestedDerived, class Scalar, class Mesh, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Div<
      Variational::ShapeFunction<NestedDerived, Variational::P0g<Math::SpatialVector<Scalar>, Mesh>, Space>>>
  {
      /// @brief Finite element space type.
      using FESType = Variational::P0g<Math::SpatialVector<Scalar>, Mesh>;
      /// @brief Shape function space the expression belongs to, trial or test.
      static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
      /// @brief Scalar value type.
      using ScalarType = Scalar;
      /// @brief Operand type.
      using OperandType = Variational::ShapeFunction<NestedDerived, FESType, Space>;
  };
}

namespace Rodin::Variational
{
  template <class Operand, class Derived>
  class DivBase;

  /**
   * @ingroup DivSpecializations
   * @brief Divergence of a P0g vector GridFunction (identically zero).
   */
  template <class Scalar, class Data, class Mesh>
  class Div<GridFunction<P0g<Math::SpatialVector<Scalar>, Mesh>, Data>> final
    : public DivBase<
        GridFunction<P0g<Math::SpatialVector<Scalar>, Mesh>, Data>,
        Div<GridFunction<P0g<Math::SpatialVector<Scalar>, Mesh>, Data>>>
  {
    public:
      /// @brief Finite element space type.
      using FESType = P0g<Math::SpatialVector<Scalar>, Mesh>;
      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      /// @brief Operand type.
      using OperandType = GridFunction<FESType, Data>;
      /// @brief Parent class type.
      using Parent = DivBase<OperandType, Div<OperandType>>;

      /// @brief Constructs the expression from its operand.
      explicit Div(const OperandType& u)
        : Parent(u)
      {}

      /// @brief Copy constructor.
      Div(const Div& other)
        : Parent(other)
      {}

      /// @brief Move constructor.
      Div(Div&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Interpolates div(u) at point p (always zero for P0g).
       */
      void interpolate(ScalarType& out, const Geometry::Point&) const
      {
        out = ScalarType(0);
      }

      constexpr
      /// @brief Returns the polynomial order used on a mesh entity.
      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return 0;
      }

      /// @brief Creates a polymorphic copy.
      Div* copy() const noexcept override
      {
        return new Div(*this);
      }
  };

  /**
   * @ingroup DivSpecializations
   * @brief Divergence of a P0g vector ShapeFunction (identically zero).
   */
  template <class NestedDerived, class Scalar, class Mesh, ShapeFunctionSpaceType Space>
  class Div<ShapeFunction<NestedDerived, P0g<Math::SpatialVector<Scalar>, Mesh>, Space>> final
    : public ShapeFunctionBase<
        Div<ShapeFunction<NestedDerived, P0g<Math::SpatialVector<Scalar>, Mesh>, Space>>,
        P0g<Math::SpatialVector<Scalar>, Mesh>,
        Space>
  {
    public:
      /// @brief Finite element space type.
      using FESType = P0g<Math::SpatialVector<Scalar>, Mesh>;
      /// @brief Shape function space the expression belongs to, trial or test.
      static constexpr ShapeFunctionSpaceType SpaceType = Space;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      /// @brief Operand type.
      using OperandType = ShapeFunction<NestedDerived, FESType, SpaceType>;

      /// @brief Parent class type.
      using Parent = ShapeFunctionBase<Div<OperandType>, FESType, SpaceType>;

      /// @brief Constructs the expression from its operand.
      explicit Div(const OperandType& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u),
          m_ip(nullptr),
          m_zero(ScalarType(0))
      {}

      /// @brief Copy constructor.
      Div(const Div& other)
        : Parent(other),
          m_u(other.m_u),
          m_ip(nullptr),
          m_zero(other.m_zero)
      {}

      /// @brief Move constructor.
      Div(Div&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_ip(std::exchange(other.m_ip, nullptr)),
          m_zero(std::exchange(other.m_zero, ScalarType(0)))
      {}

      constexpr
      /// @brief Gets the operand function.
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      constexpr
      /// @brief Gets the global DOF indices for a polytope.
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      constexpr
      /// @brief Gets the integration point the expression is evaluated at.
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      /// @brief Sets the integration point the expression is evaluated at.
      Div& setIntegrationPoint(const IntegrationPoint& ip)
      {
        // keep operand aligned
        m_u.get().setIntegrationPoint(ip);
        m_ip = &ip;
        m_zero = ScalarType(0);
        return *this;
      }

      /**
       * @brief Returns div(phi_local) (always zero).
       */
      constexpr
      ScalarType getBasis(size_t local) const
      {
        (void) local;
        assert(m_ip);
        return m_zero;
      }

      constexpr
      /// @brief Returns the polynomial order used on a mesh entity.
      Optional<size_t> getOrder(const Geometry::Polytope&) const noexcept
      {
        return 0;
      }

      Div* copy() const noexcept override
      {
        return new Div(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
      const IntegrationPoint* m_ip;
      ScalarType m_zero;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Div of a P0g GridFunction
   */
  template <class Scalar, class Data, class Mesh>
  Div(const GridFunction<P0g<Math::SpatialVector<Scalar>, Mesh>, Data>&)
    -> Div<GridFunction<P0g<Math::SpatialVector<Scalar>, Mesh>, Data>>;

  /// @brief Deduction guide for @c Div.
  template <class NestedDerived, class Scalar, class Mesh, ShapeFunctionSpaceType Space>
  Div(const ShapeFunction<NestedDerived, P0g<Math::SpatialVector<Scalar>, Mesh>, Space>&)
    -> Div<ShapeFunction<NestedDerived, P0g<Math::SpatialVector<Scalar>, Mesh>, Space>>;
}

#endif
