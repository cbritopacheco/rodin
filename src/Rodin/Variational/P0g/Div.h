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
  /**
   * @brief Type traits for @c Div over a grid function: exposes the finite element
   * space, the scalar type and the operand type.
   */
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

  /**
   * @brief Type traits for @c Div over a shape function: exposes the finite element
   * space, the shape function space, the scalar type and the operand type.
   */
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

      /**
       * @brief Constructs the expression from its operand.
       * @param u Operand expression.
       */
      explicit Div(const OperandType& u)
        : Parent(u)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      Div(const Div& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      Div(Div&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Interpolates div(u) at point p (always zero for P0g).
       * @param out Storage for the computed result.
       * @param point Evaluation point; the result is independent of this argument.
       */
      void interpolate(
        ScalarType& out, [[maybe_unused]] const Geometry::Point& point) const
      {
        out = ScalarType(0);
      }

      /**
       * @brief Returns the polynomial order used on a mesh entity.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       * @param polytope Mesh entity; the reported order is independent of this argument.
       */
      constexpr Optional<size_t> getOrder(
        [[maybe_unused]] const Geometry::Polytope& polytope) const noexcept
      {
        return 0;
      }

      /**
       * @brief Creates a polymorphic copy.
       * @returns Pointer to a newly allocated copy; the caller owns the returned object.
       */
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

      /**
       * @brief Constructs the expression from its operand.
       * @param u Operand expression.
       */
      explicit Div(const OperandType& u)
        : Parent(u.getFiniteElementSpace()),
          m_u(u),
          m_ip(nullptr),
          m_zero(ScalarType(0))
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      Div(const Div& other)
        : Parent(other),
          m_u(other.m_u),
          m_ip(nullptr),
          m_zero(other.m_zero)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      Div(Div&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u)),
          m_ip(std::exchange(other.m_ip, nullptr)),
          m_zero(std::exchange(other.m_zero, ScalarType(0)))
      {}

      /**
       * @brief Gets the operand function.
       * @returns The operand function.
       */
      constexpr
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      /**
       * @brief Returns the number of local basis functions for a polytope.
       * @param element Finite element used by the operation.
       * @returns Number of local basis functions on the selected entity.
       */
      constexpr
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      /**
       * @brief Gets the integration point the expression is evaluated at.
       * @returns The integration point the expression is evaluated at.
       */
      constexpr
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      /**
       * @brief Sets the integration point the expression is evaluated at.
       * @param ip Integration point at which the expression is evaluated.
       * @returns Reference to this object after the operation.
       */
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
       * @param local Index in the local numbering.
       * @returns Value of the selected local basis function at the evaluation point.
       */
      constexpr
      ScalarType getBasis(size_t local) const
      {
        (void) local;
        assert(m_ip);
        return m_zero;
      }

      /**
       * @brief Returns the polynomial order used on a mesh entity.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       * @param polytope Mesh entity; the reported order is independent of this argument.
       */
      constexpr Optional<size_t> getOrder(
        [[maybe_unused]] const Geometry::Polytope& polytope) const noexcept
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
   * @param u Operand expression.
   */
  template <class Scalar, class Data, class Mesh>
  Div(const GridFunction<P0g<Math::SpatialVector<Scalar>, Mesh>, Data>& u)
    -> Div<GridFunction<P0g<Math::SpatialVector<Scalar>, Mesh>, Data>>;

  /**
   * @brief Deduction guide for @c Div.
   * @param u Operand expression.
   */
  template <class NestedDerived, class Scalar, class Mesh, ShapeFunctionSpaceType Space>
  Div(
    const ShapeFunction<NestedDerived, P0g<Math::SpatialVector<Scalar>, Mesh>, Space>& u)
    -> Div<ShapeFunction<NestedDerived, P0g<Math::SpatialVector<Scalar>, Mesh>, Space>>;
}

#endif
