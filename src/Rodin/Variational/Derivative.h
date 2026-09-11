/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Derivative.h
 * @brief Directional derivative operator for scalar functions.
 *
 * This file defines the Derivative class, which computes directional derivatives
 * of scalar functions. This generalizes the gradient concept to derivatives
 * in specific directions.
 *
 * ## Mathematical Foundation
 * The directional derivative of a function @f$ u @f$ in direction @f$ \mathbf{v} @f$ is:
 * @f[
 *   D_{\mathbf{v}} u = \nabla u \cdot \mathbf{v} = \lim_{h \to 0} \frac{u(x + h\mathbf{v}) - u(x)}{h}
 * @f]
 *
 * ## Coordinate Derivatives
 * Special cases are partial derivatives:
 * - @f$ \frac{\partial u}{\partial x_i} = D_{\mathbf{e}_i} u @f$
 * where @f$ \mathbf{e}_i @f$ is the @f$ i @f$-th coordinate direction.
 *
 * ## Applications
 * - Normal derivatives: @f$ \frac{\partial u}{\partial n} = \nabla u \cdot \mathbf{n} @f$
 * - Material derivatives in fluid dynamics
 * - Characteristic methods for PDEs
 * - Shape sensitivity analysis
 *
 * ## Usage Example
 * ```cpp
 * // Normal derivative on boundary
 * auto n = BoundaryNormal();
 * auto normal_deriv = Derivative(u, n);  // \partialu/\partialn = \nabla u\cdotn
 * ```
 */
#ifndef RODIN_VARIATIONAL_DERIVATIVE_H
#define RODIN_VARIATIONAL_DERIVATIVE_H

#include <cassert>
#include <cstdlib>

#include "ForwardDecls.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "ShapeFunction.h"

namespace Rodin::FormLanguage
{
  /// @brief Type traits for @c Derivative over a shape function: exposes the shape
  /// function space, the finite element space, the scalar type and the operand type.
  template <class NestedDerived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<Variational::Derivative<Variational::ShapeFunction<NestedDerived, FES, Space>>>
  {
    /// @brief Shape function space the expression belongs to, trial or test.
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;
    /// @brief Finite element space type.
    using FESType = FES;
    /// @brief Scalar value type.
    using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;
    /// @brief Operand type.
    using OperandType = Variational::ShapeFunction<NestedDerived, FESType, Space>;
  };
}

/**
 * @defgroup DerivativeSpecializations Derivative Template Specializations
 * @brief Template specializations of the Derivative class.
 * @see @ref Derivative
 *
 * | Specialization | Description |
 * |----------------|-------------|
 * | @ref Derivative "Derivative<H1<K, Scalar, Mesh>, ShapeFunction<NestedDerived, H1<K, Scalar, Mesh>, Space>>" | Directional derivative of an H1 shape function. |
 * | @ref Derivative "Derivative<P1<Range, Mesh>, GridFunction<P1<Range, Mesh>, Data>>" | Directional derivative of a P1 grid function. |
 */

namespace Rodin::Variational
{
  /**
   * @ingroup RodinVariational
   * @brief Base class for directional derivative operators.
   *
   * DerivativeBase provides the foundation for computing directional derivatives
   * of scalar functions in specified directions.
   *
   * @tparam Operand Type of the function being differentiated
   * @tparam Derived Derived class (CRTP pattern)
   */
  template <class Operand, class Derived>
  class DerivativeBase;

  /**
   * @ingroup GradSpecializations
   */
  template <class FES, class Data, class Derived>
  class DerivativeBase<GridFunction<FES, Data>, Derived>
    : public ScalarFunctionBase<
        typename FormLanguage::Traits<FES>::ScalarType, DerivativeBase<GridFunction<FES, Data>, Derived>>
  {
    public:
      /// @brief Finite element space type.
      using FESType = FES;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      /// @brief Operand type.
      using OperandType = GridFunction<FESType, Data>;

      /// @brief Parent class type.
      using Parent = ScalarFunctionBase<ScalarType, DerivativeBase<OperandType, Derived>>;

      /// @brief Constructs the expression from its operand.
      DerivativeBase(const OperandType& u)
        : m_u(u)
      {
        assert(u.getFiniteElementSpace().getVectorDimension() == 1);
      }

      /**
       * @brief Copy constructor
       */
      DerivativeBase(const DerivativeBase& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      /**
       * @brief Move constructor
       */
      DerivativeBase(DerivativeBase&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u))
      {}

      constexpr
      size_t getDimension() const
      {
        return m_u.get().getFiniteElementSpace().getMesh().getSpaceDimension();
      }

    public:
      /**
       * @brief Evaluates the partial derivative at a Point.
       *
       * Resolves mesh ownership and dispatches to the derived class's
       * @c interpolate. Falls back to inclusion / submesh restriction
       * when the polytope's mesh is not the FES mesh.
       */
      ScalarType getValue(const Geometry::Point& p) const
      {
        const auto& fes = getOperand().getFiniteElementSpace();
        const auto& fesMesh = fes.getMesh();

        ScalarType value = ScalarType(0);
        if (fesMesh.isLocalPoint(p))
        {
          this->interpolate(value, p);
        }
        else if (const auto inclusion = fesMesh.inclusion(p))
        {
          this->interpolate(value, *inclusion);
        }
        else if (fesMesh.isSubMesh())
        {
          const auto& submesh = fesMesh.asSubMesh();
          const auto restriction = submesh.restriction(p);
          this->interpolate(value, *restriction);
        }
        else
        {
          assert(false);
        }
        return value;
      }

      /**
       * @brief Evaluates the partial derivative at an IntegrationPoint.
       *
       * If the polytope is owned by the FES mesh, dispatches to
       * @c interpolate(out, ip). Otherwise falls back to inclusion / submesh
       * restriction.
       */
      ScalarType getValue(const IntegrationPoint& ip) const
      {
        const auto& p = ip.getPoint();
        const auto& fes = getOperand().getFiniteElementSpace();
        const auto& fesMesh = fes.getMesh();

        if (fesMesh.isLocalPoint(p))
        {
          ScalarType value = ScalarType(0);
          this->interpolate(value, ip);
          return value;
        }

        ScalarType value = ScalarType(0);
        if (const auto inclusion = fesMesh.inclusion(p))
        {
          this->interpolate(value, *inclusion);
        }
        else if (fesMesh.isSubMesh())
        {
          const auto& submesh = fesMesh.asSubMesh();
          const auto restriction = submesh.restriction(p);
          this->interpolate(value, *restriction);
        }
        else
        {
          assert(false);
        }
        return value;
      }

      /**
       * @brief Interpolation function to be overriden in Derived type.
       */
      constexpr
      void interpolate(ScalarType& out, const Geometry::Point& p) const
      {
        static_cast<const Derived&>(*this).interpolate(out, p);
      }

      constexpr
      /// @brief Interpolates at an integration point.
      void interpolate(ScalarType& out, const IntegrationPoint& ip) const
      {
        if constexpr (requires (const Derived& f, ScalarType& r, const IntegrationPoint& q) { f.interpolate(r, q); })
          static_cast<const Derived&>(*this).interpolate(out, ip);
        else
          static_cast<const Derived&>(*this).interpolate(out, ip.getPoint());
      }

      constexpr
      /// @brief Gets the operand function.
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      /**
       * @brief Copy function to be overriden in Derived type.
       */
      DerivativeBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
  };

  /// @brief Partial derivative of a shape function.
  template <class NestedDerived, class FES, ShapeFunctionSpaceType SpaceType>
  class Derivative<ShapeFunction<NestedDerived, FES, SpaceType>> final
    : public ShapeFunctionBase<Derivative<ShapeFunction<NestedDerived, FES, SpaceType>>>
  {
    public:
      /// Finite element space type
      using FESType = FES;
      /// @brief Shape function space the expression belongs to, trial or test.
      static constexpr ShapeFunctionSpaceType Space = SpaceType;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      /// Operand type
      using OperandType = ShapeFunction<NestedDerived, FESType, Space>;

      /// Parent class
      using Parent = ShapeFunctionBase<Derivative<OperandType>, FESType, Space>;

      Derivative(size_t i, const OperandType& u)
        : Parent(u.getFiniteElementSpace()),
          m_i(i),
          m_u(u)
      {}

      /// @brief Copy constructor.
      Derivative(const Derivative& other)
        : Parent(other),
          m_i(other.m_i),
          m_u(other.m_u)
      {}

      /// @brief Move constructor.
      Derivative(Derivative&& other)
        : Parent(std::move(other)),
          m_i(other.m_i),
          m_u(std::move(other.m_u))
      {}

      constexpr
      /// @brief Gets the operand function.
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      constexpr
      /// @brief Gets the operand in the shape function expression.
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }

      constexpr
      /// @brief Gets the global DOF indices for a polytope.
      size_t getDOFs(const Geometry::Polytope& element) const
      {
        return getOperand().getDOFs(element);
      }

      /// @brief Gets the integration point the expression is evaluated at.
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_ip);
        return *m_ip;
      }

      // Derivative& setIntegrationPoint(const Geometry::Point& p)
      // {
      //   m_ip = &ip;
      //   const auto& polytope = p.getPolytope();
      //   const size_t d = polytope.getDimension();
      //   const Index i = polytope.getIndex();
      //   const auto& fes = this->getFiniteElementSpace();
      //   const auto& fe = fes.getFiniteElement(d, i);
      //   const auto& rc = p.getReferenceCoordinates();
      //   const size_t dofs = this->getDOFs(p.getPolytope());
      //   m_gradients.resize(dofs);
      //   for (size_t local = 0; local < dofs; local++)
      //     m_gradients[local] = p.getJacobianInverse().transpose() * fe.getGradient(local)(rc);
      //   return *this;
      // }

      /// @brief Gets the basis function of a local degree of freedom.
      decltype(auto) getBasis(size_t local) const
      {
        return m_gradients[local](m_i);
      }

      Derivative* copy() const noexcept override
      {
        return new Derivative(*this);
      }

    private:
      size_t m_i;
      std::reference_wrapper<const OperandType> m_u;

      const IntegrationPoint* m_ip;

      std::vector<Math::SpatialVector<Real>> m_gradients;
  };

  /// @brief Deduction guide for @c Derivative.
  template <class NestedDerived, class FES, ShapeFunctionSpaceType SpaceType>
  Derivative(size_t i, const ShapeFunction<NestedDerived, FES, SpaceType>& u)
    -> Derivative<ShapeFunction<NestedDerived, FES, SpaceType>>;

  /**
   * @brief %Utility function for computing @f$ \partial_x u @f$
   * @param[in] u GridFunction instance
   *
   * Given a scalar function @f$ u : \mathbb{R}^s \rightarrow \mathbb{R} @f$,
   * this function constructs the derivative in the @f$ x @f$ direction
   * @f$
   *   \dfrac{\partial u}{\partial x}
   * @f$
   */
  template <class Operand>
  auto Dx(const Operand& u)
  {
    return Derivative(0, u);
  }

  /**
   * @brief %Utility function for computing @f$ \partial_y u @f$
   * @param[in] u GridFunction instance
   *
   * Given a scalar function @f$ u : \mathbb{R}^s \rightarrow \mathbb{R} @f$,
   * this function constructs the derivative in the @f$ y @f$ direction
   * @f$
   *   \dfrac{\partial u}{\partial y}
   * @f$
   */
  template <class Operand>
  auto Dy(const Operand& u)
  {
    return Derivative(1, u);
  }

  /**
   * @brief %Utility function for computing @f$ \partial_z u @f$
   * @param[in] u GridFunction instance
   *
   * Given a scalar function @f$ u : \mathbb{R}^s \rightarrow \mathbb{R} @f$,
   * this function constructs the derivative in the @f$ y @f$ direction
   * @f$
   *   \dfrac{\partial u}{\partial z}
   * @f$
   */
  template <class Operand>
  auto Dz(const Operand& u)
  {
    return Derivative(2, u);
  }
}

#endif
