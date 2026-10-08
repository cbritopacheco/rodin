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
#include "Grad.h"
#include "Rodin/Variational/IntegrationPoint.h"
#include "ShapeFunction.h"

namespace Rodin::FormLanguage
{
  /**
   * @brief Type traits for @c Derivative over a shape function: exposes the shape
   * function space, the finite element space, the scalar type and the operand type.
   */
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

namespace Rodin::Variational
{
  /**
   * @defgroup DerivativeSpecializations Derivative Template Specializations
   * @brief Template specializations of the Derivative class.
   * @see @ref Derivative
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref Derivative "Derivative<GridFunction<MatrixFES, Data>>" | Matrix-valued grid-function operator for all supported spaces and backends. |
   * | @ref Derivative "Derivative<ShapeFunction<Derived, MatrixFES, Space>>" | Matrix shape-function operator with the scalar family's geometry and trace semantics. |
   * | @ref Derivative "Derivative<H1<K, Scalar, Mesh>, ShapeFunction<NestedDerived, H1<K, Scalar, Mesh>, Space>>" | Directional derivative of an H1 shape function. |
   * | @ref Derivative "Derivative<P1<Range, Mesh>, GridFunction<P1<Range, Mesh>, Data>>" | Directional derivative of a P1 grid function. |
   */

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
   * @brief CRTP base for the partial derivative of a grid function.
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

      /**
       * @brief Constructs the expression from its operand.
       * @param u Operand expression.
       */
      DerivativeBase(const OperandType& u)
        : m_u(u)
      {
        assert(u.getFiniteElementSpace().getVectorDimension() == 1);
      }

      /**
       * @brief Copy constructor
       * @param other Object to copy from.
       */
      DerivativeBase(const DerivativeBase& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      /**
       * @brief Move constructor
       * @param other Object to move from.
       */
      DerivativeBase(DerivativeBase&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u))
      {}

      /**
       * @brief Gets the topological dimension.
       * @returns The topological dimension.
       */
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
       * @param p Point at which the operation is evaluated.
       * @returns Value of the expression at the supplied evaluation point.
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
       * @param ip Integration point at which the expression is evaluated.
       * @returns Value of the expression at the supplied evaluation point.
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
       * @param out Storage for the computed result.
       * @param p Point at which the operation is evaluated.
       */
      constexpr
      void interpolate(ScalarType& out, const Geometry::Point& p) const
      {
        static_cast<const Derived&>(*this).interpolate(out, p);
      }

      /**
       * @brief Interpolates at an integration point.
       * @param out Storage for the computed result.
       * @param ip Integration point at which the expression is evaluated.
       */
      constexpr
      void interpolate(ScalarType& out, const IntegrationPoint& ip) const
      {
        if constexpr (requires (const Derived& f, ScalarType& r, const IntegrationPoint& q) { f.interpolate(r, q); })
          static_cast<const Derived&>(*this).interpolate(out, ip);
        else
          static_cast<const Derived&>(*this).interpolate(out, ip.getPoint());
      }

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
       * @brief Copy function to be overriden in Derived type.
       * @returns Pointer to a newly allocated copy; the caller owns the returned object.
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
    requires(
      !FormLanguage::IsMatrixRange<typename FormLanguage::Traits<FES>::RangeType>::Value)
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

      /**
       * @brief Constructs the partial derivative of an operand along a direction.
       * @param i Index of the requested entry.
       * @param u Operand expression.
       */
      Derivative(size_t i, const OperandType& u)
        : Parent(u.getFiniteElementSpace()),
          m_i(i),
          m_u(u)
      {}

      /**
       * @brief Copy constructor.
       * @param other Object to copy from.
       */
      Derivative(const Derivative& other)
        : Parent(other),
          m_i(other.m_i),
          m_u(other.m_u)
      {}

      /**
       * @brief Move constructor.
       * @param other Object to move from.
       */
      Derivative(Derivative&& other)
        : Parent(std::move(other)),
          m_i(other.m_i),
          m_u(std::move(other.m_u))
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
       * @brief Gets the operand in the shape function expression.
       * @returns The operand in the shape function expression.
       */
      constexpr
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
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

      /**
       * @brief Gets the basis function of a local degree of freedom.
       * @param local Index in the local numbering.
       * @returns Value of the selected local basis function at the evaluation point.
       */
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

  /**
   * @brief Deduction guide for @c Derivative.
   * @param i Index of the requested entry.
   * @param u Shape-function operand.
   */
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
   * @returns Expression for the derivative in the first coordinate direction.
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
   * @returns Expression for the derivative in the second coordinate direction.
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
   * @returns Expression for the derivative in the third coordinate direction.
   */
  template <class Operand>
  auto Dz(const Operand& u)
  {
    return Derivative(2, u);
  }
}

namespace Rodin::FormLanguage
{
  /// @brief Type traits for the matrix or tensor expression specialization.
  template <class FES, class Data>
    requires IsMatrixRange<typename Traits<FES>::RangeType>::Value
  struct Traits<Variational::Derivative<Variational::GridFunction<FES, Data>>>
  {
      /// @brief Finite element space type.
      using FESType = FES;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename Traits<FES>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialMatrix<ScalarType>;
  };
  /// @brief Type traits for the matrix or tensor expression specialization.
  template <class Derived, class FES, Variational::ShapeFunctionSpaceType Space>
    requires IsMatrixRange<typename Traits<FES>::RangeType>::Value
  struct Traits<Variational::Derivative<Variational::ShapeFunction<Derived, FES, Space>>>
  {
      /// @brief Finite element space type.
      using FESType = FES;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename Traits<FES>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialMatrix<ScalarType>;
      /// @brief Trial or test shape-function space.
      static constexpr auto SpaceType = Space;
  };
}
namespace Rodin::Variational
{
  /**
   * @ingroup DerivativeSpecializations
   * @brief Matrix-field derivative, @f$ D_k A_{ij}=\partial_k A_{ij} @f$.
   */
  template <class FES, class Data>
    requires FormLanguage::IsMatrixRange<
      typename FormLanguage::Traits<FES>::RangeType>::Value
  class Derivative<GridFunction<FES, Data>> final
    : public FunctionBase<Derivative<GridFunction<FES, Data>>>
  {
    public:
      /// @brief CRTP or finite element base class.
      using Parent = FunctionBase<Derivative>;
      /// @brief Cloned or referenced expression operand type.
      using OperandType = GridFunction<FES, Data>;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialMatrix<ScalarType>;
      /**
       * @brief Constructs a derivative along the selected ambient coordinate.
       * @param direction Direction in which the derivative is evaluated.
       * @param operand Operand expression.
       */
      Derivative(size_t direction, const OperandType& operand)
        : m_gradient(operand),
          m_direction(direction)
      {}
      /**
       * @brief Constructs a derivative along the selected ambient coordinate.
       * @param other Object to copy from.
       */
      Derivative(const Derivative& other)
        : Parent(other),
          m_gradient(other.m_gradient),
          m_direction(other.m_direction)
      {}
      /**
       * @brief Constructs a derivative along the selected ambient coordinate.
       * @param other Object to move from.
       */
      Derivative(Derivative&& other)
        : Parent(std::move(other)),
          m_gradient(std::move(other.m_gradient)),
          m_direction(other.m_direction)
      {}
      /**
       * @brief Returns the differentiated or indexed operand.
       * @returns The differentiated or indexed operand.
       */
      const OperandType& getOperand() const
      {
        return m_gradient.getOperand();
      }
      /**
       * @brief Evaluates the expression at the supplied physical or integration point.
       * @param point Point at which the operation is evaluated.
       * @returns Value of the expression at the supplied evaluation point.
       */
      template <class Point>
      RangeType getValue(const Point& point) const
      {
        auto gradientExpression = m_gradient;
        gradientExpression.traceOf(this->getTraceDomain());
        const auto gradient = gradientExpression.getValue(point);
        if (m_direction >= gradient.getDimension(2))
          Alert::Exception()
            << "Partial derivative direction exceeds the spatial dimension."
            << Alert::Raise;
        RangeType value(gradient.getDimension(0), gradient.getDimension(1));
        for (size_t row = 0; row < value.rows(); ++row)
          for (size_t col = 0; col < value.cols(); ++col)
            value(row, col) = gradient(row, col, m_direction);
        return value;
      }
      /**
       * @brief Returns the polynomial order when it is known.
       * @param poly Mesh entity used by this operation.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       */
      Optional<size_t> getOrder(const Geometry::Polytope& poly) const noexcept
      {
        return m_gradient.getOrder(poly);
      }
      Derivative* copy() const noexcept override
      {
        return new Derivative(*this);
      }

    private:
      Grad<OperandType> m_gradient;
      size_t m_direction;
  };

  /**
   * @ingroup DerivativeSpecializations
   * @brief Matrix-basis derivative, @f$ D_k A_{ij}=\partial_k A_{ij} @f$.
   */
  template <class Derived, class FES, ShapeFunctionSpaceType Space>
    requires FormLanguage::IsMatrixRange<
      typename FormLanguage::Traits<FES>::RangeType>::Value
  class Derivative<ShapeFunction<Derived, FES, Space>> final
    : public ShapeFunctionBase<Derivative<ShapeFunction<Derived, FES, Space>>, FES, Space>
  {
    public:
      /// @brief CRTP or finite element base class.
      using Parent = ShapeFunctionBase<Derivative, FES, Space>;
      /// @brief Cloned or referenced expression operand type.
      using OperandType = ShapeFunction<Derived, FES, Space>;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialMatrix<ScalarType>;
      /**
       * @brief Constructs a derivative along the selected ambient coordinate.
       * @param direction Direction in which the derivative is evaluated.
       * @param operand Operand expression.
       */
      Derivative(size_t direction, const OperandType& operand)
        : Parent(operand.getFiniteElementSpace()),
          m_gradient(operand),
          m_direction(direction)
      {}
      /**
       * @brief Constructs a derivative along the selected ambient coordinate.
       * @param other Object to copy from.
       */
      Derivative(const Derivative& other)
        : Parent(other),
          m_gradient(other.m_gradient),
          m_direction(other.m_direction)
      {}
      /**
       * @brief Constructs a derivative along the selected ambient coordinate.
       * @param other Object to move from.
       */
      Derivative(Derivative&& other)
        : Parent(std::move(other)),
          m_gradient(std::move(other.m_gradient)),
          m_direction(other.m_direction)
      {}
      /**
       * @brief Returns the differentiated or indexed operand.
       * @returns The differentiated or indexed operand.
       */
      const OperandType& getOperand() const
      {
        return m_gradient.getOperand();
      }
      /**
       * @brief Returns the leaf shape function used for assembly.
       * @returns The leaf shape function used for assembly.
       */
      const auto& getLeaf() const
      {
        return m_gradient.getLeaf();
      }
      /**
       * @brief Returns the local basis count for the selected polytope.
       * @param poly Mesh entity used by this operation.
       * @returns Number of local basis functions on the selected entity.
       */
      size_t getDOFs(const Geometry::Polytope& poly) const
      {
        return m_gradient.getDOFs(poly);
      }
      /**
       * @brief Returns the currently bound integration point.
       * @returns The currently bound integration point.
       */
      const IntegrationPoint& getIntegrationPoint() const
      {
        return m_gradient.getIntegrationPoint();
      }
      /**
       * @brief Binds the integration point and prepares local basis values.
       * @param point Point at which the operation is evaluated.
       * @returns Reference to this object after the operation.
       */
      Derivative& setIntegrationPoint(const IntegrationPoint& point)
      {
        m_gradient.setIntegrationPoint(point);
        return *this;
      }
      /**
       * @brief Returns a basis value at the bound integration point.
       * @param local Index in the local numbering.
       * @returns Value of the selected local basis function at the evaluation point.
       */
      RangeType getBasis(size_t local) const
      {
        const auto gradient = m_gradient.getBasis(local);
        if (m_direction >= gradient.getDimension(2))
          Alert::Exception()
            << "Partial derivative direction exceeds the spatial dimension."
            << Alert::Raise;
        RangeType value(gradient.getDimension(0), gradient.getDimension(1));
        for (size_t row = 0; row < value.rows(); ++row)
          for (size_t col = 0; col < value.cols(); ++col)
            value(row, col) = gradient(row, col, m_direction);
        return value;
      }
      /**
       * @brief Returns the polynomial order when it is known.
       * @param poly Mesh entity used by this operation.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       */
      Optional<size_t> getOrder(const Geometry::Polytope& poly) const noexcept
      {
        return m_gradient.getOrder(poly);
      }
      Derivative* copy() const noexcept override
      {
        return new Derivative(*this);
      }

    private:
      Grad<OperandType> m_gradient;
      size_t m_direction;
  };
  /**
   * @brief Deduces the matrix space or coefficient type from constructor arguments.
   * @param direction Direction in which the derivative is evaluated.
   * @param operand Operand expression.
   */
  template <class FES, class Data>
    requires FormLanguage::IsMatrixRange<
               typename FormLanguage::Traits<FES>::RangeType>::Value
  Derivative(size_t direction,
    const GridFunction<FES, Data>& operand) -> Derivative<GridFunction<FES, Data>>;
  /**
   * @brief Deduces the matrix space or coefficient type from constructor arguments.
   * @param direction Direction in which the derivative is evaluated.
   * @param operand Operand expression.
   */
  template <class Derived, class FES, ShapeFunctionSpaceType Space>
    requires FormLanguage::IsMatrixRange<
               typename FormLanguage::Traits<FES>::RangeType>::Value
  Derivative(size_t direction, const ShapeFunction<Derived, FES, Space>& operand)
    -> Derivative<ShapeFunction<Derived, FES, Space>>;
}

#endif
