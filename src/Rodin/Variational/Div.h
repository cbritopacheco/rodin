/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Div.h
 * @brief Divergence operator for vector-valued functions.
 *
 * This file defines the Div class, which computes the divergence of
 * vector-valued functions in variational formulations. The divergence is
 * a fundamental differential operator that measures the "outflow" of a vector field.
 *
 * ## Mathematical Foundation
 * For a vector field @f$ \mathbf{u} : \Omega \subset \mathbb{R}^d \to \mathbb{R}^d @f$,
 * the divergence is defined as:
 * @f[
 *   \nabla \cdot \mathbf{u} = \sum_{i=1}^d \frac{\partial u_i}{\partial x_i}
 * @f]
 *
 * ## Applications
 * - Incompressibility constraint: @f$ \nabla \cdot \mathbf{u} = 0 @f$
 * - Conservation laws: @f$ \nabla \cdot \mathbf{F} = 0 @f$
 * - Mixed formulations for elliptic problems
 */
#ifndef RODIN_VARIATIONAL_DIV_H
#define RODIN_VARIATIONAL_DIV_H

#include "ForwardDecls.h"
#include "Grad.h"

#include "Jacobian.h"
#include "GridFunction.h"
#include "Rodin/Types.h"
#include "TestFunction.h"
#include "TrialFunction.h"
#include "RealFunction.h"

namespace Rodin::Variational
{
  /**
    * @defgroup DivSpecializations Div Template Specializations
    * @brief Template specializations of the Div class.
    * @see @ref Div
    *
    * | Specialization | Description |
    * |----------------|-------------|
   * | @ref Div "Div<GridFunction<MatrixFES, Data>>" | Matrix-valued grid-function operator for all supported spaces and backends. |
   * | @ref Div "Div<ShapeFunction<Derived, MatrixFES, Space>>" | Matrix shape-function operator with the scalar family's geometry and trace semantics. |
    * | @ref DivBase "DivBase<GridFunction<FES, Data>, Derived>" | Generic divergence base for vector-valued grid functions. |
    * | @ref Div "Div<P0g<Scalar, Mesh>, GridFunction<P0g<Scalar, Mesh>, Data>>" | Divergence of a discontinuous P0g grid function. |
    * | @ref Div "Div<P0g<Scalar, Mesh>, ShapeFunction<NestedDerived, P0g<Scalar, Mesh>, Space>>" | Divergence of a P0g shape-function expression. |
    * | @ref Div "Div<H1<K, Scalar, Mesh>, GridFunction<H1<K, Scalar, Mesh>, Data>>" | Divergence of an H1 grid function. |
    * | @ref Div "Div<H1<K, Scalar, Mesh>, ShapeFunction<NestedDerived, H1<K, Scalar, Mesh>, Space>>" | Divergence of an H1 shape-function expression. |
    * | @ref Div "Div<P1<Scalar, Mesh>, GridFunction<P1<Scalar, Mesh>, Data>>" | Divergence of a P1 grid function. |
    * | @ref Div "Div<P1<Scalar, Mesh>, ShapeFunction<NestedDerived, P1<Scalar, Mesh>, Space>>" | Divergence of a P1 shape-function expression. |
    */

  /**
   * @ingroup RodinVariational
   * @brief Base class for divergence operator implementations.
   *
   * DivBase provides the foundation for computing divergences of vector-valued
   * functions. The divergence operator maps vector fields to scalar fields.
   *
   * @tparam Operand Type of the vector function
   * @tparam Derived Derived class (CRTP pattern)
   */
  template <class Operand, class Derived>
  class DivBase;

  /**
   * @ingroup DivSpecializations
   * @brief Divergence of a P1 GridFunction
   */
  template <class FES, class Data, class Derived>
  class DivBase<GridFunction<FES, Data>, Derived>
    : public ScalarFunctionBase<typename FormLanguage::Traits<FES>::ScalarType, DivBase<GridFunction<FES, Data>, Derived>>
  {
    public:
      /// @brief Finite element space type.
      using FESType = FES;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      /// @brief Operand type.
      using OperandType = GridFunction<FES, Data>;

      /// Parent class
      using Parent = ScalarFunctionBase<ScalarType, DivBase<OperandType, Derived>>;

      /**
       * @brief Constructs the Div of a @f$ \mathbb{P}_1 @f$ function @f$ u
       * @f$.
       * @param[in] u P1 GridFunction
       */
      DivBase(const OperandType& u)
        : m_u(u)
      {}

      /**
       * @brief Copy constructor
       * @param other Object to copy from.
       */
      DivBase(const DivBase& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      /**
       * @brief Move constructor
       * @param other Object to move from.
       */
      DivBase(DivBase&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u))
      {}

    public:
      /**
       * @brief Evaluates the divergence at a Point.
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
       * @brief Evaluates the divergence at an IntegrationPoint.
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
       * @brief Gets the operand grid function.
       * @return Reference to the vector-valued grid function
       */
      constexpr
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      /**
       * @brief Interpolates the divergence at a point (to be overridden in derived class).
       * @param[out] out Output scalar for divergence result
       * @param[in] p Point at which to interpolate
       *
       * This virtual function is overridden in derived classes (e.g., P1::Div)
       * to provide finite element-specific divergence computation.
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
       * @brief Returns the polynomial order used on a mesh entity.
       * @param poly Mesh entity used by this operation.
       * @returns Polynomial order on the entity, or an empty optional when no order is available.
       */
      Optional<size_t> getOrder(const Geometry::Polytope& poly) const noexcept
      {
        return static_cast<const Derived&>(*this).getOrder(poly);
      }

      /**
       * @brief Copy function to be overriden in Derived type.
       * @returns Pointer to a newly allocated copy; the caller owns the returned object.
       */
      DivBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
  };
}

namespace Rodin::FormLanguage
{
  /// @brief Type traits for the matrix or tensor expression specialization.
  template <class FES, class Data>
    requires IsMatrixRange<typename Traits<FES>::RangeType>::Value
  struct Traits<Variational::Div<Variational::GridFunction<FES, Data>>>
  {
      /// @brief Finite element space type.
      using FESType = FES;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename Traits<FES>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialVector<ScalarType>;
  };
  /// @brief Type traits for the matrix or tensor expression specialization.
  template <class Derived, class FES, Variational::ShapeFunctionSpaceType Space>
    requires IsMatrixRange<typename Traits<FES>::RangeType>::Value
  struct Traits<Variational::Div<Variational::ShapeFunction<Derived, FES, Space>>>
  {
      /// @brief Finite element space type.
      using FESType = FES;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename Traits<FES>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialVector<ScalarType>;
      /// @brief Trial or test shape-function space.
      static constexpr auto SpaceType = Space;
  };
}
namespace Rodin::Variational
{
  /**
   * @ingroup DivSpecializations
   * @brief Matrix-field div, @f$ (\operatorname{div}A)_i=\sum_j\partial_j A_{ij} @f$.
   */
  template <class FES, class Data>
    requires FormLanguage::IsMatrixRange<
      typename FormLanguage::Traits<FES>::RangeType>::Value
  class Div<GridFunction<FES, Data>> final
    : public FunctionBase<Div<GridFunction<FES, Data>>>
  {
    public:
      /// @brief CRTP or finite element base class.
      using Parent = FunctionBase<Div>;
      /// @brief Cloned or referenced expression operand type.
      using OperandType = GridFunction<FES, Data>;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialVector<ScalarType>;
      /**
       * @brief Constructs row-wise divergence of a matrix field.
       * @param operand Operand expression.
       */
      Div(const OperandType& operand)
        : m_gradient(operand)
      {}
      /**
       * @brief Constructs row-wise divergence of a matrix field.
       * @param other Object to copy from.
       */
      Div(const Div& other)
        : Parent(other),
          m_gradient(other.m_gradient)
      {}
      /**
       * @brief Constructs row-wise divergence of a matrix field.
       * @param other Object to move from.
       */
      Div(Div&& other)
        : Parent(std::move(other)),
          m_gradient(std::move(other.m_gradient))
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
        RangeType value(gradient.getDimension(0));
        value.setZero();
        if (gradient.getDimension(1) != gradient.getDimension(2))
          Alert::Exception()
            << "Matrix divergence requires columns equal to the spatial dimension."
            << Alert::Raise;
        for (size_t row = 0; row < gradient.getDimension(0); ++row)
          for (size_t k = 0; k < gradient.getDimension(2); ++k)
            value(row) += gradient(row, k, k);
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
      Div* copy() const noexcept override
      {
        return new Div(*this);
      }

    private:
      Grad<OperandType> m_gradient;
  };

  /**
   * @ingroup DivSpecializations
   * @brief Matrix-basis div, @f$ (\operatorname{div}A)_i=\sum_j\partial_j A_{ij} @f$.
   */
  template <class Derived, class FES, ShapeFunctionSpaceType Space>
    requires FormLanguage::IsMatrixRange<
      typename FormLanguage::Traits<FES>::RangeType>::Value
  class Div<ShapeFunction<Derived, FES, Space>> final
    : public ShapeFunctionBase<Div<ShapeFunction<Derived, FES, Space>>, FES, Space>
  {
    public:
      /// @brief CRTP or finite element base class.
      using Parent = ShapeFunctionBase<Div, FES, Space>;
      /// @brief Cloned or referenced expression operand type.
      using OperandType = ShapeFunction<Derived, FES, Space>;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialVector<ScalarType>;
      /**
       * @brief Constructs row-wise divergence of a matrix field.
       * @param operand Operand expression.
       */
      Div(const OperandType& operand)
        : Parent(operand.getFiniteElementSpace()),
          m_gradient(operand)
      {}
      /**
       * @brief Constructs row-wise divergence of a matrix field.
       * @param other Object to copy from.
       */
      Div(const Div& other)
        : Parent(other),
          m_gradient(other.m_gradient)
      {}
      /**
       * @brief Constructs row-wise divergence of a matrix field.
       * @param other Object to move from.
       */
      Div(Div&& other)
        : Parent(std::move(other)),
          m_gradient(std::move(other.m_gradient))
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
      Div& setIntegrationPoint(const IntegrationPoint& point)
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
        RangeType value(gradient.getDimension(0));
        value.setZero();
        if (gradient.getDimension(1) != gradient.getDimension(2))
          Alert::Exception()
            << "Matrix divergence requires columns equal to the spatial dimension."
            << Alert::Raise;
        for (size_t row = 0; row < gradient.getDimension(0); ++row)
          for (size_t k = 0; k < gradient.getDimension(2); ++k)
            value(row) += gradient(row, k, k);
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
      Div* copy() const noexcept override
      {
        return new Div(*this);
      }

    private:
      Grad<OperandType> m_gradient;
  };
  /// @brief Deduces the matrix space or coefficient type from constructor arguments.
  template <class FES, class Data>
    requires FormLanguage::IsMatrixRange<
               typename FormLanguage::Traits<FES>::RangeType>::Value
  Div(const GridFunction<FES, Data>&) -> Div<GridFunction<FES, Data>>;
  /// @brief Deduces the matrix space or coefficient type from constructor arguments.
  template <class Derived, class FES, ShapeFunctionSpaceType Space>
    requires FormLanguage::IsMatrixRange<
               typename FormLanguage::Traits<FES>::RangeType>::Value
  Div(
    const ShapeFunction<Derived, FES, Space>&) -> Div<ShapeFunction<Derived, FES, Space>>;
}

#endif
