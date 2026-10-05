/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Jacobian.h
 * @brief Jacobian matrix operator for vector-valued functions.
 *
 * This file defines the Jacobian class, which computes the Jacobian matrix
 * (matrix of all first-order partial derivatives) of vector-valued functions
 * in variational formulations.
 *
 * ## Mathematical Foundation
 * For a vector-valued function @f$ \mathbf{u} : \Omega \subset \mathbb{R}^d \to \mathbb{R}^n @f$,
 * the Jacobian matrix is defined as:
 * @f[
 *   J_{ij} = \frac{\partial u_i}{\partial x_j}
 * @f]
 * resulting in an @f$ n \times d @f$ matrix.
 *
 * ## Special Cases
 * - When @f$ n = d @f$, the determinant @f$ \det(J) @f$ appears in change of variables
 * - For @f$ d = n = 2,3 @f$, related to deformation gradient in mechanics
 * - Transpose @f$ J^T @f$ gives the gradient of components
 *
 * ## Applications
 * - Nonlinear elasticity: deformation gradient tensor
 * - Fluid dynamics: velocity gradient tensor
 * - Differential geometry: metric tensor computations
 * - Coordinate transformations
 *
 * ## Usage Example
 * ```cpp
 * // Displacement gradient for elasticity
 * P1 Vh(mesh, mesh.getSpaceDimension());
 * GridFunction<P1> u(Vh);
 * auto F = Jacobian(u);  // Deformation gradient F = \nabla u
 * ```
 */
#ifndef RODIN_VARIATIONAL_JACOBIAN_H
#define RODIN_VARIATIONAL_JACOBIAN_H

#include "ForwardDecls.h"
#include "Grad.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "MatrixFunction.h"
#include "IntegrationPoint.h"

namespace Rodin::Variational
{
  /**
   * @defgroup JacobianSpecializations Jacobian Template Specializations
   * @brief Template specializations of the Jacobian class.
   * @see @ref Jacobian
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref Jacobian "Jacobian<GridFunction<MatrixFES, Data>>" | Matrix-valued grid-function operator for all supported spaces and backends. |
   * | @ref Jacobian "Jacobian<ShapeFunction<Derived, MatrixFES, Space>>" | Matrix shape-function operator with the scalar family's geometry and trace semantics. |
   * | @ref JacobianBase "JacobianBase<GridFunction<FES, Data>, Derived>" | Generic Jacobian base for vector-valued grid functions. |
   * | @ref Jacobian "Jacobian<P0g<Scalar, Mesh>, GridFunction<P0g<Scalar, Mesh>, Data>>" | Jacobian of a discontinuous P0g grid function. |
   * | @ref Jacobian "Jacobian<P0g<Scalar, Mesh>, ShapeFunction<NestedDerived, P0g<Scalar, Mesh>, Space>>" | Jacobian of a P0g shape-function expression. |
   * | @ref Jacobian "Jacobian<H1<K, Scalar, Mesh>, GridFunction<H1<K, Scalar, Mesh>, Data>>" | Jacobian of an H1 grid function. |
   * | @ref Jacobian "Jacobian<H1<K, Scalar, Mesh>, ShapeFunction<NestedDerived, H1<K, Scalar, Mesh>, Space>>" | Jacobian of an H1 shape-function expression. |
   * | @ref Jacobian "Jacobian<P1<Range, Mesh>, GridFunction<P1<Range, Mesh>, Data>>" | Jacobian of a P1 grid function. |
   * | @ref Jacobian "Jacobian<P1<Range, Mesh>, ShapeFunction<NestedDerived, P1<Range, Mesh>, Space>>" | Jacobian of a P1 shape-function expression. |
   */

  /**
   * @ingroup RodinVariational
   * @brief Base class for Jacobian matrix operator implementations.
   *
   * JacobianBase provides the foundation for computing Jacobian matrices of
   * vector-valued functions.
   *
   * @tparam Operand Type of the vector function
   * @tparam Derived Derived class (CRTP pattern)
   */
  template <class Operand, class Derived>
  class JacobianBase;

  /**
   * @ingroup JacobianSpecializations
   * @brief Jacobian of a P1 GridFunction
   */
  template <class FES, class Data, class Derived>
  class JacobianBase<GridFunction<FES, Data>, Derived>
    : public MatrixFunctionBase<
        typename FormLanguage::Traits<FES>::ScalarType, JacobianBase<GridFunction<FES, Data>, Derived>>
  {
    public:
      /// @brief Finite element space type.
      using FESType = FES;

      /// @brief Scalar value type.
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      /// @brief Range (evaluation value) type.
      using RangeType = Math::SpatialMatrix<ScalarType>;

      /// @brief Small spatial matrix value type.
      using SpatialMatrixType = Math::SpatialMatrix<ScalarType>;

      /// @brief Operand type.
      using OperandType = GridFunction<FESType, Data>;

      /// @brief Parent class type.
      using Parent =
        MatrixFunctionBase<ScalarType, JacobianBase<OperandType, Derived>>;

      /**
       * @brief Constructs the Jacobian of a P1 grid function.
       * @param[in] u P1 GridFunction to differentiate
       *
       * Creates the Jacobian matrix operator @f$ J(\mathbf{u}) @f$ where
       * @f$ J_{ij} = \frac{\partial u_i}{\partial x_j} @f$.
       */
      JacobianBase(const OperandType& u)
        : m_u(u)
      {}

      /**
       * @brief Copy constructor.
       * @param[in] other Jacobian to copy
       */
      JacobianBase(const JacobianBase& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      /**
       * @brief Move constructor.
       * @param[in] other Jacobian to move from
       */
      JacobianBase(JacobianBase&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u))
      {}

      /**
       * @brief Gets the number of rows in the Jacobian matrix.
       * @return Number of components in the vector function
       *
       * For a vector function @f$ \mathbf{u}: \mathbb{R}^d \to \mathbb{R}^m @f$,
       * the Jacobian has @f$ m @f$ rows.
       */
      constexpr
      size_t getRows() const
      {
        return getOperand().getFiniteElementSpace().getVectorDimension();
      }

      /**
       * @brief Gets the number of columns in the Jacobian matrix.
       * @return Spatial dimension of the domain
       *
       * For a vector function @f$ \mathbf{u}: \mathbb{R}^d \to \mathbb{R}^m @f$,
       * the Jacobian has @f$ d @f$ columns.
       */
      constexpr
      size_t getColumns() const
      {
        return getOperand().getFiniteElementSpace().getMesh().getSpaceDimension();
      }

    public:
      /**
       * @brief Evaluates the Jacobian matrix at a Point.
       *
       * Resolves mesh ownership and dispatches to the derived class's
       * @c interpolate. Falls back to inclusion / submesh restriction
       * when the polytope's mesh is not the FES mesh.
       * @param p Point at which the operation is evaluated.
       * @returns Value of the expression at the supplied evaluation point.
       */
      SpatialMatrixType getValue(const Geometry::Point& p) const
      {
        const auto& fes = getOperand().getFiniteElementSpace();
        const auto& fesMesh = fes.getMesh();

        SpatialMatrixType value;
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
       * @brief Evaluates the Jacobian matrix at an IntegrationPoint.
       *
       * If the polytope is owned by the FES mesh, dispatches to
       * @c interpolate(out, ip). Otherwise falls back to inclusion / submesh
       * restriction when the polytope's mesh is not the FES mesh.
       * @param ip Integration point at which the expression is evaluated.
       * @returns Value of the expression at the supplied evaluation point.
       */
      SpatialMatrixType getValue(const IntegrationPoint& ip) const
      {
        const auto& p = ip.getPoint();
        const auto& fes = getOperand().getFiniteElementSpace();
        const auto& fesMesh = fes.getMesh();

        if (fesMesh.isLocalPoint(p))
        {
          SpatialMatrixType value;
          this->interpolate(value, ip);
          return value;
        }

        SpatialMatrixType value;
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
       * @brief Interpolates the Jacobian at a point (to be overridden in derived class).
       * @param out Storage for the computed result.
       * @param p Point at which the operation is evaluated.
       */
      constexpr
      void interpolate(SpatialMatrixType& out, const Geometry::Point& p) const
      {
        static_cast<const Derived&>(*this).interpolate(out, p);
      }

      /**
       * @brief Interpolates the Jacobian at an integration point.
       *
       * Forwards to the derived class's @c interpolate(IntegrationPoint)
       * if it provides one; otherwise falls back to the Point overload
       * via @c ip.getPoint().
       * @param out Storage for the computed result.
       * @param ip Integration point at which the expression is evaluated.
       */
      constexpr
      void interpolate(SpatialMatrixType& out, const IntegrationPoint& ip) const
      {
        if constexpr (requires (const Derived& f, SpatialMatrixType& r, const IntegrationPoint& q) { f.interpolate(r, q); })
          static_cast<const Derived&>(*this).interpolate(out, ip);
        else
          static_cast<const Derived&>(*this).interpolate(out, ip.getPoint());
      }

      /// @brief Returns the polynomial order used on a mesh entity.
      /// @param polytope Mesh entity used by this operation.
      /// @returns Polynomial order on the entity, or an empty optional when no order is available.
      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& polytope) const noexcept
      {
        return static_cast<const Derived&>(*this).getOrder(polytope);
      }

      /**
       * @brief Creates a polymorphic copy (to be overridden in derived class).
       * @return Pointer to a new copy
       */
      JacobianBase* copy() const noexcept override
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
  struct Traits<Variational::Jacobian<Variational::GridFunction<FES, Data>>>
  {
      /// @brief Finite element space type.
      using FESType = FES;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename Traits<FES>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialTensor<ScalarType>;
  };
  /// @brief Type traits for the matrix or tensor expression specialization.
  template <class Derived, class FES, Variational::ShapeFunctionSpaceType Space>
    requires IsMatrixRange<typename Traits<FES>::RangeType>::Value
  struct Traits<Variational::Jacobian<Variational::ShapeFunction<Derived, FES, Space>>>
  {
      /// @brief Finite element space type.
      using FESType = FES;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename Traits<FES>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialTensor<ScalarType>;
      /// @brief Trial or test shape-function space.
      static constexpr auto SpaceType = Space;
  };
}
namespace Rodin::Variational
{
  /**
   * @ingroup JacobianSpecializations
   * @brief Matrix-field jacobian, @f$ J_{ijk}=\partial_k A_{ij} @f$.
   */
  template <class FES, class Data>
    requires FormLanguage::IsMatrixRange<
      typename FormLanguage::Traits<FES>::RangeType>::Value
  class Jacobian<GridFunction<FES, Data>> final
    : public FunctionBase<Jacobian<GridFunction<FES, Data>>>
  {
    public:
      /// @brief CRTP or finite element base class.
      using Parent = FunctionBase<Jacobian>;
      /// @brief Cloned or referenced expression operand type.
      using OperandType = GridFunction<FES, Data>;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialTensor<ScalarType>;
      /// @brief Constructs a rank-three matrix Jacobian.
      /// @param operand Operand expression.
      Jacobian(const OperandType& operand)
        : m_gradient(operand)
      {}
      /// @brief Constructs a rank-three matrix Jacobian.
      /// @param other Object to copy from.
      Jacobian(const Jacobian& other)
        : Parent(other),
          m_gradient(other.m_gradient)
      {}
      /// @brief Constructs a rank-three matrix Jacobian.
      /// @param other Object to move from.
      Jacobian(Jacobian&& other)
        : Parent(std::move(other)),
          m_gradient(std::move(other.m_gradient))
      {}
      /// @brief Returns the differentiated or indexed operand.
      /// @returns The differentiated or indexed operand.
      const OperandType& getOperand() const
      {
        return m_gradient.getOperand();
      }
      /// @brief Evaluates the expression at the supplied physical or integration point.
      /// @param point Point at which the operation is evaluated.
      /// @returns Value of the expression at the supplied evaluation point.
      template <class Point>
      RangeType getValue(const Point& point) const
      {
        auto gradientExpression = m_gradient;
        gradientExpression.traceOf(this->getTraceDomain());
        const auto gradient = gradientExpression.getValue(point);
        return gradient;
      }
      /// @brief Returns the polynomial order when it is known.
      /// @param poly Mesh entity used by this operation.
      /// @returns Polynomial order on the entity, or an empty optional when no order is available.
      Optional<size_t> getOrder(const Geometry::Polytope& poly) const noexcept
      {
        return m_gradient.getOrder(poly);
      }
      Jacobian* copy() const noexcept override
      {
        return new Jacobian(*this);
      }

    private:
      Grad<OperandType> m_gradient;
  };

  /**
   * @ingroup JacobianSpecializations
   * @brief Matrix-basis jacobian, @f$ J_{ijk}=\partial_k A_{ij} @f$.
   */
  template <class Derived, class FES, ShapeFunctionSpaceType Space>
    requires FormLanguage::IsMatrixRange<
      typename FormLanguage::Traits<FES>::RangeType>::Value
  class Jacobian<ShapeFunction<Derived, FES, Space>> final
    : public ShapeFunctionBase<Jacobian<ShapeFunction<Derived, FES, Space>>, FES, Space>
  {
    public:
      /// @brief CRTP or finite element base class.
      using Parent = ShapeFunctionBase<Jacobian, FES, Space>;
      /// @brief Cloned or referenced expression operand type.
      using OperandType = ShapeFunction<Derived, FES, Space>;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialTensor<ScalarType>;
      /// @brief Constructs a rank-three matrix Jacobian.
      /// @param operand Operand expression.
      Jacobian(const OperandType& operand)
        : Parent(operand.getFiniteElementSpace()),
          m_gradient(operand)
      {}
      /// @brief Constructs a rank-three matrix Jacobian.
      /// @param other Object to copy from.
      Jacobian(const Jacobian& other)
        : Parent(other),
          m_gradient(other.m_gradient)
      {}
      /// @brief Constructs a rank-three matrix Jacobian.
      /// @param other Object to move from.
      Jacobian(Jacobian&& other)
        : Parent(std::move(other)),
          m_gradient(std::move(other.m_gradient))
      {}
      /// @brief Returns the differentiated or indexed operand.
      /// @returns The differentiated or indexed operand.
      const OperandType& getOperand() const
      {
        return m_gradient.getOperand();
      }
      /// @brief Returns the leaf shape function used for assembly.
      /// @returns The leaf shape function used for assembly.
      const auto& getLeaf() const
      {
        return m_gradient.getLeaf();
      }
      /// @brief Returns the local basis count for the selected polytope.
      /// @param poly Mesh entity used by this operation.
      /// @returns Number of local basis functions on the selected entity.
      size_t getDOFs(const Geometry::Polytope& poly) const
      {
        return m_gradient.getDOFs(poly);
      }
      /// @brief Returns the currently bound integration point.
      /// @returns The currently bound integration point.
      const IntegrationPoint& getIntegrationPoint() const
      {
        return m_gradient.getIntegrationPoint();
      }
      /// @brief Binds the integration point and prepares local basis values.
      /// @param point Point at which the operation is evaluated.
      /// @returns Reference to this object after the operation.
      Jacobian& setIntegrationPoint(const IntegrationPoint& point)
      {
        m_gradient.setIntegrationPoint(point);
        return *this;
      }
      /// @brief Returns a basis value at the bound integration point.
      /// @param local Index in the local numbering.
      /// @returns Value of the selected local basis function at the evaluation point.
      RangeType getBasis(size_t local) const
      {
        const auto gradient = m_gradient.getBasis(local);
        return gradient;
      }
      /// @brief Returns the polynomial order when it is known.
      /// @param poly Mesh entity used by this operation.
      /// @returns Polynomial order on the entity, or an empty optional when no order is available.
      Optional<size_t> getOrder(const Geometry::Polytope& poly) const noexcept
      {
        return m_gradient.getOrder(poly);
      }
      Jacobian* copy() const noexcept override
      {
        return new Jacobian(*this);
      }

    private:
      Grad<OperandType> m_gradient;
  };
  /// @brief Deduces the matrix space or coefficient type from constructor arguments.
  template <class FES, class Data>
    requires FormLanguage::IsMatrixRange<
               typename FormLanguage::Traits<FES>::RangeType>::Value
  Jacobian(const GridFunction<FES, Data>&) -> Jacobian<GridFunction<FES, Data>>;
  /// @brief Deduces the matrix space or coefficient type from constructor arguments.
  template <class Derived, class FES, ShapeFunctionSpaceType Space>
    requires FormLanguage::IsMatrixRange<
               typename FormLanguage::Traits<FES>::RangeType>::Value
  Jacobian(const ShapeFunction<Derived, FES, Space>&)
    -> Jacobian<ShapeFunction<Derived, FES, Space>>;
}

#endif
