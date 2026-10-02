/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file Grad.h
 * @brief Gradient operator for scalar functions and shape functions.
 *
 * This file defines the Grad class, which computes the gradient (spatial 
 * derivative) of scalar functions in variational formulations. The gradient
 * is a fundamental differential operator in finite element analysis.
 */
#ifndef RODIN_VARIATIONAL_GRAD_H
#define RODIN_VARIATIONAL_GRAD_H

#include "ForwardDecls.h"

#include "Rodin/Math/SpatialVector.h"
#include "Rodin/Math/SpatialTensor.h"

#include "VectorFunction.h"
#include "IntegrationPoint.h"

namespace Rodin::FormLanguage
{
  /// @brief Type traits for @c Grad over a grid function: exposes the finite element
  /// space, the operand type and the range type.
  template <class FES, class Data>
  struct Traits<Variational::Grad<Variational::GridFunction<FES, Data>>>
  {
      /// @brief Finite element space type.
      using FESType = FES;

      /// @brief Operand type.
      using OperandType = Variational::GridFunction<FESType, Data>;

      /// @brief Range (evaluation value) type.
      using RangeType =
        Math::SpatialVector<typename FormLanguage::Traits<FESType>::ScalarType>;
  };

  /// @brief Type traits for @c Grad over a shape function: exposes the finite element
  /// space, the shape function space, the operand type and the range type.
  template <class NestedDerived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Grad<Variational::ShapeFunction<NestedDerived, FES, Space>>>
  {
      /// @brief Finite element space type.
      using FESType = FES;
      /// @brief Shape function space the expression belongs to, trial or test.
      static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;

      /// @brief Operand type.
      using OperandType = Variational::ShapeFunction<NestedDerived, FESType, SpaceType>;

      /// @brief Range (evaluation value) type.
      using RangeType =
        Math::SpatialVector<typename FormLanguage::Traits<FESType>::ScalarType>;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup GradSpecializations Grad Template Specializations
   * @brief Template specializations of the Grad class.
   * @see @ref Grad
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref Grad "Grad<GridFunction<MatrixFES, Data>>" | Matrix-valued grid-function operator for all supported spaces and backends. |
   * | @ref Grad "Grad<ShapeFunction<Derived, MatrixFES, Space>>" | Matrix shape-function operator with the scalar family's geometry and trace semantics. |
   * | @ref GradBase "GradBase<GridFunction<FES, Data>, Derived>" | Generic gradient base for scalar grid functions. |
   * | @ref Grad "Grad<H1<K, Scalar, Mesh>, GridFunction<H1<K, Scalar, Mesh>, Data>>" | Gradient of an H1 grid function. |
   * | @ref Grad "Grad<H1<K, Scalar, Mesh>, ShapeFunction<NestedDerived, H1<K, Scalar, Mesh>, Space>>" | Gradient of an H1 shape-function expression. |
   * | @ref Grad "Grad<P1<Range, Mesh>, GridFunction<P1<Range, Mesh>, Data>>" | Gradient of a P1 grid function. |
   * | @ref Grad "Grad<P1<Range, Mesh>, ShapeFunction<NestedDerived, P1<Range, Mesh>, Space>>" | Gradient of a P1 shape-function expression. |
   * | @ref Grad "Grad<P0<Range, Mesh>, GridFunction<P0<Range, Mesh>, Data>>" | Gradient of a P0 grid function. |
   * | @ref Grad "Grad<P0<Range, Mesh>, ShapeFunction<NestedDerived, P0<Range, Mesh>, Space>>" | Gradient of a P0 shape-function expression. |
   */

  /**
   * @brief Base class for gradient operator implementations.
   *
   * GradBase provides the foundation for computing gradients of different
   * function types (grid functions, shape functions, etc.).
   *
   * @tparam Operand Type of the function being differentiated
   * @tparam Derived Derived class (CRTP pattern)
   */
  template <class Operand, class Derived>
  class GradBase;

  /**
   * @ingroup GradSpecializations
   * @brief Gradient operator for grid functions.
   *
   * Computes the spatial gradient of a scalar grid function:
   * @f[
   *   \nabla u = \left(\frac{\partial u}{\partial x_1}, \frac{\partial u}{\partial x_2}, \ldots, \frac{\partial u}{\partial x_d}\right)^T
   * @f]
   * where @f$ u @f$ is a scalar function and @f$ d @f$ is the spatial dimension.
   *
   * ## Mathematical Foundation
   * For a scalar function @f$ u : \Omega \subset \mathbb{R}^d \to \mathbb{R} @f$,
   * the gradient is a vector field:
   * @f[
   *   \nabla u : \Omega \to \mathbb{R}^d
   * @f]
   *
   * In the finite element context, for @f$ u_h = \sum_i u_i \phi_i @f$:
   * @f[
   *   \nabla u_h = \sum_i u_i \nabla \phi_i
   * @f]
   *
   * ## Usage Example
   * ```cpp
   * // In a variational formulation
   * auto stiffness = Integral(Grad(u), Grad(v));  // Laplacian term
   * ```
   *
   * @tparam FES Finite element space type
   * @tparam Data Data storage type
   * @tparam Derived Derived class for CRTP
   */
  template <class FES, class Data, class Derived>
  class GradBase<GridFunction<FES, Data>, Derived>
    : public VectorFunctionBase<
        typename FormLanguage::Traits<FES>::ScalarType, GradBase<GridFunction<FES, Data>, Derived>>
  {
    public:
      /// @brief Finite element space type
      using FESType = FES;

      /// @brief Type of the output range
      using RangeType = Math::SpatialVector<typename FormLanguage::Traits<FESType>::ScalarType>;

      /// @brief Scalar type for computations
      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      /// @brief Spatial vector type
      using SpatialVectorType = Math::SpatialVector<ScalarType>;

      /// @brief Type of the operand (grid function)
      using OperandType = GridFunction<FESType, Data>;

      /// @brief Parent class type
      using Parent = VectorFunctionBase<ScalarType, GradBase<OperandType, Derived>>;

      /**
       * @brief Constructs the gradient operator for a grid function.
       * @param[in] u Scalar grid function to differentiate
       *
       * @pre The grid function must be scalar-valued (vector dimension = 1)
       */
      GradBase(const OperandType& u)
        : m_u(u)
      {
        assert(u.getFiniteElementSpace().getVectorDimension() == 1);
      }

      /**
       * @brief Copy constructor
       */
      GradBase(const GradBase& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      /**
       * @brief Move constructor
       */
      GradBase(GradBase&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u))
      {}

      /**
       * @brief Gets the spatial dimension of the gradient.
       * @return Dimension of the spatial domain
       *
       * Returns the number of spatial dimensions @f$ d @f$ where the gradient
       * @f$ \nabla u \in \mathbb{R}^d @f$.
       */
      constexpr
      size_t getDimension() const
      {
        return m_u.get().getFiniteElementSpace().getMesh().getSpaceDimension();
      }

    public:
      /**
       * @brief Evaluates the gradient at a Point.
       *
       * Resolves mesh ownership and dispatches to the derived class's
       * @c interpolate. Falls back to inclusion / submesh restriction
       * when the polytope's mesh is not the FES mesh.
       */
      SpatialVectorType getValue(const Geometry::Point& p) const
      {
        const auto& fes = getOperand().getFiniteElementSpace();
        const auto& fesMesh = fes.getMesh();

        SpatialVectorType value;
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
       * @brief Evaluates the gradient at an IntegrationPoint.
       *
       * If the polytope is owned by the FES mesh, dispatches to
       * @c interpolate(out, ip). Otherwise falls back to inclusion / submesh
       * restriction.
       */
      SpatialVectorType getValue(const IntegrationPoint& ip) const
      {
        const auto& p = ip.getPoint();
        const auto& fes = getOperand().getFiniteElementSpace();
        const auto& fesMesh = fes.getMesh();

        if (fesMesh.isLocalPoint(p))
        {
          SpatialVectorType value;
          this->interpolate(value, ip);
          return value;
        }

        SpatialVectorType value;
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
       * @brief Interpolates the gradient at a point (to be overridden in derived class).
       * @param[out] out Output vector for gradient result
       * @param[in] p Point at which to interpolate
       *
       * This virtual function is overridden in derived classes (e.g., P1::Grad)
       * to provide finite element-specific gradient interpolation.
       */
      constexpr
      void interpolate(SpatialVectorType& out, const Geometry::Point& p) const
      {
        static_cast<const Derived&>(*this).interpolate(out, p);
      }

      /// @brief Interpolates at an integration point.
      constexpr
      void interpolate(SpatialVectorType& out, const IntegrationPoint& ip) const
      {
        if constexpr (requires (const Derived& f, SpatialVectorType& r, const IntegrationPoint& q) { f.interpolate(r, q); })
          static_cast<const Derived&>(*this).interpolate(out, ip);
        else
          static_cast<const Derived&>(*this).interpolate(out, ip.getPoint());
      }

      /**
       * @brief Gets the operand grid function.
       * @return Reference to the grid function being differentiated
       */
      constexpr
      const OperandType& getOperand() const
      {
        return m_u.get();
      }

      /// @brief Returns the polynomial order used on a mesh entity.
      constexpr
      Optional<size_t> getOrder(const Geometry::Polytope& polytope) const noexcept
      {
        return static_cast<const Derived&>(*this).getOrder(polytope);
      }

      /**
       * @brief Copy function to be overriden in Derived type.
       */
      GradBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
  };

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Grad of a GridFunction
   */
  template <class FES, class Data>
  Grad(const GridFunction<FES, Data>&) -> Grad<GridFunction<FES, Data>>;

  /**
   * @ingroup RodinCTAD
   * @brief CTAD for Grad of a ShapeFunction
   */
  template <class NestedDerived, class FES, ShapeFunctionSpaceType Space>
  Grad(const ShapeFunction<NestedDerived, FES, Space>&)
    -> Grad<ShapeFunction<NestedDerived, FES, Space>>;
}

namespace Rodin::FormLanguage
{
  template <class FES, class Data>
    requires IsMatrixRange<typename Traits<FES>::RangeType>::Value
  struct Traits<Variational::Grad<Variational::GridFunction<FES, Data>>>
  {
      /// @brief Finite element space type.
      using FESType = FES;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename Traits<FES>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialTensor<ScalarType>;
      /// @brief Cloned or referenced expression operand type.
      using OperandType = Variational::GridFunction<FES, Data>;
  };
  /// @brief Type traits for the matrix or tensor expression specialization.
  template <class Derived, class FES, Variational::ShapeFunctionSpaceType Space>
    requires IsMatrixRange<typename Traits<FES>::RangeType>::Value
  struct Traits<Variational::Grad<Variational::ShapeFunction<Derived, FES, Space>>>
  {
      /// @brief Finite element space type.
      using FESType = FES;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename Traits<FES>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialTensor<ScalarType>;
      /// @brief Cloned or referenced expression operand type.
      using OperandType = Variational::ShapeFunction<Derived, FES, Space>;
      /// @brief Trial or test shape-function space.
      static constexpr auto SpaceType = Space;
  };
}
namespace Rodin::Variational
{
  /**
   * @ingroup GradSpecializations
   * @brief Gradient of a matrix solution, @f$ G_{ijk}=\partial_k A_{ij} @f$.
   * Uses the space's basis and physical Jacobian for every geometry/backend.
   */
  template <class FES, class Data>
    requires FormLanguage::IsMatrixRange<
      typename FormLanguage::Traits<FES>::RangeType>::Value
  class Grad<GridFunction<FES, Data>> final
    : public FunctionBase<Grad<GridFunction<FES, Data>>>
  {
    public:
      /// @brief CRTP or finite element base class.
      using Parent = FunctionBase<Grad>;
      /// @brief Cloned or referenced expression operand type.
      using OperandType = GridFunction<FES, Data>;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialTensor<ScalarType>;
      /// @brief Constructs the matrix gradient with the derivative axis last.
      explicit Grad(const OperandType& operand)
        : m_operand(operand)
      {}
      /// @brief Constructs the matrix gradient with the derivative axis last.
      Grad(const Grad& other)
        : Parent(other),
          m_operand(other.m_operand)
      {}
      /// @brief Constructs the matrix gradient with the derivative axis last.
      Grad(Grad&& other)
        : Parent(std::move(other)),
          m_operand(other.m_operand)
      {}
      /// @brief Returns the differentiated or indexed operand.
      const OperandType& getOperand() const
      {
        return m_operand.get();
      }
      /// @brief Evaluates the expression at the supplied physical or integration point.
      RangeType getValue(const Geometry::Point& point) const
      {
        const auto& fes = getOperand().getFiniteElementSpace();
        const auto p = fes.getDerivativePoint(point, this->getTraceDomain());
        const auto& poly = p.getPolytope();
        const auto& dofs = fes.getDOFs(poly.getDimension(), poly.getIndex());
        RangeType value(
          fes.getRows(), fes.getColumns(), fes.getMesh().getSpaceDimension());
        value.setZero();
        for (size_t a = 0; a < static_cast<size_t>(dofs.size()); ++a)
          value += getOperand()[dofs[a]] * fes.getGradientBasis(a, p);
        return value;
      }
      /// @brief Evaluates the expression at the supplied physical or integration point.
      RangeType getValue(const IntegrationPoint& point) const
      {
        return getValue(point.getPoint());
      }
      /// @brief Returns the polynomial order when it is known.
      Optional<size_t> getOrder(const Geometry::Polytope& poly) const noexcept
      {
        const auto order = getOperand().getOrder(poly);
        return order ? Optional<size_t>(*order ? *order - 1 : 0) : std::nullopt;
      }
      Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }

    private:
      std::reference_wrapper<const OperandType> m_operand;
  };

  /**
   * @ingroup GradSpecializations
   * @brief Matrix basis gradients as rank-three tensors.
   */
  template <class Derived, class FES, ShapeFunctionSpaceType Space>
    requires FormLanguage::IsMatrixRange<
      typename FormLanguage::Traits<FES>::RangeType>::Value
  class Grad<ShapeFunction<Derived, FES, Space>> final
    : public ShapeFunctionBase<Grad<ShapeFunction<Derived, FES, Space>>, FES, Space>
  {
    public:
      /// @brief CRTP or finite element base class.
      using Parent = ShapeFunctionBase<Grad, FES, Space>;
      /// @brief Cloned or referenced expression operand type.
      using OperandType = ShapeFunction<Derived, FES, Space>;
      /// @brief Scalar type of matrix or tensor entries.
      using ScalarType = typename FormLanguage::Traits<FES>::ScalarType;
      /// @brief Evaluated matrix, tensor, or scalar range type.
      using RangeType = Math::SpatialTensor<ScalarType>;
      /// @brief Constructs the matrix gradient with the derivative axis last.
      explicit Grad(const OperandType& operand)
        : Parent(operand.getFiniteElementSpace()),
          m_operand(operand.copy())
      {}
      /// @brief Constructs the matrix gradient with the derivative axis last.
      Grad(const Grad& other)
        : Parent(other),
          m_operand(other.m_operand->copy())
      {}
      /// @brief Constructs the matrix gradient with the derivative axis last.
      Grad(Grad&& other)
        : Parent(std::move(other)),
          m_operand(std::move(other.m_operand))
      {}
      /// @brief Returns the differentiated or indexed operand.
      const OperandType& getOperand() const
      {
        return *m_operand;
      }
      /// @brief Returns the leaf shape function used for assembly.
      const auto& getLeaf() const
      {
        return getOperand().getLeaf();
      }
      /// @brief Returns the local basis count for the selected polytope.
      size_t getDOFs(const Geometry::Polytope& poly) const
      {
        return getOperand().getDOFs(poly);
      }
      /// @brief Returns the currently bound integration point.
      const IntegrationPoint& getIntegrationPoint() const
      {
        assert(m_point);
        return *m_point;
      }
      /// @brief Binds the integration point and prepares local basis values.
      Grad& setIntegrationPoint(const IntegrationPoint& point)
      {
        m_point = &point;
        const auto& fes = this->getFiniteElementSpace();
        const auto& p = point.getPoint();
        const auto& poly = p.getPolytope();
        const size_t count =
          fes.getFiniteElement(poly.getDimension(), poly.getIndex()).getCount();
        m_basis.resize(count);
        for (size_t a = 0; a < count; ++a)
          m_basis[a] = fes.getGradientBasis(a, p);
        return *this;
      }
      /// @brief Returns a basis value at the bound integration point.
      const RangeType& getBasis(size_t local) const
      {
        return m_basis.at(local);
      }
      /// @brief Returns the polynomial order when it is known.
      Optional<size_t> getOrder(const Geometry::Polytope& poly) const noexcept
      {
        const auto order = getOperand().getOrder(poly);
        return order ? Optional<size_t>(*order ? *order - 1 : 0) : std::nullopt;
      }
      Grad* copy() const noexcept override
      {
        return new Grad(*this);
      }

    private:
      std::unique_ptr<OperandType> m_operand;
      const IntegrationPoint* m_point = nullptr;
      std::vector<RangeType> m_basis;
  };
}

#endif
