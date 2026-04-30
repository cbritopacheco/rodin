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

#include "VectorFunction.h"
#include "IntegrationPoint.h"

namespace Rodin::FormLanguage
{
  template <class FES, class Data>
  struct Traits<Variational::Grad<Variational::GridFunction<FES, Data>>>
  {
    using FESType = FES;

    using OperandType = Variational::GridFunction<FESType, Data>;

    using RangeType = Math::SpatialVector<typename FormLanguage::Traits<FESType>::ScalarType>;
  };

  template <class NestedDerived, class FES, Variational::ShapeFunctionSpaceType Space>
  struct Traits<
    Variational::Grad<Variational::ShapeFunction<NestedDerived, FES, Space>>>
  {
    using FESType = FES;
    static constexpr Variational::ShapeFunctionSpaceType SpaceType = Space;

    using OperandType = Variational::ShapeFunction<NestedDerived, FESType, SpaceType>;

    using RangeType = Math::SpatialVector<typename FormLanguage::Traits<FESType>::ScalarType>;
  };
}

namespace Rodin::Variational
{
  /**
   * @defgroup GradSpecializations Grad Template Specializations
   * @brief Template specializations of the Grad class.
   * @see Grad
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

    protected:
      /**
       * @brief Evaluation cache keyed by (mesh, polytope, quadrature point).
       *
       * Populated by getValue(IntegrationPoint) on a hit-miss path. Invalidated
       * whenever evaluation falls back to inclusion / submesh restriction.
       */
      struct Cache
      {
        struct Key
        {
          const void* mesh = nullptr;
          Geometry::Polytope::Type geom = Geometry::Polytope::Type::Point;
          size_t dim = 0;
          Index cell = 0;
          const QF::QuadratureFormulaBase* qf = nullptr;
          size_t qp = 0;
          bool valid = false;

          bool operator==(const Key& o) const noexcept
          {
            if (!valid || !o.valid)
              return false;
            return mesh == o.mesh
                && geom == o.geom
                && dim  == o.dim
                && cell == o.cell
                && qf   == o.qf
                && qp   == o.qp;
          }
        };

        Key key;
        SpatialVectorType value;
      };

      mutable Cache m_cache;

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

        m_cache.key.valid = false;
        if (fesMesh.isLocalPoint(p))
        {
          this->interpolate(m_cache.value, p);
        }
        else if (const auto inclusion = fesMesh.inclusion(p))
        {
          this->interpolate(m_cache.value, *inclusion);
        }
        else if (fesMesh.isSubMesh())
        {
          const auto& submesh = fesMesh.asSubMesh();
          const auto restriction = submesh.restriction(p);
          this->interpolate(m_cache.value, *restriction);
        }
        else
        {
          assert(false);
        }
        return m_cache.value;
      }

      /**
       * @brief Evaluates the gradient at an IntegrationPoint.
       *
       * If the polytope is owned by the FES mesh, performs a cache lookup
       * keyed by (mesh, polytope, qf, qp); on a miss it dispatches to
       * @c interpolate(out, ip) and caches the result. Falls back to
       * inclusion / submesh restriction otherwise.
       */
      SpatialVectorType getValue(const IntegrationPoint& ip) const
      {
        const auto& p = ip.getPoint();
        const auto& polytope = p.getPolytope();
        const auto& fes = getOperand().getFiniteElementSpace();
        const auto& fesMesh = fes.getMesh();

        if (fesMesh.isLocalPoint(p))
        {
          typename Cache::Key key;
          key.mesh = static_cast<const void*>(&fesMesh);
          key.geom = polytope.getGeometry();
          key.dim  = polytope.getDimension();
          key.cell = polytope.getIndex();
          key.qf   = &ip.getQuadratureFormula();
          key.qp   = ip.getIndex();
          key.valid = true;

          if (m_cache.key == key)
            return m_cache.value;

          m_cache.key = key;
          this->interpolate(m_cache.value, ip);
          return m_cache.value;
        }

        m_cache.key.valid = false;
        if (const auto inclusion = fesMesh.inclusion(p))
        {
          this->interpolate(m_cache.value, *inclusion);
        }
        else if (fesMesh.isSubMesh())
        {
          const auto& submesh = fesMesh.asSubMesh();
          const auto restriction = submesh.restriction(p);
          this->interpolate(m_cache.value, *restriction);
        }
        else
        {
          assert(false);
        }
        return m_cache.value;
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

#endif
