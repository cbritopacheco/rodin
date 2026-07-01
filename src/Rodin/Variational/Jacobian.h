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
 * auto F = Jacobian(u);  // Deformation gradient F = \nablau
 * ```
 */
#ifndef RODIN_VARIATIONAL_JACOBIAN_H
#define RODIN_VARIATIONAL_JACOBIAN_H

#include "ForwardDecls.h"
#include "Rodin/Math/SpatialMatrix.h"
#include "MatrixFunction.h"
#include "IntegrationPoint.h"

namespace Rodin::Variational
{
  /**
   * @defgroup JacobianSpecializations Jacobian Template Specializations
   * @brief Template specializations of the Jacobian class.
   * @see Jacobian
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
      using FESType = FES;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

      using RangeType = Math::SpatialMatrix<ScalarType>;

      using SpatialMatrixType = Math::SpatialMatrix<ScalarType>;

      using OperandType = GridFunction<FESType, Data>;

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

    protected:
      /**
       * @brief Evaluation cache keyed by (mesh, polytope, quadrature point).
       *
       * Populated on a hit-miss path by getValue(IntegrationPoint). The
       * matching cell-level key is invalidated whenever evaluation falls
       * back to inclusion/restriction.
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
        SpatialMatrixType value;
      };

      mutable Cache m_cache;

    public:
      /**
       * @brief Evaluates the Jacobian matrix at a Point.
       *
       * Resolves mesh ownership and dispatches to the derived class's
       * @c interpolate. Falls back to inclusion / submesh restriction
       * when the polytope's mesh is not the FES mesh.
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
       * If the polytope is owned by the FES mesh, performs a cache lookup
       * keyed by (mesh, polytope, qf, qp); on a miss it dispatches to
       * @c interpolate(out, ip) and caches the result. Falls back to
       * inclusion / submesh restriction when the polytope's mesh is not
       * the FES mesh.
       */
      SpatialMatrixType getValue(const IntegrationPoint& ip) const
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
          key.qf   = ip.getQuadratureFormula();
          key.qp   = ip.getIndex();
          key.valid = true;

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
       */
      constexpr
      void interpolate(SpatialMatrixType& out, const IntegrationPoint& ip) const
      {
        if constexpr (requires (const Derived& f, SpatialMatrixType& r, const IntegrationPoint& q) { f.interpolate(r, q); })
          static_cast<const Derived&>(*this).interpolate(out, ip);
        else
          static_cast<const Derived&>(*this).interpolate(out, ip.getPoint());
      }

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

#endif
