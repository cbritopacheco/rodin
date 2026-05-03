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
    * @see Div
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
      using FESType = FES;

      using ScalarType = typename FormLanguage::Traits<FESType>::ScalarType;

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
       */
      DivBase(const DivBase& other)
        : Parent(other),
          m_u(other.m_u)
      {}

      /**
       * @brief Move constructor
       */
      DivBase(DivBase&& other)
        : Parent(std::move(other)),
          m_u(std::move(other.m_u))
      {}

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
        ScalarType value = ScalarType(0);
      };

      mutable Cache m_cache;

    public:
      /**
       * @brief Evaluates the divergence at a Point.
       *
       * Resolves mesh ownership and dispatches to the derived class's
       * @c interpolate. Falls back to inclusion / submesh restriction
       * when the polytope's mesh is not the FES mesh.
       */
      ScalarType getValue(const Geometry::Point& p) const
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
       * @brief Evaluates the divergence at an IntegrationPoint.
       *
       * If the polytope is owned by the FES mesh, performs a cache lookup
       * keyed by (mesh, polytope, qf, qp); on a miss it dispatches to
       * @c interpolate(out, ip) and caches the result. Falls back to
       * inclusion / submesh restriction otherwise.
       */
      ScalarType getValue(const IntegrationPoint& ip) const
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

      constexpr
      void interpolate(ScalarType& out, const IntegrationPoint& ip) const
      {
        if constexpr (requires (const Derived& f, ScalarType& r, const IntegrationPoint& q) { f.interpolate(r, q); })
          static_cast<const Derived&>(*this).interpolate(out, ip);
        else
          static_cast<const Derived&>(*this).interpolate(out, ip.getPoint());
      }

      Optional<size_t> getOrder(const Geometry::Polytope& poly) const noexcept
      {
        return static_cast<const Derived&>(*this).getOrder(poly);
      }

      /**
       * @brief Copy function to be overriden in Derived type.
       */
      DivBase* copy() const noexcept override
      {
        return static_cast<const Derived&>(*this).copy();
      }

    private:
      std::reference_wrapper<const OperandType> m_u;
  };
}

#endif
