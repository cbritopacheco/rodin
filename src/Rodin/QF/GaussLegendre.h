/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_QF_GAUSSLEGENDRE_H
#define RODIN_VARIATIONAL_QF_GAUSSLEGENDRE_H

/**
 * @file
 * @brief Defines the GaussLegendre quadrature formula.
 */

#include <vector>
#include <cassert>

#include "QuadratureFormula.h"


namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   * @brief Gauss-Legendre quadrature formula on reference polytopes.
   *
   * The Gauss-Legendre quadrature rule provides optimal integration accuracy
   * for polynomials. In one dimension, an @f$ n @f$-point Gauss-Legendre rule
   * is exact for polynomials of degree up to @f$ 2n-1 @f$. The quadrature
   * points are the zeros of the Legendre polynomial @f$ P_n(x) @f$, and the
   * weights are chosen to maximize the degree of exactness.
   *
   * ## Mathematical Foundation
   *
   * For the unit interval @f$ [0,1] @f$, the rule approximates:
   * @f[
   *   \int_0^1 f(x) \, dx \approx \sum_{i=1}^n w_i f(x_i)
   * @f]
   * where the points @f$ x_i @f$ and weights @f$ w_i @f$ are uniquely
   * determined by the requirement that the formula is exact for all
   * polynomials of degree at most @f$ 2n-1 @f$.
   *
   * ## Higher Dimensions
   *
   * For tensor-product geometries (quadrilaterals, wedges), the multi-dimensional
   * rule is constructed via tensor products of 1D Gauss-Legendre rules. For
   * simplicial geometries (triangles, tetrahedra), the Duffy transformation
   * maps tensor-product rules to the simplex:
   * - Triangle: @f$ (r,s) = (u, (1-u)v) @f$ with Jacobian @f$ (1-u) @f$
   * - Tetrahedron: @f$ (r,s,t) = (u, (1-u)v, (1-u)(1-v)w) @f$ with Jacobian
   *   @f$ (1-u)^2(1-v) @f$
   *
   * ## Supported Geometries
   * - Point
   * - Segment
   * - Triangle
   * - Quadrilateral
   * - Tetrahedron
   * - Wedge
   * - Hexahedron
   */
  class GaussLegendre final : public QuadratureFormulaBase
  {
    public:
      /// Parent class type
      using Parent = QuadratureFormulaBase;

      /**
       * @brief Constructs Gauss-Legendre quadrature with uniform order.
       * @param g Geometry type
       * @param order Number of Gauss-Legendre points per 1D direction
       *
       * Creates a quadrature rule with the same order in all coordinate
       * directions. Requires @f$ \text{order} \geq 1 @f$.
       */
      GaussLegendre(Geometry::Polytope::Type g, size_t order)
        : Parent(g), m_nx(order), m_ny(order), m_nz(order)
      {
        assert(order >= 1);
        build();
      }

      /**
       * @brief Constructs Gauss-Legendre quadrature with separate 2D orders.
       * @param g Geometry type
       * @param nx Order in first direction
       * @param ny Order in second direction
       *
       * Allows different orders in each coordinate direction for 2D geometries.
       * Requires @f$ nx, ny \geq 1 @f$.
       */
      GaussLegendre(Geometry::Polytope::Type g, size_t nx, size_t ny)
        : Parent(g), m_nx(nx), m_ny(ny), m_nz(ny)
      {
        assert(nx >= 1 && ny >= 1);
        build();
      }

      /**
       * @brief Constructs Gauss-Legendre quadrature with separate 3D orders.
       * @param g Geometry type
       * @param nu Order in first direction
       * @param nv Order in second direction
       * @param nw Order in third direction
       *
       * Allows different orders in each coordinate direction for 3D geometries.
       * Requires @f$ nu, nv, nw \geq 1 @f$.
       */
      GaussLegendre(Geometry::Polytope::Type g, size_t nu, size_t nv, size_t nw)
        : Parent(g), m_nx(nu), m_ny(nv), m_nz(nw)
      {
        assert(nu >= 1 && nv >= 1 && nw >= 1);
        build();
      }

      /**
       * @brief Constructs Gauss-Legendre quadrature with default order 2.
       * @param g Geometry type
       *
       * Uses 2 points per direction by default, providing exactness for
       * polynomials of degree up to 3 in 1D.
       */
      GaussLegendre(Geometry::Polytope::Type g)
        : Parent(g), m_nx(2), m_ny(2), m_nz(2)
      {
        build();
      }

      /**
       * @brief Gets the number of quadrature points.
       * @return Total number of quadrature points
       */
      size_t getSize() const override
      {
        return m_points.size();
      }

      /**
       * @brief Gets the coordinates of the i-th quadrature point.
       * @param i Index of the quadrature point
       * @return Reference to the point coordinates in reference space
       */
      const Math::SpatialVector<Real>& getPoint(size_t i) const override
      {
        assert(i < m_points.size());
        return m_points[i];
      }

      /**
       * @brief Gets the weight of the i-th quadrature point.
       * @param i Index of the quadrature point
       * @return Weight @f$ w_i @f$ associated with the point
       */
      Real getWeight(size_t i) const override
      {
        assert(i < static_cast<size_t>(m_weights.size()));
        return m_weights.coeff(i);
      }

      /**
       * @brief Creates a copy of this quadrature formula.
       * @return Pointer to a new GaussLegendre instance
       */
      GaussLegendre* copy() const noexcept override
      {
        return new GaussLegendre(*this);
      }

    private:
      /**
       * @brief Computes 1D Gauss-Legendre nodes and weights on [0,1].
       * @param n Number of quadrature points
       * @param x Output vector for quadrature point coordinates
       * @param w Output vector for quadrature weights
       * @param maxIt Maximum Newton iterations for finding nodes (default: 100)
       * @param tol Convergence tolerance for Newton iteration (default: 1e-14)
       *
       * Computes the @f$ n @f$ Gauss-Legendre quadrature nodes and weights
       * on the unit interval @f$ [0,1] @f$. Nodes are the zeros of the
       * Legendre polynomial @f$ P_n(x) @f$, found using Newton's method.
       */
      static void gl_1d_unit(size_t n, std::vector<Real>& x, std::vector<Real>& w, size_t maxIt = 100, Real tol = 1e-14);

      /**
       * @brief Builds the quadrature rule for the selected geometry.
       */
      void build();

      /**
       * @brief Builds quadrature for a point (0D).
       */
      void build_point();

      /**
       * @brief Builds 1D Gauss-Legendre quadrature on a segment.
       * @param n Number of quadrature points
       */
      void build_segment(size_t n);

      /**
       * @brief Builds Gauss-Legendre quadrature on a quadrilateral.
       * @param nx Number of points in first direction
       * @param ny Number of points in second direction
       */
      void build_quad(size_t nx, size_t ny);

      /**
       * @brief Builds Gauss-Legendre quadrature on a triangle.
       * @param nu Number of points in first parametric direction
       * @param nv Number of points in second parametric direction
       */
      void build_tri(size_t nu, size_t nv);

      /**
       * @brief Builds Gauss-Legendre quadrature on a tetrahedron.
       * @param nu Number of points in first parametric direction
       * @param nv Number of points in second parametric direction
       * @param nw Number of points in third parametric direction
       */
      void build_tet(size_t nu, size_t nv, size_t nw);

      /**
       * @brief Builds Gauss-Legendre quadrature on a wedge (prism).
       * @param ntri Number of points for the triangular cross-section
       * @param nz Number of points in the extrusion direction
       */
      void build_wedge(size_t ntri, size_t nz);

      /**
       * @brief Builds Gauss-Legendre quadrature on a hexahedron (cube).
       * @param nx Number of points in x-direction
       * @param ny Number of points in y-direction
       * @param nz Number of points in z-direction
       */
      void build_hex(size_t nx, size_t ny, size_t nz);

      size_t m_nx { 2 }, m_ny { 2 }, m_nz { 2 }; ///< Number of points per direction
      std::vector<Math::SpatialVector<Real>> m_points; ///< Quadrature point coordinates
      Math::Vector<Real> m_weights; ///< Quadrature weights
      size_t m_maxIt; ///< Maximum Newton iterations
      Real m_tol; ///< Newton convergence tolerance
  };
}

#endif
