/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_QF_GAUSSLOBATTO_H
#define RODIN_VARIATIONAL_QF_GAUSSLOBATTO_H

#include <vector>
#include <cassert>
#include <cmath>

#include "QuadratureFormula.h"

namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   * @brief Gauss-Lobatto quadrature formula on reference polytopes.
   *
   * The Gauss-Lobatto quadrature rule is similar to Gauss-Legendre
   * quadrature but includes the endpoints of the integration interval.
   * In one dimension on the interval @f$ [0,1] @f$, the quadrature points
   * include both @f$ 0 @f$ and @f$ 1 @f$, with interior points being the
   * zeros of the derivative of the Legendre polynomial @f$ P'_{n-1}(x) @f$.
   *
   * ## Mathematical Foundation
   *
   * The 1D Gauss-Lobatto rule with @f$ n @f$ points is exact for polynomials
   * of degree up to @f$ 2n-3 @f$. The weights @f$ w_i @f$ and points
   * @f$ x_i @f$ satisfy:
   * @f[
   *   \int_0^1 f(x) \, dx \approx \sum_{i=1}^n w_i f(x_i)
   * @f]
   *
   * ## Higher Dimensions
   *
   * For tensor-product geometries (quadrilaterals, wedges), the rule is
   * constructed via tensor products of 1D Gauss-Lobatto rules. For simplicial
   * geometries (triangles, tetrahedra), the Duffy transformation is applied:
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
   *
   * @note Requires at least @f$ n \geq 2 @f$ points per direction for proper
   * Lobatto quadrature (to include both endpoints).
   */
  class GaussLobatto final : public QuadratureFormulaBase
  {
    public:
      /// Parent class type
      using Parent = QuadratureFormulaBase;

      /**
       * @brief Constructs Gauss-Lobatto quadrature with uniform order.
       * @param g Geometry type
       * @param order Number of Gauss-Lobatto points per 1D direction
       *
       * Creates a quadrature rule with the same order in all coordinate
       * directions. Requires @f$ \text{order} \geq 2 @f$.
       */
      GaussLobatto(Geometry::Polytope::Type g, size_t order)
        : Parent(g), m_nx(order), m_ny(order), m_nz(order)
      { assert(order >= 2); build(); }

      /**
       * @brief Constructs Gauss-Lobatto quadrature with separate 2D orders.
       * @param g Geometry type
       * @param nx Order in first direction
       * @param ny Order in second direction
       *
       * Allows different orders in each coordinate direction for 2D geometries.
       * Requires @f$ nx, ny \geq 2 @f$.
       */
      GaussLobatto(Geometry::Polytope::Type g, size_t nx, size_t ny)
        : Parent(g), m_nx(nx), m_ny(ny), m_nz(ny)
      { assert(nx>=2 && ny>=2); build(); }

      /**
       * @brief Constructs Gauss-Lobatto quadrature with separate 3D orders.
       * @param g Geometry type
       * @param nu Order in first direction
       * @param nv Order in second direction
       * @param nw Order in third direction
       *
       * Allows different orders in each coordinate direction for 3D geometries.
       * Requires @f$ nu, nv, nw \geq 2 @f$.
       */
      GaussLobatto(Geometry::Polytope::Type g, size_t nu, size_t nv, size_t nw)
        : Parent(g), m_nx(nu), m_ny(nv), m_nz(nw)
      { assert(nu>=2 && nv>=2 && nw>=2); build(); }

      /**
       * @brief Constructs Gauss-Lobatto quadrature with default order 2.
       * @param g Geometry type
       *
       * Uses 2 points per direction, which is the minimum for Lobatto quadrature.
       */
      GaussLobatto(Geometry::Polytope::Type g)
        : Parent(g), m_nx(2), m_ny(2), m_nz(2)
      { build(); }

      /**
       * @brief Gets the number of quadrature points.
       * @return Total number of quadrature points
       */
      size_t getSize() const override { return m_points.size(); }

      /**
       * @brief Gets the coordinates of the i-th quadrature point.
       * @param i Index of the quadrature point
       * @return Reference to the point coordinates in reference space
       */
      const Math::SpatialVector<Real>& getPoint(size_t i) const override
      { return m_points[i]; }

      /**
       * @brief Gets the weight of the i-th quadrature point.
       * @param i Index of the quadrature point
       * @return Weight @f$ w_i @f$ associated with the point
       */
      Real getWeight(size_t i) const override
      { return m_weights.coeff(i); }

      /**
       * @brief Creates a copy of this quadrature formula.
       * @return Pointer to a new GaussLobatto instance
       */
      GaussLobatto* copy() const noexcept override
      { return new GaussLobatto(*this); }

    private:
      /**
       * @brief Computes 1D Gauss-Lobatto nodes and weights on [0,1].
       * @param n Number of quadrature points (must be at least 2)
       * @param x01 Output vector for quadrature point coordinates
       * @param w01 Output vector for quadrature weights
       * @param maxIt Maximum Newton iterations for finding interior nodes (default: 100)
       * @param tol Convergence tolerance for Newton iteration (default: 1e-14)
       *
       * Computes the @f$ n @f$ Gauss-Lobatto quadrature nodes and weights
       * on the unit interval @f$ [0,1] @f$. The endpoints @f$ x_1 = 0 @f$
       * and @f$ x_n = 1 @f$ are included. Interior nodes @f$ x_i @f$ for
       * @f$ i = 2, \ldots, n-1 @f$ are found as zeros of the derivative
       * @f$ P'_{n-1}(x) @f$ of the Legendre polynomial using Newton's method.
       *
       * The weights are computed using:
       * @f[
       *   w_i = \frac{2}{n(n-1) [P_{n-1}(x_i)]^2}
       * @f]
       */
      static void gll_1d_unit(size_t n, std::vector<Real>& x01, std::vector<Real>& w01,
                              size_t maxIt = 100, Real tol = 1e-14)
      {
        assert(n >= 2);
        // On [-1,1], endpoints ±1, interior are roots of P'_{n-1}(x).
        std::vector<Real> xi(n), wi(n);

        xi[0] = -1.0; xi[n-1] = 1.0;
        // endpoints weights on [-1,1]
        wi[0] = wi[n-1] = 2.0 / (static_cast<Real>(n)*(static_cast<Real>(n)-1.0));

        if (n > 2)
        {
          const size_t m = n - 2; // number of interior nodes
          for (size_t k = 0; k < m; ++k)
          {
            // good initial guess for roots of P'_{n-1} (Chebyshev-like)
            const Real kk = static_cast<Real>(k) + 1.0;
            Real xk = -std::cos(M_PI * kk / static_cast<Real>(n-1));

            // Newton on P'_{n-1}(x) = 0
            for (size_t it = 0; it < maxIt; ++it)
            {
              // Evaluate P_{n-1}(xk) and P_{n-2}(xk) by recurrence
              const size_t mdeg = n - 1;
              Real p0 = 1.0;
              Real p1 = xk;
              for (size_t j = 2; j <= mdeg; ++j)
              {
                const Real jj = static_cast<Real>(j);
                const Real pj = ((2*jj-1)*xk*p1 - (jj-1)*p0) / jj;
                p0 = p1; p1 = pj;
              }
              const Real Pn1   = p1; // P_{n-1}
              const Real Pn2   = p0; // P_{n-2}
              // P'_{n-1}(x) via (n-1)(x P_{n-1}-P_{n-2})/(x^2-1)
              const Real dPn1  = (static_cast<Real>(mdeg) * (xk*Pn1 - Pn2)) / (xk*xk - 1.0);

              // Need derivative of P'_{n-1}: use identity for second derivative
              // P''_{n-1}(x) can be derived; practical route: differentiate recurrence numerically
              // Compute also P'_{n-2}(x) to avoid cancellation
              // Here we use a stable relation:
              // d/dx P_j satisfies: P'_j = j/(x^2-1) (x P_j - P_{j-1})
              // Differentiate again (closed form). For robustness, use a small step secant on dPn1.
              Real xk_eps = xk * (1.0 + 1e-8) + ((xk==0)?1e-8:0);
              // evaluate dP at xk_eps
              {
                Real p0e = 1.0, p1e = xk_eps;
                for (size_t j = 2; j <= mdeg; ++j)
                {
                  const Real jj = static_cast<Real>(j);
                  const Real pj = ((2*jj-1)*xk_eps*p1e - (jj-1)*p0e) / jj;
                  p0e = p1e; p1e = pj;
                }
                const Real Pn1e  = p1e;
                const Real Pn2e  = p0e;
                const Real dPn1e = (static_cast<Real>(mdeg) * (xk_eps*Pn1e - Pn2e)) / (xk_eps*xk_eps - 1.0);
                const Real d2 = (dPn1e - dPn1) / (xk_eps - xk);
                const Real dx = dPn1 / d2;
                xk -= dx;
                if (std::abs(dx) < tol) break;
              }
            }

            // Final P_{n-1} value at converged xk for weight
            const size_t mdeg2 = n - 1;
            Real p0 = 1.0, p1 = xk;
            for (size_t j = 2; j <= mdeg2; ++j)
            {
              const Real jj = static_cast<Real>(j);
              const Real pj = ((2*jj-1)*xk*p1 - (jj-1)*p0) / jj;
              p0 = p1; p1 = pj;
            }
            const Real Pn1 = p1;

            xi[k+1] = xk;
            wi[k+1] = 2.0 / (static_cast<Real>(n)*(static_cast<Real>(n)-1.0) * Pn1*Pn1);
          }
        }

        // map to [0,1]
        x01.resize(n);
        w01.resize(n);
        for (size_t i=0;i<n;++i) {
          x01[i] = 0.5*(xi[i] + 1.0);
          w01[i] = 0.5*wi[i];
        }
      }

      /**
       * @brief Builds the quadrature rule for the selected geometry.
       *
       * Dispatches to the appropriate geometry-specific builder based on
       * the polytope type.
       */
      void build()
      {
        switch (getGeometry())
        {
          case Geometry::Polytope::Type::Point:         build_point(); break;
          case Geometry::Polytope::Type::Segment:       build_segment(m_nx); break;
          case Geometry::Polytope::Type::Quadrilateral: build_quad(m_nx, m_ny); break;
          case Geometry::Polytope::Type::Triangle:      build_tri(m_nx, m_nx); break;
          case Geometry::Polytope::Type::Tetrahedron:   build_tet(m_nx, m_ny, m_nz); break;
          case Geometry::Polytope::Type::Wedge:         build_wedge(m_nx, m_nz); break;
        }
      }

      /**
       * @brief Builds quadrature for a point (0D).
       *
       * Creates a single quadrature point at the origin with weight 1.
       */
      void build_point()
      {
        m_points.clear();
        m_weights.resize(1);
        Math::SpatialVector<Real> p; p.resize(1); p[0]=0.0;
        m_points.push_back(std::move(p));
        m_weights[0]=1.0;
      }

      /**
       * @brief Builds 1D Gauss-Lobatto quadrature on a segment.
       * @param n Number of quadrature points
       *
       * Constructs @f$ n @f$ Gauss-Lobatto points on the reference segment
       * @f$ [0,1] @f$ including both endpoints.
       */
      void build_segment(size_t n)
      {
        std::vector<Real> x,w; gll_1d_unit(n,x,w);
        m_points.clear(); m_points.reserve(n); m_weights.resize(n);
        for (size_t i=0;i<n;++i){
          Math::SpatialVector<Real> p; p.resize(1); p[0]=x[i];
          m_points.push_back(std::move(p));
          m_weights[i]=w[i]; // length 1
        }
      }

      /**
       * @brief Builds Gauss-Lobatto quadrature on a quadrilateral.
       * @param nx Number of points in first direction
       * @param ny Number of points in second direction
       *
       * Constructs a tensor-product quadrature rule on the reference
       * quadrilateral @f$ [0,1] \times [0,1] @f$ using 1D Gauss-Lobatto
       * rules in each direction. Total number of points: @f$ nx \times ny @f$.
       */
      void build_quad(size_t nx, size_t ny)
      {
        std::vector<Real> x,wx,y,wy;
        gll_1d_unit(nx,x,wx); gll_1d_unit(ny,y,wy);
        const size_t N = nx*ny;
        m_points.clear(); m_points.reserve(N); m_weights.resize(N);
        size_t k=0;
        for(size_t j=0;j<ny;++j)
          for(size_t i=0;i<nx;++i){
            Math::SpatialVector<Real> p; p.resize(2); p[0]=x[i]; p[1]=y[j];
            m_points.push_back(std::move(p));
            m_weights[k++] = wx[i]*wy[j]; // area 1
          }
      }

      /**
       * @brief Builds Gauss-Lobatto quadrature on a triangle.
       * @param nu Number of points in first parametric direction
       * @param nv Number of points in second parametric direction
       *
       * Uses the Duffy transformation to map tensor-product Gauss-Lobatto
       * points to the reference triangle. The transformation is:
       * @f[
       *   (r,s) = (u, (1-u)v)
       * @f]
       * with Jacobian @f$ (1-u) @f$. The reference triangle has vertices
       * at @f$ (0,0) @f$, @f$ (1,0) @f$, and @f$ (0,1) @f$ with area 1/2.
       */
      void build_tri(size_t nu, size_t nv)
      {
        // Duffy with 1D GLL in each param:
        // (r,s)=(u,(1-u)v),  J=(1-u)
        std::vector<Real> u,wu,v,wv;
        gll_1d_unit(nu,u,wu); gll_1d_unit(nv,v,wv);
        const size_t N = nu*nv;
        m_points.clear(); m_points.reserve(N); m_weights.resize(N);
        size_t k=0;
        for(size_t j=0;j<nv;++j)
          for(size_t i=0;i<nu;++i){
            const Real r = u[i];
            const Real s = (1.0 - u[i]) * v[j];
            Math::SpatialVector<Real> p; p.resize(2); p[0]=r; p[1]=s;
            m_points.push_back(std::move(p));
            m_weights[k++] = wu[i]*wv[j]*(1.0 - u[i]); // integrates to 1/2
          }
      }

      /**
       * @brief Builds Gauss-Lobatto quadrature on a tetrahedron.
       * @param nu Number of points in first parametric direction
       * @param nv Number of points in second parametric direction
       * @param nw Number of points in third parametric direction
       *
       * Uses the 3D Duffy transformation to map tensor-product Gauss-Lobatto
       * points to the reference tetrahedron. The transformation is:
       * @f[
       *   (r,s,t) = (u, (1-u)v, (1-u)(1-v)w)
       * @f]
       * with Jacobian @f$ (1-u)^2(1-v) @f$. The reference tetrahedron has
       * volume 1/6.
       */
      void build_tet(size_t nu, size_t nv, size_t nw)
      {
        // 3D Duffy with 1D GLL:
        // (r,s,t)=(u,(1-u)v,(1-u)(1-v)w),  J=(1-u)^2(1-v)
        std::vector<Real> u,wu,v,wv,w,ww;
        gll_1d_unit(nu,u,wu); gll_1d_unit(nv,v,wv); gll_1d_unit(nw,w,ww);
        const size_t N = nu*nv*nw;
        m_points.clear(); m_points.reserve(N); m_weights.resize(N);
        size_t k=0;
        for(size_t kk=0; kk<nw; ++kk)
          for(size_t j=0; j<nv; ++j)
            for(size_t i=0; i<nu; ++i){
              const Real r = u[i];
              const Real s = (1.0 - u[i]) * v[j];
              const Real t = (1.0 - u[i]) * (1.0 - v[j]) * w[kk];
              Math::SpatialVector<Real> p; p.resize(3); p[0]=r; p[1]=s; p[2]=t;
              m_points.push_back(std::move(p));
              m_weights[k++] = wu[i]*wv[j]*ww[kk] * (1.0 - u[i])*(1.0 - u[i])*(1.0 - v[j]); // volume 1/6
            }
      }

      /**
       * @brief Builds Gauss-Lobatto quadrature on a wedge (prism).
       * @param ntri Number of points for the triangular cross-section
       * @param nz Number of points in the extrusion direction
       *
       * Constructs a quadrature rule on the reference wedge by combining:
       * - A Duffy-transformed triangular rule (using @p ntri points)
       * - A 1D Gauss-Lobatto rule in the extrusion direction (using @p nz points)
       *
       * The reference wedge has volume 1/2.
       */
      void build_wedge(size_t ntri, size_t nz)
      {
        // triangle (Duffy with GLL) × segment (GLL)
        std::vector<Real> u,wu,v,wv,z,wz;
        gll_1d_unit(ntri,u,wu); gll_1d_unit(ntri,v,wv); gll_1d_unit(nz,z,wz);
        const size_t N = ntri*ntri*nz;
        m_points.clear(); m_points.reserve(N); m_weights.resize(N);
        size_t k=0;
        for(size_t kk=0; kk<nz; ++kk)
          for(size_t j=0; j<ntri; ++j)
            for(size_t i=0; i<ntri; ++i){
              const Real r = u[i];
              const Real s = (1.0 - u[i]) * v[j];
              const Real t = z[kk];
              Math::SpatialVector<Real> p; p.resize(3); p[0]=r; p[1]=s; p[2]=t;
              m_points.push_back(std::move(p));
              m_weights[k++] = wu[i]*wv[j]*(1.0 - u[i]) * wz[kk]; // wedge volume = 1/2
            }
      }

      size_t m_nx{2}, m_ny{2}, m_nz{2}; ///< Number of points per direction
      std::vector<Math::SpatialVector<Real>> m_points; ///< Quadrature point coordinates
      Math::Vector<Real> m_weights; ///< Quadrature weights
  };
}

#endif
