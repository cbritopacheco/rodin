#ifndef RODIN_VARIATIONAL_QF_GAUSSLEGENDRE_H
#define RODIN_VARIATIONAL_QF_GAUSSLEGENDRE_H

#include <vector>
#include <cassert>

#include "QuadratureFormula.h"

namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   * @brief Gauss–Legendre quadrature formula on reference polytopes.
   *
   * Supports Point, Segment, Triangle, Quadrilateral, Tetrahedron, and Wedge.
   * Order can be specified per-dimension.
   */
  class GaussLegendre final : public QuadratureFormulaBase
  {
    public:
      using Parent = QuadratureFormulaBase;

      /**
       * @brief Construct with uniform order for all dimensions.
       * @param g Geometry type
       * @param order Number of Gauss–Legendre points per 1D direction
       */
      GaussLegendre(Geometry::Polytope::Type g, size_t order)
        : Parent(g), m_nx(order), m_ny(order), m_nz(order)
      {
        assert(order >= 1);
        build();
      }

      /**
       * @brief Construct with separate orders in two directions.
       * @param g Geometry type
       * @param nx Order in first direction
       * @param ny Order in second direction
       */
      GaussLegendre(Geometry::Polytope::Type g, size_t nx, size_t ny)
        : Parent(g), m_nx(nx), m_ny(ny), m_nz(ny)
      {
        assert(nx >= 1 && ny >= 1);
        build();
      }

      /**
       * @brief Construct with separate orders in three directions.
       * @param g Geometry type
       * @param nu Order in first direction
       * @param nv Order in second direction
       * @param nw Order in third direction
       */
      GaussLegendre(Geometry::Polytope::Type g, size_t nu, size_t nv, size_t nw)
        : Parent(g), m_nx(nu), m_ny(nv), m_nz(nw)
      {
        assert(nu >= 1 && nv >= 1 && nw >= 1);
        build();
      }

      /**
       * @brief Construct with default order 2 in all directions.
       * @param g Geometry type
       */
      GaussLegendre(Geometry::Polytope::Type g)
        : Parent(g), m_nx(2), m_ny(2), m_nz(2)
      {
        build();
      }

      /**
       * @brief Number of quadrature points.
       */
      size_t getSize() const override
      {
        return m_points.size();
      }

      /**
       * @brief Coordinates of the i-th quadrature point.
       * @param i Index of point
       */
      const Math::SpatialVector<Real>& getPoint(size_t i) const override
      {
        assert(i < m_points.size());
        return m_points[i];
      }

      /**
       * @brief Weight of the i-th quadrature point.
       * @param i Index of point
       */
      Real getWeight(size_t i) const override
      {
        assert(i < static_cast<size_t>(m_weights.size()));
        return m_weights.coeff(i);
      }

      /**
       * @brief Clone this quadrature rule.
       */
      GaussLegendre* copy() const noexcept override
      {
        return new GaussLegendre(*this);
      }

    private:
      static void gl_1d_unit(size_t n, std::vector<Real>& x, std::vector<Real>& w, size_t maxIt = 100, Real tol = 1e-14);

      void build();
      void build_point();
      void build_segment(size_t n);
      void build_quad(size_t nx, size_t ny);
      void build_tri(size_t nu, size_t nv);
      void build_tet(size_t nu, size_t nv, size_t nw);
      void build_wedge(size_t ntri, size_t nz);

      size_t m_nx { 2 }, m_ny { 2 }, m_nz { 2 };
      std::vector<Math::SpatialVector<Real>> m_points;
      Math::Vector<Real> m_weights;
      size_t m_maxIt;
      Real m_tol;
  };
}

#endif
