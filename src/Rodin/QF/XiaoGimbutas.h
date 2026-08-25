/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_QF_XIAOGIMBUTAS_H
#define RODIN_QF_XIAOGIMBUTAS_H

/**
 * @file
 * @brief Defines the XiaoGimbutas simplex quadrature formula.
 */

#include "QuadratureFormula.h"

namespace Rodin::QF
{
  /**
   * @ingroup RodinQuadrature
   * @brief Near-minimal positive interior quadrature on simplices.
   *
   * Provides the tabulated positive interior rules described by Xiao and
   * Gimbutas @cite xiao2010numerical. The published coefficients are vendored
   * from the authors' triasymq distribution under its BSD 3-Clause license.
   * The triangle table reaches degree 50 and the tetrahedron table degree 15.
   * Every entry is independently tested for polynomial exactness, positive
   * weights, and interior nodes.
   *
   * Beyond the tabulated range, the default dispatcher uses a positive
   * Gauss--Jacobi conical-product rule.
   */
  class XiaoGimbutas final : public QuadratureFormulaBase
  {
    public:
      /// @brief Highest tabulated degree for @p g, or zero if unsupported.
      static size_t getMaxDegree(Geometry::Polytope::Type g);

      /// @brief Whether a rule of strength @p degree is tabulated for @p g.
      static bool isAvailable(size_t degree, Geometry::Polytope::Type g)
      {
        return degree >= 1 && degree <= getMaxDegree(g);
      }

      /**
       * @brief Constructs the rule of strength @p degree on @p g.
       * @pre isAvailable(degree, g)
       */
      XiaoGimbutas(size_t degree, Geometry::Polytope::Type g);

      /// @brief Copies a tabulated quadrature formula.
      XiaoGimbutas(const XiaoGimbutas&) = default;

      size_t getSize() const override
      {
        return m_count;
      }

      Real getWeight(size_t i) const override
      {
        assert(i < m_count);
        return m_data[i * (m_dimension + 1) + m_dimension];
      }

      const Math::SpatialVector<Real>& getPoint(size_t i) const override
      {
        assert(i < m_count);
        return m_points[i];
      }

      /// @brief The element this rule is defined on.
      Geometry::Polytope::Type getGeometry() const
      {
        return m_geometry;
      }

      XiaoGimbutas* copy() const noexcept override
      {
        return new XiaoGimbutas(*this);
      }

    private:
      Geometry::Polytope::Type m_geometry;
      size_t m_dimension;
      size_t m_count;
      const Real* m_data;
      std::vector<Math::SpatialVector<Real>> m_points;
  };
}

#endif
