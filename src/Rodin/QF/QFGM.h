/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_QF_QFGM_H
#define RODIN_VARIATIONAL_QF_QFGM_H

#define RODIN_QFGM_MAX_ORDER 32

#include "QuadratureFormula.h"

namespace Rodin::QF
{
  /**
   * @brief Grundmann-Möller Simplex Quadrature
   *
   * Implements the Grundmann-Möller integration rules as described in
   * @cite grundmann1978invariant.
   */
  class QFGM final : public QuadratureFormulaBase
  {
    static std::vector<std::array<Math::Vector, RODIN_QFGM_MAX_ORDER>> initializeWeights();
    static std::vector<std::array<Math::Matrix, RODIN_QFGM_MAX_ORDER>> initializePoints();

    public:
      /// Parent class
      using Parent = QuadratureFormulaBase;

      constexpr
      QFGM(size_t order, Geometry::Polytope::Geometry geom)
        : Parent(geom),
          m_degree((order + 1) / 2)
      {
        assert(order < RODIN_QFGM_MAX_ORDER);
        assert(Geometry::Polytope::isSimplex(geom));
      }

      constexpr
      QFGM(const QFGM& other)
        : Parent(other),
          m_degree(other.m_degree)
      {}

      constexpr
      QFGM(QFGM&& other)
        : Parent(std::move(other)),
          m_degree(std::move(other.m_degree))
      {}

      inline
      const Math::Matrix& getPoints() const override
      {
        return s_points[Geometry::Polytope::getGeometryDimension(getGeometry())][m_degree];
      }

      inline
      const Math::Vector& getWeights() const override
      {
        return s_weights[Geometry::Polytope::getGeometryDimension(getGeometry())][m_degree];
      }

    private:
      static const std::vector<std::array<Math::Vector, RODIN_QFGM_MAX_ORDER>> s_weights;
      static const std::vector<std::array<Math::Matrix, RODIN_QFGM_MAX_ORDER>> s_points;

      const size_t m_degree;
  };
}

#endif
