/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_QF_QFGM_H
#define RODIN_VARIATIONAL_QF_QFGM_H

#define RODIN_QFGM_MAX_S 8
#define RODIN_QFGM_MAX_ORDER 2 * RODIN_QFGM_MAX_S + 1

#include <boost/multi_array.hpp>

#include "QuadratureFormula.h"

namespace Rodin::QF
{
  /**
   * @brief Grundmann-Möller Simplex Quadrature
   *
   * Implements the Grundmann-Möller integration rules as described in
   * @cite grundmann1978invariant.
   */
  class QFGM : public QuadratureFormulaBase
  {
    public:
      /// Parent class
      using Parent = QuadratureFormulaBase;

      constexpr
      QFGM(size_t s, Geometry::Polytope::Type geom);

      constexpr
      QFGM(const QFGM& other)
        : Parent(other),
          m_s(other.m_s),
          m_n(other.m_n),
          m_order(other.m_order)
      {}

      constexpr
      QFGM(QFGM&& other)
        : Parent(std::move(other)),
          m_s(std::move(other.m_s)),
          m_n(std::move(other.m_n)),
          m_order(std::move(other.m_order))
      {}

      size_t getSize() const override;

      Scalar getWeight(size_t i) const override;

      const Math::SpatialVector& getPoint(size_t i) const override;

      static boost::multi_array<size_t, 2> initSizes();
      static boost::multi_array<Math::Vector, 2> initWeights();
      static boost::multi_array<std::vector<Math::SpatialVector>, 2> initPoints();
    private:

      static boost::multi_array<size_t, 2> s_sizes;
      static boost::multi_array<Math::Vector, 2> s_weights;
      static boost::multi_array<std::vector<Math::SpatialVector>, 2> s_points;

      const size_t m_s;
      const size_t m_n;
      const size_t m_order;
  };
}

#endif
