/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_QF_QFGG_H
#define RODIN_VARIATIONAL_QF_QFGG_H

#define RODIN_VARIATIONAL_QF_QFGG_DEFAULT_ORDER 1

#include "QuadratureFormula.h"

#include "QF1P1.h"

namespace Rodin::QF
{
  /**
   * @brief Quadrature on a polytope with any of @ref
   * Geometry::Polytope::Geometry "Rodin's supported geometries".
   * @see @ref Geometry::Polytope::Geometry "Polytope::Geometry"
   */
  class QFGG : public QuadratureFormulaBase
  {
    public:
      /// Parent class
      using Parent = QuadratureFormulaBase;

      constexpr
      QFGG(Geometry::Polytope::Type g)
        : QFGG(RODIN_VARIATIONAL_QF_QFGG_DEFAULT_ORDER, g)
      {}

      constexpr
      QFGG(size_t order, Geometry::Polytope::Type g)
        : Parent(g),
          m_qf1p1(g),
          m_order(order)
      {}

      constexpr
      QFGG(const QFGG& other)
        : Parent(other),
          m_qf1p1(other.m_qf1p1),
          m_order(other.m_order)
      {}

      constexpr
      QFGG(QFGG&& other)
        : Parent(std::move(other)),
          m_qf1p1(std::move(other.m_qf1p1)),
          m_order(std::move(other.m_order))
      {}

      inline
      size_t getSize() const override
      {
        return m_qf1p1.getSize();
      }

      inline
      Scalar getWeight(size_t i) const override
      {
        return m_qf1p1.getWeight(i);
      }

      inline
      const Math::SpatialVector& getPoint(size_t i) const override
      {
        return m_qf1p1.getPoint(i);
      }

    private:
      const QF1P1 m_qf1p1;
      const size_t m_order;
  };
}

#endif
