/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_INTEGRATIONPOINT_H
#define RODIN_VARIATIONAL_INTEGRATIONPOINT_H

#include <cassert>

#include "Rodin/Types.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/QF/QuadratureFormula.h"

namespace Rodin::Variational
{
  /**
   * @brief Evaluation point used by variational operators.
   *
   * An IntegrationPoint always stores a geometric point. It may optionally be
   * associated to a quadrature formula and quadrature index:
   * - `getQuadratureFormula() != nullptr`: the point is part of a quadrature
   *   loop, and `getIndex()` identifies the quadrature sample in that formula.
   * - `getQuadratureFormula() == nullptr`: the point represents a generic
   *   pointwise evaluation (outside a quadrature tabulation context). In this
   *   case, consumers are expected to evaluate using
   *   `getPoint().getReferenceCoordinates()` rather than quadrature tabulation.
   *
   * @note The quadrature formula pointer is non-owning. The referenced
   * quadrature object must outlive all IntegrationPoint uses that dereference
   * this pointer.
   */
  class IntegrationPoint
  {
    public:
      IntegrationPoint(const Geometry::Point& p)
        : m_p(p), m_qf(nullptr), m_qp(0)
      {}

      /**
       * @brief Constructs a quadrature-associated integration point.
       * @param[in] p Geometric point (usually built from quadrature sample coordinates)
       * @param[in] qf Non-owning pointer to the quadrature formula that owns
       *                the sample indexing, or @c nullptr for direct
       *                pointwise evaluation
       * @param[in] qp Quadrature sample index in @p qf
       */
      IntegrationPoint(
          const Geometry::Point& p,
          const QF::QuadratureFormulaBase* qf,
          size_t qp)
        : m_p(p), m_qf(qf), m_qp(qp)
      {}

      /**
       * @returns Underlying geometric point context.
       */
      const Geometry::Point& getPoint() const
      {
        return m_p;
      }

      /**
       * @returns Non-owning pointer to the quadrature formula, or @c nullptr
       * when this integration point is not in a quadrature context.
       */
      const QF::QuadratureFormulaBase* getQuadratureFormula() const
      {
        return m_qf;
      }

      /**
       * @returns Quadrature index associated with this point.
       * @note Meaningful only when getQuadratureFormula() is not @c nullptr.
       */
      size_t getIndex() const
      {
        return m_qp;
      }

    private:
      std::reference_wrapper<const Geometry::Point> m_p;
      const QF::QuadratureFormulaBase* m_qf;
      size_t m_qp;
  };
}

#endif
