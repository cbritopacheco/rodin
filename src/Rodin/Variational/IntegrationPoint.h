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
  class IntegrationPoint
  {
    public:
      IntegrationPoint(const Geometry::Point& p)
        : m_p(p), m_qp(0)
      {}

      IntegrationPoint(
          const Geometry::Point& p,
          const QF::QuadratureFormulaBase& qf,
          size_t qp)
        : m_p(p), m_qf(qf), m_qp(qp)
      {}

      const Geometry::Point& getPoint() const
      {
        return m_p;
      }

      bool hasQuadratureFormula() const
      {
        return m_qf.has_value();
      }

      const QF::QuadratureFormulaBase& getQuadratureFormula() const
      {
        assert(m_qf);
        return m_qf->get();
      }

      size_t getIndex() const
      {
        return m_qp;
      }

    private:
      std::reference_wrapper<const Geometry::Point> m_p;
      Optional<std::reference_wrapper<const QF::QuadratureFormulaBase>> m_qf;
      size_t m_qp;
  };
}

#endif
