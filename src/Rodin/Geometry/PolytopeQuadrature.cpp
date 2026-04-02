/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "PolytopeQuadrature.h"
#include "Rodin/QF/QuadratureFormula.h"

namespace Rodin::Geometry
{
  PolytopeQuadrature::PolytopeQuadrature(
      const Polytope& polytope, const QF::QuadratureFormulaBase& qf)
    : m_qf(&qf)
  {
    m_ps.reserve(qf.getSize());
    for (size_t qp = 0; qp < qf.getSize(); ++qp)
      m_ps.emplace_back(polytope, qf.getPoint(qp));
  }
}
