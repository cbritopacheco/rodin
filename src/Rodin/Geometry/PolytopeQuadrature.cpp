/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "PolytopeQuadrature.h"
#include "PolytopeTransformation.h"
#include "Rodin/QF/QuadratureFormula.h"

namespace Rodin::Geometry
{
  PolytopeQuadrature::PolytopeQuadrature(
    const Polytope& polytope, const QF::QuadratureFormula& qf)
    : m_qf(&qf)
  {
    std::vector<Math::SpatialMatrix<Real>> jacobians;
    polytope.getTransformation().jacobian(jacobians, qf);
    const size_t count = jacobians.size();
    assert(count == qf.getSize());
    m_ps.reserve(count);
    for (size_t qp = 0; qp < count; ++qp)
      m_ps.emplace_back(jacobians[qp], polytope, qf, qp);
  }
}
