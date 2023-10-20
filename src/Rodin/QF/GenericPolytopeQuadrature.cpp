/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "GenericPolytopeQuadrature.h"

#include "QF1P1.h"
#include "GrundmannMoller.h"

namespace Rodin::QF
{
  GenericPolytopeQuadrature::GenericPolytopeQuadrature(size_t order, Geometry::Polytope::Type g)
    : Parent(g),
      m_order(order)
  {
    switch (g)
    {
      case Geometry::Polytope::Type::Point:
      {
        m_qf = std::make_unique<QF1P1>(g);
        break;
      }
      case Geometry::Polytope::Type::Segment:
      {
        m_qf = std::make_unique<GrundmannMoller>(order, g);
        break;
      }
      case Geometry::Polytope::Type::Triangle:
      {
        m_qf = std::make_unique<GrundmannMoller>(order, g);
        break;
      }
      case Geometry::Polytope::Type::Quadrilateral:
      {
        m_qf = std::make_unique<QF1P1>(g);
        break;
      }
      case Geometry::Polytope::Type::Tetrahedron:
      {
        m_qf = std::make_unique<GrundmannMoller>(order, g);
        break;
      }
    }
  }
}
