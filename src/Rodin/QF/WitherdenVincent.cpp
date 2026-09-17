/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Coefficients of the WitherdenVincent symmetric rules.
 *
 * The published coefficients are vendored from the Witherden--Vincent tables
 * distributed by PyFR in WitherdenVincentData.h. Every table is checked
 * against an independent moment oracle and against properties that use no
 * oracle.
 *
 * @par Counts
 * Published counts are tabulated in tests/unit/Rodin/QF/PublishedCounts.h and
 * asserted by the tests.
 *
 * Each entry is a flat run of (x, y[, z], w) per node.
 */

#include "WitherdenVincent.h"
#include "WitherdenVincentData.h"

namespace Rodin::QF
{
  namespace
  {
    const std::vector<std::vector<Real>>& tableFor(Geometry::Polytope::Type g)
    {
      static const std::vector<std::vector<Real>> s_empty;
      switch (g)
      {
        case Geometry::Polytope::Type::Triangle:
          return Data::WitherdenVincent::triangle;
        case Geometry::Polytope::Type::Quadrilateral:
          return Data::WitherdenVincent::quadrilateral;
        case Geometry::Polytope::Type::Tetrahedron:
          return Data::WitherdenVincent::tetrahedron;
        case Geometry::Polytope::Type::Wedge:
          return Data::WitherdenVincent::wedge;
        case Geometry::Polytope::Type::Pyramid:
          return Data::WitherdenVincent::pyramid;
        case Geometry::Polytope::Type::Hexahedron:
          return Data::WitherdenVincent::hexahedron;
        default:
          return s_empty;
      }
    }
  }

  size_t WitherdenVincent::getMaxDegree(Geometry::Polytope::Type g)
  {
    return tableFor(g).size();
  }

  WitherdenVincent::WitherdenVincent(size_t degree, Geometry::Polytope::Type g)
    : m_geometry(g),
      m_dimension(Geometry::Polytope::Traits(g).getDimension())
  {
    assert(isAvailable(degree, g));
    const auto& entry = tableFor(g)[degree - 1];
    m_data = entry.data();
    m_count = entry.size() / (m_dimension + 1);
    m_points.reserve(m_count);
    for (size_t i = 0; i < m_count; ++i)
    {
      Math::SpatialVector<Real> p;
      p.resize(static_cast<Eigen::Index>(m_dimension));
      for (size_t k = 0; k < m_dimension; ++k)
        p[static_cast<Eigen::Index>(k)] = entry[i * (m_dimension + 1) + k];
      m_points.push_back(std::move(p));
    }
  }
}
