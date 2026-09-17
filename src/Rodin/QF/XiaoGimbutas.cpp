/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Coefficients of the XiaoGimbutas simplex rules.
 *
 * The published coefficients are vendored from the Xiao--Gimbutas triasymq
 * distribution in XiaoGimbutasData.h. Every table is checked against an
 * independent moment oracle.
 *
 * Each entry is a flat run of (x, y[, z], w) per node.
 *
 * These rules are generally asymmetric, which is a different construction
 * from the symmetric rules of WitherdenVincent.cpp and reaches different
 * counts. A table is therefore compared only with its own family.
 *
 */

#include "XiaoGimbutas.h"
#include "XiaoGimbutasData.h"

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
          return Data::XiaoGimbutas::triangle;
        case Geometry::Polytope::Type::Tetrahedron:
          return Data::XiaoGimbutas::tetrahedron;
        default:
          return s_empty;
      }
    }
  }

  size_t XiaoGimbutas::getMaxDegree(Geometry::Polytope::Type g)
  {
    return tableFor(g).size();
  }

  XiaoGimbutas::XiaoGimbutas(size_t degree, Geometry::Polytope::Type g)
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
