/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Implementation of the GenericPolytopeQuadrature dispatcher class.
 */

#include "GenericPolytopeQuadrature.h"

#include "Centroid.h"
#include "GaussLegendre.h"
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
        // For points, use the single-point centroid rule
        m_qf = std::make_unique<Centroid>(g);
        break;
      }
      case Geometry::Polytope::Type::Segment:
      {
        // For segments, use Grundmann-Möller with appropriate s parameter
        const size_t i = (m_order / 2) * 2 + 1;
        m_qf = std::make_unique<GrundmannMoller>(i / 2, g);
        break;
      }
      case Geometry::Polytope::Type::Triangle:
      {
        // For triangles, use Grundmann-Möller with appropriate s parameter
        const size_t i = (m_order / 2) * 2 + 1;
        m_qf = std::make_unique<GrundmannMoller>(i / 2, g);
        break;
      }
      case Geometry::Polytope::Type::Quadrilateral:
      {
        // For quadrilaterals, use tensor-product Gauss-Legendre
        m_qf = std::make_unique<GaussLegendre>(g);
        break;
      }
      case Geometry::Polytope::Type::Tetrahedron:
      {
        // For tetrahedra, use Grundmann-Möller with appropriate s parameter
        const size_t i = (m_order / 2) * 2 + 1;
        m_qf = std::make_unique<GrundmannMoller>(i / 2, g);
        break;
      }
    }
  }
}
