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
        const size_t n = std::max<size_t>(1, (m_order + 2) / 2);
        m_qf = std::make_unique<GaussLegendre>(g, n);
        break;
      }
      case Geometry::Polytope::Type::Triangle:
      {
        const size_t s = (m_order + 1) / 2;
        m_qf = std::make_unique<GrundmannMoller>(s, g);
        break;
      }
      case Geometry::Polytope::Type::Quadrilateral:
      {
        // For quadrilaterals, use tensor-product Gauss-Legendre
        const size_t n = std::max<size_t>(1, (m_order + 2) / 2); // ceil((order+1)/2)
        m_qf = std::make_unique<GaussLegendre>(g, n, n);
        break;
      }
      case Geometry::Polytope::Type::Tetrahedron:
      {
        // For tetrahedra, use Grundmann-Möller with appropriate s parameter
        const size_t i = (m_order / 2) * 2 + 1;
        m_qf = std::make_unique<GrundmannMoller>(i / 2, g);
        break;
      }
      case Geometry::Polytope::Type::Wedge:
      {
        // For wedges, tensor-product: triangle (Duffy) × segment
        // Use same n in all directions; 2n - 1 >= m_order
        const size_t n = std::max<size_t>(1, (m_order + 1) / 2);
        // GaussLegendre(Polytope::Type::Wedge, ntri, nz)
        m_qf = std::make_unique<GaussLegendre>(g, n, n);
        break;
      }
      case Geometry::Polytope::Type::Hexahedron:
      {
        const size_t n = std::max<size_t>(1, (m_order + 2) / 2);
        m_qf = std::make_unique<GaussLegendre>(g, n, n, n);
        break;
      }
    }
  }
}
