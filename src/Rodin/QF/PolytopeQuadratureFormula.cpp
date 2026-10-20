/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Implementation of the PolytopeQuadratureFormula dispatcher class.
 */

#include "PolytopeQuadratureFormula.h"

#include "Centroid.h"
#include "GaussLegendre.h"
#include "TensorProduct.h"
#include "WitherdenVincent.h"
#include "XiaoGimbutas.h"

namespace Rodin::QF
{
  std::unique_ptr<const QuadratureFormulaBase> PolytopeQuadratureFormula::build(
    size_t order, Geometry::Polytope::Type g)
  {
    switch (g)
    {
      case Geometry::Polytope::Type::Point:
        return std::make_unique<Centroid>(g);

      case Geometry::Polytope::Type::Segment:
      {
        const size_t n = std::max<size_t>(1, (order + 2) / 2);
        return std::make_unique<GaussLegendre>(g, n);
      }

      case Geometry::Polytope::Type::Triangle:
      {
        if (XiaoGimbutas::isAvailable(order, g))
          return std::make_unique<XiaoGimbutas>(order, g);
        const size_t n = std::max<size_t>(1, (order + 2) / 2);
        return std::make_unique<GaussLegendre>(g, n, n);
      }

      case Geometry::Polytope::Type::Quadrilateral:
      {
        const size_t n = std::max<size_t>(1, (order + 2) / 2);
        const GaussLegendre line(Geometry::Polytope::Type::Segment, n);
        return std::make_unique<TensorProduct>(g, line, line);
      }

      case Geometry::Polytope::Type::Tetrahedron:
      {
        if (XiaoGimbutas::isAvailable(order, g))
          return std::make_unique<XiaoGimbutas>(order, g);
        const size_t n = std::max<size_t>(1, (order + 2) / 2);
        return std::make_unique<GaussLegendre>(g, n, n, n);
      }

      case Geometry::Polytope::Type::Wedge:
      {
        const size_t n = std::max<size_t>(1, (order + 2) / 2);
        const GaussLegendre line(Geometry::Polytope::Type::Segment, n);
        if (XiaoGimbutas::isAvailable(order, Geometry::Polytope::Type::Triangle))
        {
          const XiaoGimbutas triangle(order, Geometry::Polytope::Type::Triangle);
          return std::make_unique<TensorProduct>(g, triangle, line);
        }
        const size_t nt = std::max<size_t>(1, (order + 2) / 2);
        const GaussLegendre triangle(Geometry::Polytope::Type::Triangle, nt, nt);
        return std::make_unique<TensorProduct>(g, triangle, line);
      }

      case Geometry::Polytope::Type::Pyramid:
      {
        if (WitherdenVincent::isAvailable(order, g))
          return std::make_unique<WitherdenVincent>(order, g);

        // The (1 - z)^2 Jacobian of the collapse now rides in the Gauss-Jacobi
        // weight rather than being evaluated at the nodes, so the integrand in
        // the collapsed coordinates is polynomial of degree p in each
        // direction and every direction needs ceil((p + 1) / 2) points.
        const size_t n = std::max<size_t>(1, (order + 2) / 2);
        return std::make_unique<GaussLegendre>(g, n, n, n);
      }

      case Geometry::Polytope::Type::Hexahedron:
      {
        const size_t n = std::max<size_t>(1, (order + 2) / 2);
        const GaussLegendre line(Geometry::Polytope::Type::Segment, n);
        return std::make_unique<TensorProduct>(g, line, line, line);
      }

      default:
      {
        assert(false);
        return {};
      }
    }
  }

  const QuadratureFormula& PolytopeQuadratureFormula::get(
    size_t order, Geometry::Polytope::Type g)
  {
    return QuadratureFormula::get<PolytopeQuadratureFormula>(order, g);
  }

  PolytopeQuadratureFormula::PolytopeQuadratureFormula(
    size_t order, Geometry::Polytope::Type g)
    : m_qf(build(order, g)),
      m_order(order),
      m_geometry(g)
  {}
}
