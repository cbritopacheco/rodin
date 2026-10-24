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

#include <type_traits>
#include <utility>

#include "Centroid.h"
#include "GaussLegendre.h"
#include "TensorProduct.h"
#include "WitherdenVincent.h"
#include "XiaoGimbutas.h"

namespace Rodin::QF
{
  template <class Pointer>
  Pointer PolytopeQuadratureFormula::buildFormula(
    size_t order, Geometry::Polytope::Type g)
  {
    const auto construct = []<class Formula>(auto&&... args) -> Pointer {
      if constexpr (std::is_same_v<Pointer, std::unique_ptr<const QuadratureFormulaBase>>)
        return std::make_unique<const Formula>(std::forward<decltype(args)>(args)...);
      else
      {
        static_assert(
          std::is_same_v<Pointer, std::shared_ptr<const QuadratureFormulaBase>>);
        return std::make_shared<const Formula>(std::forward<decltype(args)>(args)...);
      }
    };
    switch (g)
    {
      case Geometry::Polytope::Type::Point:
        return construct.template operator()<Centroid>(g);

      case Geometry::Polytope::Type::Segment:
      {
        const size_t n = std::max<size_t>(1, order / 2 + 1);
        return construct.template operator()<GaussLegendre>(g, n);
      }

      case Geometry::Polytope::Type::Triangle:
      {
        if (XiaoGimbutas::isAvailable(order, g))
          return construct.template operator()<XiaoGimbutas>(order, g);
        const size_t n = std::max<size_t>(1, order / 2 + 1);
        return construct.template operator()<GaussLegendre>(g, n, n);
      }

      case Geometry::Polytope::Type::Quadrilateral:
      {
        const size_t n = std::max<size_t>(1, order / 2 + 1);
        const GaussLegendre line(Geometry::Polytope::Type::Segment, n);
        return construct.template operator()<TensorProduct>(g, line, line);
      }

      case Geometry::Polytope::Type::Tetrahedron:
      {
        if (XiaoGimbutas::isAvailable(order, g))
          return construct.template operator()<XiaoGimbutas>(order, g);
        const size_t n = std::max<size_t>(1, order / 2 + 1);
        return construct.template operator()<GaussLegendre>(g, n, n, n);
      }

      case Geometry::Polytope::Type::Wedge:
      {
        const size_t n = std::max<size_t>(1, order / 2 + 1);
        const GaussLegendre line(Geometry::Polytope::Type::Segment, n);
        if (XiaoGimbutas::isAvailable(order, Geometry::Polytope::Type::Triangle))
        {
          const XiaoGimbutas triangle(order, Geometry::Polytope::Type::Triangle);
          return construct.template operator()<TensorProduct>(g, triangle, line);
        }
        const size_t nt = std::max<size_t>(1, order / 2 + 1);
        const GaussLegendre triangle(Geometry::Polytope::Type::Triangle, nt, nt);
        return construct.template operator()<TensorProduct>(g, triangle, line);
      }

      case Geometry::Polytope::Type::Pyramid:
      {
        if (WitherdenVincent::isAvailable(order, g))
          return construct.template operator()<WitherdenVincent>(order, g);

        // The (1 - z)^2 Jacobian of the collapse now rides in the Gauss-Jacobi
        // weight rather than being evaluated at the nodes, so the integrand in
        // the collapsed coordinates is polynomial of degree p in each
        // direction and every direction needs ceil((p + 1) / 2) points.
        const size_t n = std::max<size_t>(1, order / 2 + 1);
        return construct.template operator()<GaussLegendre>(g, n, n, n);
      }

      case Geometry::Polytope::Type::Hexahedron:
      {
        const size_t n = std::max<size_t>(1, order / 2 + 1);
        const GaussLegendre line(Geometry::Polytope::Type::Segment, n);
        return construct.template operator()<TensorProduct>(g, line, line, line);
      }

      default:
      {
        assert(false);
        return {};
      }
    }
  }

  std::unique_ptr<const QuadratureFormulaBase> PolytopeQuadratureFormula::build(
    size_t order, Geometry::Polytope::Type g)
  {
    return buildFormula<std::unique_ptr<const QuadratureFormulaBase>>(order, g);
  }

  const QuadratureFormula& PolytopeQuadratureFormula::get(
    size_t order, Geometry::Polytope::Type g)
  {
    return QuadratureFormula::get<PolytopeQuadratureFormula>(order, g);
  }

  PolytopeQuadratureFormula::PolytopeQuadratureFormula(
    size_t order, Geometry::Polytope::Type g)
    : m_qf(buildFormula<std::shared_ptr<const QuadratureFormulaBase>>(order, g)),
      m_order(order),
      m_geometry(g)
  {}
}
