/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_TESTS_QF_QUADRATURESELECTION_H
#define RODIN_TESTS_QF_QUADRATURESELECTION_H

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <memory>

#include "Rodin/QF/Centroid.h"
#include "Rodin/QF/GaussLegendre.h"
#include "Rodin/QF/TensorProduct.h"
#include "Rodin/QF/XiaoGimbutas.h"
#include "Rodin/QF/WitherdenVincent.h"

namespace Rodin::Tests::QF
{
  using namespace Rodin::QF;

  /**
   * @brief Independently constructs the documented family for numerical comparisons.
   * @param order Required degree of polynomial exactness.
   * @param g Reference geometry.
   * @returns An owned concrete rule with the documented default samples.
   */
  inline std::unique_ptr<const QuadratureFormulaBase> makeQuadratureFormula(
    size_t order, Geometry::Polytope::Type g)
  {
    switch (g)
    {
      case Geometry::Polytope::Type::Point:
        return std::make_unique<const Centroid>(g);

      case Geometry::Polytope::Type::Segment:
      {
        const size_t n = std::max<size_t>(1, order / 2 + 1);
        return std::make_unique<const GaussLegendre>(g, n);
      }

      case Geometry::Polytope::Type::Triangle:
      {
        if (XiaoGimbutas::isAvailable(order, g))
          return std::make_unique<const XiaoGimbutas>(order, g);
        const size_t n = std::max<size_t>(1, order / 2 + 1);
        return std::make_unique<const GaussLegendre>(g, n, n);
      }

      case Geometry::Polytope::Type::Quadrilateral:
      {
        const size_t n = std::max<size_t>(1, order / 2 + 1);
        const GaussLegendre line(Geometry::Polytope::Type::Segment, n);
        return std::make_unique<const TensorProduct>(g, line, line);
      }

      case Geometry::Polytope::Type::Tetrahedron:
      {
        if (XiaoGimbutas::isAvailable(order, g))
          return std::make_unique<const XiaoGimbutas>(order, g);
        const size_t n = std::max<size_t>(1, order / 2 + 1);
        return std::make_unique<const GaussLegendre>(g, n, n, n);
      }

      case Geometry::Polytope::Type::Wedge:
      {
        const size_t n = std::max<size_t>(1, order / 2 + 1);
        const GaussLegendre line(Geometry::Polytope::Type::Segment, n);
        if (XiaoGimbutas::isAvailable(order, Geometry::Polytope::Type::Triangle))
        {
          const XiaoGimbutas triangle(order, Geometry::Polytope::Type::Triangle);
          return std::make_unique<const TensorProduct>(g, triangle, line);
        }
        const size_t nt = std::max<size_t>(1, order / 2 + 1);
        const GaussLegendre triangle(Geometry::Polytope::Type::Triangle, nt, nt);
        return std::make_unique<const TensorProduct>(g, triangle, line);
      }

      case Geometry::Polytope::Type::Pyramid:
      {
        if (WitherdenVincent::isAvailable(order, g))
          return std::make_unique<const WitherdenVincent>(order, g);

        // The (1 - z)^2 Jacobian of the collapse now rides in the Gauss-Jacobi
        // weight rather than being evaluated at the nodes, so the integrand in
        // the collapsed coordinates is polynomial of degree p in each
        // direction and every direction needs ceil((p + 1) / 2) points.
        const size_t n = std::max<size_t>(1, order / 2 + 1);
        return std::make_unique<const GaussLegendre>(g, n, n, n);
      }

      case Geometry::Polytope::Type::Hexahedron:
      {
        const size_t n = std::max<size_t>(1, order / 2 + 1);
        const GaussLegendre line(Geometry::Polytope::Type::Segment, n);
        return std::make_unique<const TensorProduct>(g, line, line, line);
      }

      default:
      {
        assert(false);
        return {};
      }
    }
  }

  /**
   * @brief Copies concrete family samples without a production dispatcher wrapper.
   * @param source Owned default-family rule to copy.
   * @returns An independent concrete copy of the source rule.
   */
  inline std::unique_ptr<const QuadratureFormulaBase> copyQuadratureFormula(
    const QuadratureFormulaBase& source)
  {
    if (const auto* formula = dynamic_cast<const Centroid*>(&source))
      return std::make_unique<const Centroid>(*formula);
    if (const auto* formula = dynamic_cast<const GaussLegendre*>(&source))
      return std::make_unique<const GaussLegendre>(*formula);
    if (const auto* formula = dynamic_cast<const TensorProduct*>(&source))
      return std::make_unique<const TensorProduct>(*formula);
    if (const auto* formula = dynamic_cast<const XiaoGimbutas*>(&source))
      return std::make_unique<const XiaoGimbutas>(*formula);
    if (const auto* formula = dynamic_cast<const WitherdenVincent*>(&source))
      return std::make_unique<const WitherdenVincent>(*formula);
    std::abort();
  }
}

#endif
