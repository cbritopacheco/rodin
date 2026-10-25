/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Default concrete-family selection for persistent quadrature formulas.
 */

#include "QuadratureFormula.h"

#include <algorithm>
#include <cstdlib>

#include "Rodin/Alert/Exception.h"

#include "Centroid.h"
#include "GaussLegendre.h"
#include "TensorProduct.h"
#include "WitherdenVincent.h"
#include "XiaoGimbutas.h"

namespace Rodin::QF
{
  const QuadratureFormula& QuadratureFormula::get(
    size_t order, Geometry::Polytope::Type geometry)
  {
    using Key = std::tuple<size_t, Geometry::Polytope::Type>;
    const Key key{order, geometry};
    auto& lookup = getLookup<QuadratureFormula, Key>();
    for (const auto& entry : lookup.entries)
    {
      if (entry && entry->key == key)
        return *entry->pointer;
    }
    const auto& formula = select(order, geometry);
    using Entry = Lookup<Key>::Entry;
    lookup.entries[lookup.next] = Entry{key, &formula};
    lookup.next = (lookup.next + 1) % lookup.entries.size();
    return formula;
  }

  const QuadratureFormula& QuadratureFormula::select(
    size_t order, Geometry::Polytope::Type g)
  {
    switch (g)
    {
      case Geometry::Polytope::Type::Point:
        return get<Centroid>(g);

      case Geometry::Polytope::Type::Segment:
      {
        const size_t n = std::max<size_t>(1, order / 2 + 1);
        return get<GaussLegendre>(g, n);
      }

      case Geometry::Polytope::Type::Triangle:
      {
        if (XiaoGimbutas::isAvailable(order, g))
          return get<XiaoGimbutas>(order, g);
        const size_t n = std::max<size_t>(1, order / 2 + 1);
        return get<GaussLegendre>(g, n, n);
      }

      case Geometry::Polytope::Type::Quadrilateral:
      {
        const size_t n = std::max<size_t>(1, order / 2 + 1);
        const auto& line = get<GaussLegendre>(Geometry::Polytope::Type::Segment, n);
        return get<TensorProduct>(g, line, line);
      }

      case Geometry::Polytope::Type::Tetrahedron:
      {
        if (XiaoGimbutas::isAvailable(order, g))
          return get<XiaoGimbutas>(order, g);
        const size_t n = std::max<size_t>(1, order / 2 + 1);
        return get<GaussLegendre>(g, n, n, n);
      }

      case Geometry::Polytope::Type::Wedge:
      {
        const size_t n = std::max<size_t>(1, order / 2 + 1);
        const auto& line = get<GaussLegendre>(Geometry::Polytope::Type::Segment, n);
        if (XiaoGimbutas::isAvailable(order, Geometry::Polytope::Type::Triangle))
        {
          const auto& triangle =
            get<XiaoGimbutas>(order, Geometry::Polytope::Type::Triangle);
          return get<TensorProduct>(g, triangle, line);
        }
        const size_t nt = std::max<size_t>(1, order / 2 + 1);
        const auto& triangle =
          get<GaussLegendre>(Geometry::Polytope::Type::Triangle, nt, nt);
        return get<TensorProduct>(g, triangle, line);
      }

      case Geometry::Polytope::Type::Pyramid:
      {
        if (WitherdenVincent::isAvailable(order, g))
          return get<WitherdenVincent>(order, g);

        // The (1 - z)^2 Jacobian of the collapse now rides in the Gauss-Jacobi
        // weight rather than being evaluated at the nodes, so the integrand in
        // the collapsed coordinates is polynomial of degree p in each
        // direction and every direction needs ceil((p + 1) / 2) points.
        const size_t n = std::max<size_t>(1, order / 2 + 1);
        return get<GaussLegendre>(g, n, n, n);
      }

      case Geometry::Polytope::Type::Hexahedron:
      {
        const size_t n = std::max<size_t>(1, order / 2 + 1);
        const auto& line = get<GaussLegendre>(Geometry::Polytope::Type::Segment, n);
        return get<TensorProduct>(g, line, line, line);
      }

      default:
      {
        Alert::Exception() << __func__ << ": Unsupported quadrature geometry."
                           << Alert::Raise;
        std::abort();
      }
    }
  }

}
