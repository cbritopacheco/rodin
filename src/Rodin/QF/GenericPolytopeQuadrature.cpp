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
  struct Key
  {
    size_t order;

    Geometry::Polytope::Type g;

    friend bool operator<(const Key& a, const Key& b)
    {
      if (a.g != b.g) return a.g < b.g;
      return a.order < b.order;
    }

    friend bool operator==(const Key& a, const Key& b)
    {
      return a.order == b.order && a.g == b.g;
    }
  };

  std::unique_ptr<const QuadratureFormulaBase>
  GenericPolytopeQuadrature::build(size_t order, Geometry::Polytope::Type g)
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
        const size_t s = (order + 1) / 2;
        return std::make_unique<GrundmannMoller>(s, g);
      }

      case Geometry::Polytope::Type::Quadrilateral:
      {
        const size_t n = std::max<size_t>(1, (order + 2) / 2);
        return std::make_unique<GaussLegendre>(g, n, n);
      }

      case Geometry::Polytope::Type::Tetrahedron:
      {
        const size_t i = (order / 2) * 2 + 1;
        return std::make_unique<GrundmannMoller>(i / 2, g);
      }

      case Geometry::Polytope::Type::Wedge:
      {
        const size_t n = std::max<size_t>(1, (order + 2) / 2);
        return std::make_unique<GaussLegendre>(g, n, n);
      }

      case Geometry::Polytope::Type::Hexahedron:
      {
        const size_t n = std::max<size_t>(1, (order + 2) / 2);
        return std::make_unique<GaussLegendre>(g, n, n, n);
      }

      default:
      {
        assert(false);
        return {};
      }
    }
  }

  const QuadratureFormulaBase& GenericPolytopeQuadrature::get(size_t order, Geometry::Polytope::Type g)
  {
    static thread_local FlatMap<Key, std::unique_ptr<const QuadratureFormulaBase>> pool;

    const Key key{order, g};

    auto it = pool.find(key);
    if (it == pool.end())
      it = pool.emplace(key, build(order, g)).first;

    return *it->second;
  }

  // Optional: keep the old behavior, but delegate to build()
  GenericPolytopeQuadrature::GenericPolytopeQuadrature(size_t order, Geometry::Polytope::Type g)
    : Parent(g), m_qf(build(order, g)), m_order(order)
  {}
}
