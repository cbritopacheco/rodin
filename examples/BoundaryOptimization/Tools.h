#ifndef RODIN_EXAMPLES_BOUNDARYOPTIMIZATION_TOOLS_H
#define RODIN_EXAMPLES_BOUNDARYOPTIMIZATION_TOOLS_H

#include <vector>
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/Polytope.h"

namespace Rodin::Examples::BoundaryOptimization
{
  template <class TopologicalFunction>
  std::vector<Geometry::Point> locations(const TopologicalFunction& tf)
  {
    const Real tc = tf.max();
    std::vector<Geometry::Point> cs;
    for (auto it = tf.getFiniteElementSpace().getMesh().getVertex(); !it.end(); ++it)
    {
      const Geometry::Point p(*it, it->getTransformation(),
          Geometry::Polytope::getVertices(
            Geometry::Polytope::Type::Point).col(0), it->getCoordinates());
      const Real tp = tf(p);
      if (tp > 1e-12 && (tp / tc) > (1 - 1e-12))
      {
        cs.emplace_back(std::move(p));
        break;
      }
    }
    return cs;
  }

  template <class DistanceFunction>
  void holes(Real radius, DistanceFunction& dist, const std::vector<Geometry::Point>& cs)
  {
    if (cs.size())
    {
      auto insert =
        [&](const Geometry::Point& v)
        {
          Real d = dist(v);
          for (const auto& c : cs)
          {
            const Real dd = (v - c).norm() - radius;
            d = std::min(d, dd);
          }
          return d;
        };
      dist = insert;
    }
  }

  inline
  size_t rmc(Geometry::MeshBase& mesh,
      const FlatSet<Geometry::Attribute>& attrs, Geometry::Attribute a, Real tol = 1e-5)
  {
    const size_t per = mesh.getPerimeter();
    const size_t D = mesh.getDimension();
    auto ccl = mesh.ccl(
        [](const Geometry::Polytope& p1, const Geometry::Polytope& p2)
        {
          return p1.getAttribute() == p2.getAttribute();
        }, D - 1, attrs );
    size_t ccs = ccl.getCount();

    for (const auto& cc : ccl)
    {
      Real area = 0;
      for (const Index i : cc)
        area += mesh.getFace(i)->getMeasure();
      if ((area / per) < tol)
      {
        for (const Index i : cc)
        {
          if (attrs.contains(mesh.getFace(i)->getAttribute()))
            mesh.setAttribute({ D - 1, i }, a);
        }
        ccs--;
      }
    }
    return ccs;
  }
}

#endif
