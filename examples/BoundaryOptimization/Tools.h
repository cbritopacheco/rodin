#ifndef RODIN_EXAMPLES_BOUNDARYOPTIMIZATION_TOOLS_H
#define RODIN_EXAMPLES_BOUNDARYOPTIMIZATION_TOOLS_H

#include <vector>
#include "Rodin/Geometry/Mesh.h"
#include "Rodin/Geometry/Point.h"
#include "Rodin/Geometry/Polytope.h"

namespace Rodin::Examples::BoundaryOptimization
{
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
