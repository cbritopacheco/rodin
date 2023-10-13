/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

inline
Scalar K(const Geometry::Point& x, const Geometry::Point& y)
{
  return log((x.asVector() - y.asVector()).norm());
}


int main(int, char**)
{
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, 16, 16);
  P1 fes(mesh);

  ScalarFunction f = [](const Geometry::Point& p) { return 1; };
  Potential sl(K, f);

  GridFunction u(fes);
  u = sl;

  auto it = mesh.getElement();
  const auto& el = *it;
  const Math::SpatialVector rc{{0, 0}};
  Geometry::Point p(el, el.getTransformation(), rc);

  u.save("u.gf");
  mesh.save("u.mesh");

  return 0;
}

