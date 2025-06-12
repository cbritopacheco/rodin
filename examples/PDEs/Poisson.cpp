/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  // Build a mesh
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Quadrilateral, { 32, 32 });

  mesh.scale(1. / 31.0);
  mesh.getConnectivity().compute(1, 2);

  P1 vh(mesh);
  ScalarFunction f =
    [](const Point& p) {
      return cos(10 * M_PI * p.x());
    };

  TrialFunction u(vh);
  TestFunction  v(vh);

  // Define problem
  Problem poisson(u, v);
  poisson = Integral(Grad(u), Grad(v))
          - Integral(f, v)
          + DirichletBC(u, Zero());

  // Solve
  CG(poisson).solve();

  // Save solution
  u.getSolution().save("Poisson.gf");
  mesh.save("Poisson.mesh");

  return 0;
}
