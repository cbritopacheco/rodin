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

static constexpr Geometry::Attribute Boundary = 2;

static constexpr Scalar m = 1;

int main(int, char**)
{
  const size_t n = 200;

  // Build a mesh
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
  mesh.getConnectivity().compute(1, 2);
  mesh.scale(1.0 / (n - 1));
  mesh.scale(2 * M_PI);

  // Functions
  P1 vh(mesh);

  TrialFunction u(vh);
  TestFunction  v(vh);

  // Define problem
  ScalarFunction f =
    [](const Point& p)
    {
      return sin(p.x() + M_PI / 4) * cos(p.y() + M_PI / 4);
    };

  ScalarFunction g = 0.0;

  IndexMap<IndexSet> dofs;
  for (Index i = 0; i < vh.getSize(); i += n)
    dofs[i].insert(i + n - 1);

  for (Index i = 0; i < n; i++)
    dofs[i].insert(i + n * (n - 1));

  Problem poisson(u, v);
  poisson = Integral(Grad(u), Grad(v))
          - Integral(f, v)
          + PeriodicBC(u, dofs);
  poisson.assemble();

  // Solve the problem
  Solver::SparseLU solver;
  poisson.solve(solver);

  // Save solution
  u.getSolution().save("Periodic.gf");
  mesh.save("Periodic.mesh");

  return 0;
}


