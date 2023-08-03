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

static constexpr Attribute dQ = 2;
static constexpr Attribute dR = 3;

int main(int, char**)
{
  // Build a mesh
  Mesh mesh;
  mesh.load("Q.medit.mesh", IO::FileFormat::MEDIT);
  mesh.getConnectivity().compute(1, 2);

  // Functions
  P1 vh(mesh);

  TrialFunction u(vh);
  TestFunction  v(vh);

  // Define problem
  ScalarFunction f = 1.0;
  ScalarFunction zero = 0.0;

  ScalarFunction gamma =
    [](const Geometry::Point& p)
    {
      if (p.getPolytope().getAttribute() == 1)
        return 100.0;
      else
        return 1.0;
    };

  Problem poisson(u, v);
  poisson = Integral(gamma * Grad(u), Grad(v))
          - Integral(f, v)
          + DirichletBC(u, zero).on(dQ);

  // Solve the problem
  Solver::SparseLU solver;
  poisson.solve(solver);

  P1 th(mesh, 2);
  GridFunction grad(th);
  grad = Grad(u.getSolution());
  grad.save("Gradient.gf");

  // Save solution
  u.getSolution().save("RS2023.gf");
  mesh.save("RS2023.mesh");

  return 0;
}

