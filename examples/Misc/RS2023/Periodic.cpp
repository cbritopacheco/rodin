/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Math.h>
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

static constexpr Scalar m = 10;
static constexpr Scalar pi = Math::Constants::pi();
static constexpr Attribute dQ = 2;

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
  ScalarFunction phi = [](const Point& p) { return p.x(); };

  ScalarFunction gamma =
    [](const Point& p)
    {
      return 2 + sin(pi * m * p.x()) * cos(pi * m * p.y());
    };

  GridFunction gf(vh);
  gf = gamma;
  mesh.save("Q.mesh", IO::FileFormat::MFEM);
  gf.save("gamma.gf");

  Problem poisson(u, v);
  poisson = Integral(gamma * Grad(u), Grad(v))
          - Integral(v)
          + DirichletBC(u, phi).on(dQ);

  // Solve the problem
  Solver::SparseLU solver;
  poisson.solve(solver);

  // Save solution
  u.getSolution().save("RS2023.gf");
  mesh.save("RS2023.mesh");

  return 0;
}


