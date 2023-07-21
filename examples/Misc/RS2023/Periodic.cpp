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
  mesh.save("Q.mesh", IO::FileFormat::MFEM);

  // Functions
  P1 vh(mesh);

  TrialFunction u(vh);
  TestFunction  v(vh);

  // Define problem
  ScalarFunction phi = 5;

  ScalarFunction gamma =
    [](const Point& p)
    {
      return 2 + sin(pi * m * p.x()) * cos(pi * m * p.y());
    };

  GridFunction conductivity(vh);
  conductivity = gamma;
  conductivity.save("Conductivity.gf");

  Problem poisson(u, v);
  poisson = Integral(gamma * Grad(u), Grad(v))
          + DirichletBC(u, phi).on(dQ);

  // Solve the problem
  Solver::SparseLU solver;
  poisson.solve(solver);

  u.getSolution().save("Periodic.gf");

  P1 th(mesh, 2);
  GridFunction grad(th);
  grad = Grad(u.getSolution());
  grad.save("Gradient.gf");


  return 0;
}


