/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <Rodin/IO/XDMF.h>

using namespace Rodin;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main()
{
  Mesh mesh = Mesh().UniformGrid(Polytope::Type::Triangle, {16, 16});
  mesh.scale(1.0 / 15.0);
  mesh.getConnectivity().compute(1, 2);

  P1 vh(mesh);
  TrialFunction u(vh);
  TestFunction v(vh);
  RealFunction f = 1;

  Problem poisson(u, v);
  poisson = Integral(Grad(u), Grad(v)) - Integral(f, v) + DirichletBC(u, Zero());
  CG(poisson).solve();

  IO::XDMF xdmf("Poisson");
  xdmf.grid().setMesh(mesh).add("u", u.getSolution());
  xdmf.write().flush();
  return 0;
}
