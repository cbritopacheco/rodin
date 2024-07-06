/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Variational/LinearElasticity.h>

#include <Rodin/Assembly/Sequential.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int argc, char** argv)
{
  // Load mesh
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 16, 16 });
  mesh.scale(1. / (15));
  mesh.displace(VectorFunction{-1, -1});
  mesh.getConnectivity().compute(1, 2);

  // Functions
  P1 vh(mesh, 2);
  P0 ph(mesh);

  TrialFunction u(vh);
  TestFunction  v(vh);

  TrialFunction p(ph);
  TestFunction  q(ph);

  Solver::CG cg;

  RealFunction f = [](const Geometry::Point& p) { return -exp(p.x()) * sin(p.y()); };

  auto n = BoundaryNormal(mesh);

  Problem darcy(u, p, v, q);
  darcy = Integral(u, v)
        - Integral(p, Div(v))
        - Integral(Div(u), q)
        - BoundaryIntegral(f, Dot(n, v));
  darcy.solve(cg);

  u.getSolution().save("u.gf");
  mesh.save("u.mesh");

  return 0;
}


