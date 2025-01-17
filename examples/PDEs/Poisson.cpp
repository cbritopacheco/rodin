/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <RodinExternal/MMG.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::External;

int main(int, char**)
{
  MMG::Mesh mesh;
  mesh.load("oscar.mesh", IO::FileFormat::MEDIT);
  mesh.getConnectivity().compute(1, 2);
  // std::set<size_t> attrs;
  // size_t count54 = 0;
  // size_t count120 = 0;
  // for (auto edge = mesh.getFace(); edge; ++edge)
  // {
  //   if (edge->getAttribute() == 54)
  //     count54++;
  //   else if (edge->getAttribute() == 120)
  //     count120++;
  // }
  // std::cout << count54 << " " << count120 << std::endl;

  mesh.trace({{6, 666}});
  mesh.trace({{{6, 7}, 666}});

  mesh.save("test.mesh", IO::FileFormat::MEDIT);

  // // Build a mesh
  // Mesh mesh;
  // mesh = mesh.UniformGrid(Polytope::Type::Quadrilateral, { 32, 32 });
  // mesh.scale(1.0 / 31);
  // mesh.getConnectivity().compute(1, 2);

  // // Functions
  // P1 vh(mesh);

  // // auto f = cos(2 * M_PI * F::x) * sin(2 * M_PI * F::y);
  // ScalarFunction f(1.0);

  // TrialFunction u(vh);
  // TestFunction  v(vh);

  // // Define problem
  // Problem poisson(u, v);
  // poisson = Integral(Grad(u), Grad(v))
  //         - Integral(f, v)
  //         + DirichletBC(u, Zero());

  // // Solve
  // CG(poisson).solve();

  // // Save solution
  // u.getSolution().save("Poisson.gf");
  // mesh.save("Poisson.mesh");

  return 0;
}
