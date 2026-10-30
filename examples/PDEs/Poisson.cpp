/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [poisson-source]
#include <cstdlib>

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
  //! [poisson-mesh]
  constexpr size_t nodes = 16;
  Mesh mesh = Mesh().UniformGrid(Polytope::Type::Triangle, {nodes, nodes});
  mesh.scale(1.0 / (nodes - 1));
  mesh.getConnectivity().compute(1, 2);
  //! [poisson-mesh]

  //! [poisson-form]
  P1 vh(mesh);
  TrialFunction u(vh);
  TestFunction v(vh);
  RealFunction f = 1;

  Problem poisson(u, v);
  poisson = Integral(Grad(u), Grad(v)) - Integral(f, v)
             + DirichletBC(u, Zero());
  //! [poisson-form]

  //! [poisson-solve]
  CG solver(poisson);
  solver.solve();
  if (!solver.success())
    return EXIT_FAILURE;
  //! [poisson-solve]

  //! [poisson-export]
  IO::XDMF xdmf("Poisson");
  xdmf.grid().setMesh(mesh).add("u", u.getSolution());
  xdmf.write().flush();
  //! [poisson-export]
  return 0;
}
//! [poisson-source]
