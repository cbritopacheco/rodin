/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/IO/ForwardDecls.h"
#include "Rodin/IO/XDMF.h"
#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <type_traits>

using namespace Rodin;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  Mesh mesh;
  mesh.load("malla_merge.mesh", IO::FileFormat::MEDIT);
  // mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 16, 16 });
  mesh.getConnectivity().compute(2, 3); // Compute boundary
  mesh.getConnectivity().compute(2, 1); // Compute boundary
  mesh.getConnectivity().compute(1, 0); // Compute boundary

  H1 vh(std::integral_constant<size_t, 2>{}, mesh);

  TrialFunction u(vh);
  TestFunction  v(vh);

  RealFunction f = 1;

  // Apply Dirichlet conditions on the entire boundary.
  Problem poisson(u, v);
  poisson = Integral(Grad(u), Grad(v))
          // - Integral(f, v)
          + DirichletBC(u, f).on(17);
  CG(poisson).solve();

  mesh.save("Poisson.mesh", IO::FileFormat::MEDIT);
  u.getSolution().save("Poisson.gf", IO::FileFormat::MFEM);

  // Save solution
  IO::XDMF xdmf("Poisson");
  xdmf.grid().setMesh(mesh).add("u", u.getSolution());
  xdmf.write();
  xdmf.close();

  return 0;
}
