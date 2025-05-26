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
  // mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 32, 32 });
  mesh.load("SphereTetra.mesh", IO::FileFormat::MEDIT);
  // mesh.scale(1.0 / 31);
  // mesh.getConnectivity().compute(1, 2);
  mesh.save("Poisson.geo", IO::FileFormat::ENSIGHT6);
  mesh.save("Poisson.mesh", IO::FileFormat::MFEM);

  // Functions
  P1 vh(mesh, 3);

  GridFunction gf(vh);
  gf = [](Math::Vector<Real>& x, const Point& p)
  {
    x.resize(3);
    // x[0] = -p.x() * p.z();
    // x[1] = -p.y() * p.z();
    // x[2] = 1 - p.z() * p.z();
    x[0] = -p.y();
    x[1] = p.x();
    x[2] = 0.0;
  };
  gf.save("Poisson.vct", IO::FileFormat::ENSIGHT6);
  gf.save("Poisson.gf", IO::FileFormat::MFEM);

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
