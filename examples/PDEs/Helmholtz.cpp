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
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

static constexpr Real waveLength = 8;
const Real waveNumber = 2 * Math::Constants::pi() / waveLength;

int main(int, char**)
{
  // Build a mesh
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 16, 16 });
  mesh.getConnectivity().compute(1, 2);

  // Functions
  P1<Complex> vh(mesh);

  TrialFunction u(vh);
  TestFunction  v(vh);

  ComplexFunction f = [](const Point& p) { return Complex(1, 1); };
  Problem helmholtz(u, v);
  helmholtz = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());

  helmholtz.assemble();

  CG(helmholtz).solve();

  P1<Real> rh(mesh);
  GridFunction uRe(rh);
  GridFunction uIm(rh);

  uRe = Re(u.getSolution());
  uIm = Im(u.getSolution());

  // Save solution
  uRe.save("uRe.gf");
  uIm.save("uIm.gf");
  mesh.save("Grid.mesh");

  return 0;
}

