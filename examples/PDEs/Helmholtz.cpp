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

static constexpr Real waveLength = 0.5;
static constexpr Real pi = Math::Constants::pi();
const Real waveNumber = 2 * Math::Constants::pi() / waveLength;

int main(int, char**)
{
  // Build a mesh
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 64, 64 });
  mesh.getConnectivity().compute(1, 2);
  mesh.scale(1.0 / 63.0);

  // Functions
  P1<Complex> vh(mesh);

  GridFunction gf(vh);
  gf = [&](const Point& p)
  {
    return
      Complex(p.x(), p.x() * p.y());// * cos(waveNumber * p.x()) * cos(waveNumber * p.y());
  };

  P1<Real> rh(mesh);
  GridFunction gfRe(rh);
  GridFunction gfIm(rh);
  gfRe = Re(gf);
  gfIm = Im(gf);
  gfRe.save("gfRe.gf");
  gfIm.save("gfIm.gf");
  mesh.save("Grid.mesh");

  ComplexFunction f =
    [&](const Point& p)
    {
      return Complex(1, 1) * waveNumber * waveNumber * cos(waveNumber * p.x()
          ) * cos(waveNumber * p.y());
    };

  TrialFunction u(vh);
  TestFunction  v(vh);

  Problem helmholtz(u, v);
  helmholtz = Integral(Grad(u), Grad(v))
            - waveNumber * waveNumber * Integral(u, v)
            - Integral(f, v)
            ;

  CG(helmholtz).solve();

  GridFunction uRe(rh);
  GridFunction uIm(rh);

  uRe = Re(u.getSolution());
  uIm = Im(u.getSolution());

  // Save solution
  uRe.save("uRe.gf");
  uIm.save("uIm.gf");

  GridFunction miaow(rh);
  miaow = [&](const Point& p)
  {
    return
      (Complex(p.x(), 1) * cos(waveNumber * p.x()) * cos(waveNumber * p.y())).real();
  };

  miaow.save("exRe.gf");

  miaow = [&](const Point& p)
  {
    return
      (Complex(p.x(), 1) * cos(waveNumber * p.x()) * cos(waveNumber * p.y())).imag();
  };

  miaow.save("exIm.gf");

  mesh.save("Grid.mesh");

  return 0;
}

