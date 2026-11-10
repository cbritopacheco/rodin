/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [helmholtz-source]
#include <cmath>
#include <cstdlib>

#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <Rodin/IO/XDMF.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main()
{
  //! [helmholtz-mesh]
  constexpr size_t nodes = 97;
  const Real pi = Math::Constants::pi();
  const Real waveNumber = 4 * pi;
  Mesh mesh = Mesh().UniformGrid(Polytope::Type::Triangle, {nodes, nodes});
  mesh.scale(1.0 / (nodes - 1));
  mesh.getConnectivity().compute(1, 2);
  //! [helmholtz-mesh]

  //! [helmholtz-form]
  P1<Complex> vh(mesh);
  TrialFunction u(vh);
  TestFunction v(vh);
  ComplexFunction exact = [=](const Point& point) {
    return std::exp(Complex(0, waveNumber * point.x())) * std::sin(pi * point.y());
  };
  auto source = pi * pi * exact;

  Problem helmholtz(u, v);
  helmholtz = Integral(Grad(u), Grad(v)) - waveNumber * waveNumber * Integral(u, v) -
    Integral(source, v) + DirichletBC(u, exact);
  //! [helmholtz-form]

  //! [helmholtz-solve]
  Solver::SparseLU solver(helmholtz);
  solver.solve();
  if (!solver.success())
    return EXIT_FAILURE;
  //! [helmholtz-solve]

  //! [helmholtz-export]
  P1 realSpace(mesh);
  GridFunction realPart(realSpace);
  GridFunction imaginaryPart(realSpace);
  realPart = Re(u.getSolution());
  imaginaryPart = Im(u.getSolution());
  IO::XDMF xdmf("Helmholtz");
  xdmf.grid().setMesh(mesh).add("real", realPart).add("imaginary", imaginaryPart);
  xdmf.write().flush();
  //! [helmholtz-export]
  return 0;
}
//! [helmholtz-source]
