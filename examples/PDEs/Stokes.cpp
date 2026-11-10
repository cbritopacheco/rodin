/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [stokes-source]
#include <cmath>
#include <cstdlib>
#include <iostream>

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
  //! [stokes-mesh]
  constexpr size_t nodes = 17;
  constexpr Attribute lid = 2;
  constexpr Attribute walls = 3;
  constexpr Real boundaryTolerance = 1e-12; // Coordinate tolerance on this unit grid.
  Mesh mesh = Mesh().UniformGrid(Polytope::Type::Triangle, {nodes, nodes});
  mesh.scale(1.0 / (nodes - 1));
  mesh.getConnectivity().compute(1, 2);
  for (auto edge = mesh.getFace(); !edge.end(); ++edge)
  {
    if (!edge->isBoundary())
      continue;
    bool onTop = true;
    for (auto vertex = edge->getVertex(); !vertex.end(); ++vertex)
      onTop = onTop && std::abs(vertex->getCoordinates().y() - 1) < boundaryTolerance;
    mesh.setAttribute({1, edge->getIndex()}, onTop ? lid : walls);
  }
  //! [stokes-mesh]

  //! [stokes-form]
  H1 velocitySpace(std::integral_constant<size_t, 2>{}, mesh, 2);
  H1 pressureSpace(std::integral_constant<size_t, 1>{}, mesh);
  P0g meanSpace(mesh);
  TrialFunction u(velocitySpace);
  TrialFunction p(pressureSpace);
  TrialFunction multiplier(meanSpace);
  TestFunction v(velocitySpace);
  TestFunction q(pressureSpace);
  TestFunction meanTest(meanSpace);
  constexpr Real viscosity = 1;
  VectorFunction lidVelocity{16 * Pow(F::x, 2) * Pow(1 - F::x, 2), Zero()};
  Problem stokes(u, p, multiplier, v, q, meanTest);
  stokes = viscosity * Integral(Jacobian(u), Jacobian(v)) - Integral(p, Div(v)) +
    Integral(Div(u), q) + Integral(multiplier, q) + Integral(p, meanTest) +
    DirichletBC(u, lidVelocity).on(lid) + DirichletBC(u, Zero(2)).on(walls);
  //! [stokes-form]

  //! [stokes-solve]
  Solver::SparseLU solver(stokes);
  solver.solve();
  if (!solver.success())
    return EXIT_FAILURE;
  const auto& system = stokes.getLinearSystem();
  const Real residual =
    (system.getOperator() * system.getSolution() - system.getVector()).norm();
  std::cout << "Algebraic residual: " << residual << std::endl;
  IO::XDMF xdmf("Stokes");
  xdmf.grid()
    .setMesh(mesh)
    .add("velocity", u.getSolution())
    .add("pressure", p.getSolution());
  xdmf.write().flush();
  //! [stokes-solve]
  return 0;
}
//! [stokes-source]
