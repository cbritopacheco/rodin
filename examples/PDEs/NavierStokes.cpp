/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [navierstokes-source]
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
  //! [navierstokes-mesh]
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
  //! [navierstokes-mesh]

  //! [navierstokes-fields]
  H1 velocitySpace(std::integral_constant<size_t, 2>{}, mesh, 2);
  H1 pressureSpace(std::integral_constant<size_t, 1>{}, mesh);
  P0g meanSpace(mesh);
  TrialFunction u(velocitySpace);
  TrialFunction p(pressureSpace);
  TrialFunction multiplier(meanSpace);
  TestFunction v(velocitySpace);
  TestFunction q(pressureSpace);
  TestFunction meanTest(meanSpace);
  constexpr Real viscosity =
    0.01; // Reynolds number 100 for unit length and peak lid speed.
  constexpr Real dt = 0.05;
  constexpr size_t steps = 200; // Final time 10, for the transient demonstration.
  constexpr Real rampTime = 1; // Ramp the lid over one convective time unit.
  GridFunction velocity(velocitySpace);
  GridFunction pressure(pressureSpace);
  velocity = Zero(2);
  pressure = 0.0;
  IO::XDMF xdmf("NavierStokes");
  xdmf.grid().setMesh(mesh).add("velocity", velocity).add("pressure", pressure);
  xdmf.write(0).flush();
  //! [navierstokes-fields]

  //! [navierstokes-step]
  for (size_t step = 1; step <= steps; ++step)
  {
    const Real time = step * dt;
    const Real ramp = std::min<Real>(1, time / rampTime);
    VectorFunction lidVelocity{ramp * 16 * Pow(F::x, 2) * Pow(1 - F::x, 2), Zero()};
    Problem flow(u, p, multiplier, v, q, meanTest);
    flow = (1 / dt) * Integral(u, v) + Integral(Dot(Mult(Jacobian(u), velocity), v)) +
      viscosity * Integral(Jacobian(u), Jacobian(v)) - Integral(p, Div(v)) +
      Integral(Div(u), q) + Integral(multiplier, q) + Integral(p, meanTest) -
      (1 / dt) * Integral(velocity, v) + DirichletBC(u, lidVelocity).on(lid) +
      DirichletBC(u, Zero(2)).on(walls);
    Solver::SparseLU solver(flow);
    solver.solve();
    if (!solver.success())
      return EXIT_FAILURE;
    velocity.setData(u.getSolution().getData());
    pressure.setData(p.getSolution().getData());
    xdmf.write(time).flush();
    std::cout << "Step " << step << " / " << steps << ", time " << time << std::endl;
  }
  //! [navierstokes-step]
  return 0;
}
//! [navierstokes-source]
