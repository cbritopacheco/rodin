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
  constexpr size_t nodes = 33;
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
  constexpr Real peakLidSpeed = 2; // Doubled driving speed for the stronger transient.
  constexpr Real viscosity = 0.01; // Reynolds number 200 for unit cavity width.
  constexpr Real dt = 0.01;
  constexpr size_t steps = 2000; // Final time 20, for the transient demonstration.
  constexpr size_t outputEvery = 5; // Export every 0.05 time units for visualization.
  constexpr Real rampTime = 1; // Ramp the lid over one simulation time unit.
  constexpr Real travelAmplitude = 0.22; // Center travels between x=0.28 and x=0.72.
  constexpr Real travelPeriod = 8; // One horizontal round trip in simulation time.
  constexpr Real halfWidth = 0.25; // Compact support stays clear of both corners.
  const Real pi = Math::Constants::pi();
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
    const Real center = 0.5 - travelAmplitude * std::cos(2 * pi * time / travelPeriod);
    RealFunction lidSpeed = [=](const Point& point) {
      const Real offset = (point.x() - center) / halfWidth;
      const Real profile = std::max<Real>(0, 1 - offset * offset);
      return peakLidSpeed * ramp * profile * profile;
    };
    VectorFunction lidVelocity{lidSpeed, Zero()};
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
    if (step % outputEvery == 0)
      xdmf.write(time).flush();
    std::cout << "Step " << step << " / " << steps << ", time " << time << ", lid center "
              << center << std::endl;
  }
  //! [navierstokes-step]
  return 0;
}
//! [navierstokes-source]
