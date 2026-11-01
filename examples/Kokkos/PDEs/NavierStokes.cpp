/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>

#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Kokkos.h>

#ifdef KOKKOS_ENABLE_CUDA
using Backend = Rodin::Context::CUDA;
#else
using Backend = ::Kokkos::DefaultHostExecutionSpace;
#endif

// A global using directive on Rodin would make Kokkos ambiguous in the Kokkos headers.
using DeviceContext = Rodin::Context::Kokkos<Rodin::Context::Local, Backend>;

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  ::Kokkos::ScopeGuard runtime(argc, argv);

  constexpr size_t nodes = 33;
  constexpr Attribute lid = 2;
  constexpr Attribute walls = 3;
  constexpr Real boundaryTolerance = 1e-12; // Coordinate tolerance on this unit grid.
  Mesh host = Mesh().UniformGrid(Polytope::Type::Triangle, {nodes, nodes});
  host.scale(1.0 / (nodes - 1));
  host.getConnectivity().compute(1, 2);
  for (auto edge = host.getFace(); !edge.end(); ++edge)
  {
    if (!edge->isBoundary())
      continue;
    bool onTop = true;
    for (auto vertex = edge->getVertex(); !vertex.end(); ++vertex)
      onTop = onTop && std::abs(vertex->getCoordinates().y() - 1) < boundaryTolerance;
    host.setAttribute({1, edge->getIndex()}, onTop ? lid : walls);
  }
  const Mesh<DeviceContext> mesh(DeviceContext{}, host);

  // Taylor-Hood pair, with a global constant multiplier fixing the mean pressure.
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
  // Final time 20 by default, for the transient demonstration.
  const size_t steps = argc > 1 ? std::stoul(argv[1]) : 2000;
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

  // The problem and the solver outlive the time loop: the sparsity pattern is
  // the same at every step, so cuDSS keeps its symbolic analysis and only
  // refactorizes the matrix. With cuDSS the matrix stays in the memory of the
  // backend, where the solver reads it in place.
#ifdef RODIN_USE_CUDSS
  using System =
    Math::LinearSystem<Rodin::Kokkos::SparseMatrix<Real, Backend>, Math::Vector<Real>>;
  Problem<System, decltype(u), decltype(p), decltype(multiplier), decltype(v),
    decltype(q), decltype(meanTest)>
    flow(u, p, multiplier, v, q, meanTest);
  Solver::CuDSS solver(flow);
#else
  Problem flow(u, p, multiplier, v, q, meanTest);
  Solver::SparseLU solver(flow);
#endif
  std::chrono::duration<double> assembly{0}, solve{0};
  for (size_t step = 1; step <= steps; ++step)
  {
    const Real time = step * dt;
    const Real ramp = std::min<Real>(1, time / rampTime);
    const Real center = 0.5 - travelAmplitude * std::cos(2 * pi * time / travelPeriod);
    const Real speed = peakLidSpeed * ramp;
    const Real width = halfWidth;
    VectorFunction lidVelocity(
      KOKKOS_LAMBDA(const Rodin::Kokkos::Geometry::PointView& point)->Real {
        const Real offset = (point.x() - center) / width;
        const Real profile = offset * offset < 1 ? 1 - offset * offset : Real(0);
        return speed * profile * profile;
      },
      Real(0));
    flow = (1 / dt) * Integral(u, v) + Integral(Dot(Mult(Jacobian(u), velocity), v)) +
      viscosity * Integral(Jacobian(u), Jacobian(v)) - Integral(p, Div(v)) +
      Integral(Div(u), q) + Integral(multiplier, q) + Integral(p, meanTest) -
      (1 / dt) * Integral(velocity, v) + DirichletBC(u, lidVelocity).on(lid) +
      DirichletBC(u, Zero(2)).on(walls);
    auto t0 = std::chrono::steady_clock::now();
    flow.assemble();
    auto t1 = std::chrono::steady_clock::now();
#ifdef RODIN_USE_CUDSS
    solver.factorize(flow.getLinearSystem());
    if (!solver.success())
      return EXIT_FAILURE;
#endif
    solver.solve();
    if (!solver.success())
      return EXIT_FAILURE;
    assembly += t1 - t0;
    solve += std::chrono::steady_clock::now() - t1;
    velocity.setData(u.getSolution().getData());
    pressure.setData(p.getSolution().getData());
    if (step % outputEvery == 0)
      xdmf.write(time).flush();
    std::cout << "Step " << step << " / " << steps << ", time " << time << ", lid center "
              << center << std::endl;
  }
  std::cout << "Assembly " << assembly.count() << " s, solve " << solve.count()
            << " s in " << steps << " steps" << std::endl;
  return 0;
}
