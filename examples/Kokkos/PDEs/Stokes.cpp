/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cmath>
#include <cstdlib>
#include <iostream>

#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Kokkos.h>

#if defined(KOKKOS_ENABLE_CUDA)
using Backend = Rodin::Context::CUDA;

#elif defined(KOKKOS_ENABLE_HIP)
using Backend = Rodin::Context::HIP;
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

  constexpr size_t nodes = 17;
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
  P1 pressureSpace(mesh);
  P0g meanSpace(mesh);
  TrialFunction u(velocitySpace);
  TrialFunction p(pressureSpace);
  TrialFunction multiplier(meanSpace);
  TestFunction v(velocitySpace);
  TestFunction q(pressureSpace);
  TestFunction meanTest(meanSpace);
  constexpr Real viscosity = 1;
  VectorFunction lidVelocity(
    KOKKOS_LAMBDA(const Rodin::Kokkos::Geometry::PointView& x)->Real {
      return 16 * x.x() * x.x() * (1 - x.x()) * (1 - x.x());
    },
    Real(0));
  Problem stokes(u, p, multiplier, v, q, meanTest);
  stokes = viscosity * Integral(Jacobian(u), Jacobian(v)) - Integral(p, Div(v)) +
    Integral(Div(u), q) + Integral(multiplier, q) + Integral(p, meanTest) +
    DirichletBC(u, lidVelocity).on(lid) + DirichletBC(u, Zero(2)).on(walls);

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
  return 0;
}
