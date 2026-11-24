/*
 *          Copyright Jean-Baptiste KECK 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cmath>
#include <cstdlib>

#include <Rodin/Solver.h>
#include <Rodin/Kokkos/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Kokkos.h>

#if defined(KOKKOS_ENABLE_CUDA)
using Backend = Rodin::Context::CUDA;

#elif defined(KOKKOS_ENABLE_HIP)
using Backend = Rodin::Context::HIP;
#elif defined(KOKKOS_ENABLE_SYCL)
using Backend = Rodin::Context::SYCL;
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

  constexpr Real spacing = 0.1;
  constexpr Real length = 4;
  constexpr Real boundaryTolerance = 1e-12; // Coordinate tolerance for this grid.
  constexpr Attribute fixed = 2;
  constexpr Attribute loaded = 3;
  constexpr Attribute free = 4;

  Mesh host = Mesh().UniformGrid(Polytope::Type::Triangle, {41, 11});
  host.scale(spacing);
  host.getConnectivity().compute(1, 2);

  for (auto face = host.getFace(); !face.end(); ++face)
  {
    if (!face->isBoundary())
      continue;
    bool onLeft = true;
    bool onRight = true;
    for (auto vertex = face->getVertex(); !vertex.end(); ++vertex)
    {
      const Real x = vertex->getCoordinates().x();
      onLeft = onLeft && std::abs(x) < boundaryTolerance;
      onRight = onRight && std::abs(x - length) < boundaryTolerance;
    }
    host.setAttribute({1, face->getIndex()}, onLeft ? fixed : onRight ? loaded : free);
  }
  const Mesh<DeviceContext> mesh(DeviceContext{}, host);

  P1 vh(mesh, mesh.getSpaceDimension());
  TrialFunction u(vh);
  TestFunction v(vh);

  constexpr Real lambda = 0.5769;
  constexpr Real mu = 0.3846;
  VectorFunction traction{0, -1e-4};

  Problem elasticity(u, v);
  elasticity = LinearElasticityIntegral(u, v)(lambda, mu) -
    BoundaryIntegral(traction, v).over(loaded) + DirichletBC(u, Zero(2)).on(fixed);

  Solver::CG solver(elasticity);
  solver.setTolerance(1e-10).solve();
  if (!solver.success())
    return EXIT_FAILURE;

  IO::XDMF xdmf("Elasticity");
  xdmf.grid().setMesh(mesh).add("displacement", u.getSolution());
  xdmf.write().flush();
  return 0;
}
