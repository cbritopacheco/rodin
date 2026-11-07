/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cstdlib>

#include <Rodin/Types.h>
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
  using namespace Rodin::Solver;
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  ::Kokkos::ScopeGuard runtime(argc, argv);

  constexpr size_t nodes = 16;
  Mesh host = Mesh().UniformGrid(Polytope::Type::Triangle, {nodes, nodes});
  host.scale(1.0 / (nodes - 1));
  host.getConnectivity().compute(1, 2);
  const Mesh<DeviceContext> mesh(DeviceContext{}, host);

  P1 vh(mesh);
  TrialFunction u(vh);
  TestFunction v(vh);
  RealFunction f = 1;

  Problem poisson(u, v);
  poisson = Integral(Grad(u), Grad(v)) - Integral(f, v) + DirichletBC(u, Zero());

  CG solver(poisson);
  solver.solve();
  if (!solver.success())
    return EXIT_FAILURE;

  IO::XDMF xdmf("Poisson");
  xdmf.grid().setMesh(mesh).add("u", u.getSolution());
  xdmf.write().flush();
  return 0;
}
