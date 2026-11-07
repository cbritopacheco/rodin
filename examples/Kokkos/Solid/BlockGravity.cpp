/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file BlockGravity.cpp
 * @brief Hyperelastic block under time-varying gravity, assembled on the device.
 *
 * Device port of examples/Solid/BlockGravity.cpp. The Neo-Hookean tangent and
 * residual are assembled by a device kernel at each Newton iteration.
 */
#include <cstddef>
#include <cstdlib>

#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <Rodin/Solid.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Solver/SparseLU.h>
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

  constexpr size_t nc = 17;
  Mesh host = Mesh().UniformGrid(Polytope::Type::Triangle, {nc, nc});
  host.scale(1.0 / static_cast<Real>(nc - 1));
  host.getConnectivity().compute(1, 2);

  constexpr Attribute bottomBC = 1;
  constexpr Real eps = 1e-10;
  for (auto it = host.getBoundary(); !it.end(); ++it)
  {
    const auto& verts = it->getVertices();
    Real ySum = 0;
    for (size_t i = 0; i < verts.size(); ++i)
      ySum += host.getVertexCoordinates(verts[i])(1);
    if (ySum / static_cast<Real>(verts.size()) < eps)
      host.setAttribute({1, it->getIndex()}, bottomBC);
  }
  const Mesh<DeviceContext> mesh(DeviceContext{}, host);

  P1 vh(mesh, mesh.getSpaceDimension());

  const Real E = 200.0;
  const Real nu = 0.3;
  const Real lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  const Real mu = E / (2.0 * (1.0 + nu));
  Solid::NeoHookean law(lambda, mu);

  GridFunction u(vh);
  u.getData().setZero();

  IO::XDMF xdmf("BlockGravity");
  xdmf.grid().setMesh(mesh).add("Displacement", u);

  constexpr size_t nSteps = 15;
  constexpr Real maxGravity = -10.0;

  TrialFunction du(vh);
  TestFunction v(vh);

  for (size_t step = 0; step <= nSteps; ++step)
  {
    Real gy = 0;
    if (step <= 5)
      gy = maxGravity * static_cast<Real>(step) / 5.0;
    else if (step <= 10)
      gy = maxGravity;
    else
      gy = maxGravity * static_cast<Real>(nSteps - step) / 5.0;

    auto ivw = Solid::InternalVirtualWork(law, u);

    Problem newton(du, v);
    newton = ivw(du, v) - Integral(VectorFunction{0, gy}, v) +
      DirichletBC(du, Zero(2)).on(bottomBC);

    Solver::SparseLU linearSolver(newton);
    Solver::NewtonSolver solver(linearSolver);
    solver.setMaxIterations(50).setAbsoluteTolerance(1e-10).setRelativeTolerance(1e-8);
    solver.solve(u);
    if (!solver.converged())
      return EXIT_FAILURE;

    xdmf.write(static_cast<Real>(step)).flush();
  }

  return 0;
}
