/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ActiveContraction.cpp
 * @brief Quasi-static active-contraction demo using FiberActiveI4Contraction.
 *
 * A 2D rectangular slab (aspect ratio 4:1) is clamped on the left edge.
 * The material is modeled as a NeoHookean solid wrapped with a
 * FiberActiveI4Contraction decorator that adds a fiber-aligned active stress
 * contribution.  Fiber direction is constant and aligned with the x-axis.
 *
 * The scalar activation level is ramped smoothly from 0 to 1 and back to 0
 * over nSteps quasi-static time steps.  At each step a Newton–Raphson solve
 * is performed to find the equilibrium displacement, driving an end-to-end
 * demonstration of the constitutive law, integrators, and solver chain.
 *
 * Output is written to XDMF for visualization (apply "Warp by Vector" in
 * ParaView to see the deformed shape).
 */
#include <cstddef>
#include <cmath>

#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <Rodin/Solid.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Solver/SparseLU.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

int main(int, char**)
{
  // ---- geometry -----------------------------------------------------------
  constexpr size_t nx = 65;
  constexpr size_t ny = 17;

  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { nx, ny });
  mesh.scale(1.0 / static_cast<Real>(ny - 1));
  mesh.getConnectivity().compute(1, 2);

  // Clamp left boundary.
  constexpr Attribute leftBC = 1;
  constexpr Real eps = 1e-10;
  for (auto it = mesh.getBoundary(); !it.end(); ++it)
  {
    const auto& verts = it->getVertices();
    const size_t nv = verts.size();

    Real xSum = 0;
    for (size_t i = 0; i < nv; ++i)
      xSum += mesh.getVertexCoordinates(verts[i])(0);

    if (xSum / static_cast<Real>(nv) < eps)
      mesh.setAttribute({ 1, it->getIndex() }, leftBC);
  }

  const size_t dim = mesh.getSpaceDimension();

  // ---- FE space -----------------------------------------------------------
  P1 Vh(mesh, dim);

  // ---- constitutive law ---------------------------------------------------
  const Real E      = 120.0;
  const Real nu     = 0.35;
  const Real lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  const Real mu     = E / (2.0 * (1.0 + nu));

  Solid::NeoHookean passive(lambda, mu);

  // Active tension scale and reference fiber stretch
  const Real activeTensionScale    = 10.0;
  const Real referenceFiberStretch = 1.0;
  Solid::FiberActiveI4Contraction<Solid::NeoHookean> law(
      passive, activeTensionScale, referenceFiberStretch);

  // Fiber direction: constant, aligned with the x-axis
  Math::SpatialVector<Real> fiberDir(static_cast<std::uint8_t>(dim));
  fiberDir.setZero();
  fiberDir(0) = 1.0;

  // ---- solution -----------------------------------------------------------
  GridFunction u(Vh);
  u.setName("Displacement");
  u = VectorFunction{ Zero(), Zero() };

  // ---- output -------------------------------------------------------------
  IO::XDMF xdmf("ActiveContraction");
  auto grid = xdmf.grid();
  grid.setMesh(mesh);
  grid.add(u);
  xdmf.write(0.0);

  // ---- Newton trial / test functions (reused across time steps) -----------
  TrialFunction du(Vh);
  TestFunction  v(Vh);
  auto zero = VectorFunction{ Zero(), Zero() };

  // ---- quasi-static time stepping -----------------------------------------
  constexpr size_t nSteps = 20;
  const Real pi = std::acos(static_cast<Real>(-1));

  for (size_t step = 1; step <= nSteps; ++step)
  {
    const Real t = static_cast<Real>(step) / static_cast<Real>(nSteps);

    // Activation: smooth ramp up then down over [0, 1]
    const Real activation = 0.5 * (1.0 - std::cos(pi * t));

    // Input: inject fiber direction and activation at every quadrature point
    auto input = [&](Solid::ConstitutivePoint& cp)
    {
      cp.set<Solid::Tags::FiberDirection>(fiberDir);
      cp.set<Solid::Tags::Activation>(activation);
    };

    // Nonlinear integrators at current linearization point u
    Solid::MaterialTangent tangent(law, du, v);
    tangent.setDisplacement(u);
    tangent.setInput(input);

    Solid::InternalForce internal(law, v);
    internal.setDisplacement(u);
    internal.setInput(input);

    // Newton problem:  K(u) du = -R(u)
    Problem newton(du, v);
    newton = tangent
           + internal
           + DirichletBC(du, zero).on(leftBC);

    SparseLU linearSolver(newton);
    NewtonSolver solver(linearSolver);
    solver.setMaxIterations(50)
          .setAbsoluteTolerance(1e-10)
          .setRelativeTolerance(1e-8);
    solver.solve(u);

    xdmf.write(t);
  }

  xdmf.close();
  return 0;
}
