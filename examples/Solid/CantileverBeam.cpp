/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file DampedTransientCantileverBeam.cpp
 * @brief Transient hyperelastic cantilever beam with inertia, damping, and release.
 *
 * A rectangular 2D beam (4 × 1) is clamped on the left edge and subject to
 * a downward traction on the right edge. The traction ramps up linearly over
 * the first half of the simulation and is then suddenly removed.
 *
 * The beam is modeled as a compressible NeoHookean solid with inertia and
 * simple mass-proportional viscous damping:
 *
 *   rho * u_tt + c * u_t - div(P(F(u))) = 0
 *
 * Time integration is performed with the Newmark-beta method using the
 * average-acceleration choice beta = 1/4, gamma = 1/2.
 *
 * Output is written to XDMF for visualization in ParaView.
 */
#include "Rodin/Math/Vector.h"
#include <cstddef>

#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <Rodin/Solid.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Solver/SparseLU.h>
#include <Rodin/Solver/CG.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

int main(int, char**)
{
  // ---- geometry -----------------------------------------------------------
  constexpr size_t nx = 33;
  constexpr size_t ny = 9;
  constexpr Real Lx = static_cast<Real>(nx - 1) / static_cast<Real>(ny - 1);

  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { nx, ny });
  mesh.scale(1.0 / static_cast<Real>(ny - 1));
  mesh.getConnectivity().compute(1, 2);

  // ---- label boundary edges -----------------------------------------------
  constexpr Attribute leftBC  = 1;
  constexpr Attribute rightBC = 2;
  constexpr Real eps = 1e-10;

  for (auto it = mesh.getBoundary(); !it.end(); ++it)
  {
    const auto& verts = it->getVertices();
    const size_t nv = verts.size();

    Real xSum = 0;
    for (size_t i = 0; i < nv; ++i)
      xSum += mesh.getVertexCoordinates(verts[i])(0);
    const Real xMid = xSum / static_cast<Real>(nv);

    if (xMid < eps)
      mesh.setAttribute({ 1, it->getIndex() }, leftBC);
    else if (xMid > Lx - eps)
      mesh.setAttribute({ 1, it->getIndex() }, rightBC);
  }

  // ---- finite-element space -----------------------------------------------
  const size_t dim = mesh.getSpaceDimension();
  P1 Vh(mesh, dim);

  // ---- material -----------------------------------------------------------
  const Real E  = 200.0;
  const Real nu = 0.3;
  const Real lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  const Real mu     = E / (2.0 * (1.0 + nu));
  Solid::NeoHookean law(lambda, mu);

  // ---- dynamics parameters ------------------------------------------------
  const Real rho = 1.0;
  const Real c   = 0.5;   // mass-proportional viscous damping coefficient
  const Real dt  = 1e-5;
  const size_t nSteps = 200;

  // Newmark average acceleration
  const Real beta  = 0.25;
  const Real gamma = 0.5;

  // Effective coefficients
  const Real aMass = rho / (beta * dt * dt);
  const Real aDamp = c * gamma / (beta * dt);

  // ---- fields -------------------------------------------------------------
  GridFunction u(Vh);      // displacement
  GridFunction vel(Vh);    // velocity
  GridFunction acc(Vh);    // acceleration

  u.setName("Displacement");
  vel.setName("Velocity");
  acc.setName("Acceleration");

  auto zero = VectorFunction{ Zero(), Zero() };
  u   = zero;
  vel = zero;
  acc = zero;

  GridFunction uPred(Vh);
  GridFunction vPred(Vh);
  GridFunction aNew(Vh);
  GridFunction rhsDamp(Vh);

  uPred.setName("PredictedDisplacement");
  vPred.setName("PredictedVelocity");
  aNew.setName("AccelerationNew");
  rhsDamp.setName("DampingRHS");

  // ---- output -------------------------------------------------------------
  IO::XDMF xdmf("CantileverBeam");
  auto grid = xdmf.grid();
  grid.setMesh(mesh);
  grid.add(u);
  grid.add(vel);
  grid.add(acc);

  xdmf.write(0.0);

  // ---- time loop ----------------------------------------------------------
  constexpr Real maxTraction = -0.1; // downward

  TrialFunction du(Vh);
  TestFunction  w(Vh);

  for (size_t step = 1; step <= nSteps; ++step)
  {
    const Real t = step * dt;

    Real ty = 0.0;
    if (step <= nSteps / 2)
      ty = maxTraction * static_cast<Real>(step) / static_cast<Real>(nSteps / 2);

    auto traction = VectorFunction{ Zero(), RealFunction(ty) };

    // ---- Newmark predictors -----------------------------------------------
    uPred = u + dt * vel + (dt * dt * (0.5 - beta)) * acc;
    vPred = vel + (dt * (1.0 - gamma)) * acc;

    // Damping RHS term:
    // c * M * (vPred - gamma/(beta dt) * uPred)
    rhsDamp = vPred - (gamma / (beta * dt)) * uPred;

    // Use predictor as initial guess for u^{n+1}
    u = uPred;

    // ---- nonlinear solid operators ----------------------------------------
    Solid::MaterialTangent tangent(law, du, w);
    tangent.setLinearizationPoint(u.getData());

    Solid::InternalForce internal(law, w);
    internal.setLinearizationPoint(u.getData());

    // Effective nonlinear problem:
    //
    //   F_int(u_{n+1})
    // + aMass * M u_{n+1}
    // + aDamp * M u_{n+1}
    // = F_ext(t_{n+1})
    // + aMass * M uPred
    // - c * M rhsDamp
    //
    // where rhsDamp = gamma/(beta dt) * uPred - vPred with the chosen sign
    // convention below:
    //
    //   rhsDamp = vPred - gamma/(beta dt) * uPred
    //
    // so the right-hand side contribution is:
    //
    //   + c * M rhsDamp
    //
    Problem newton(du, w);
    newton = tangent
           + aMass * Integral(du, w)
           + aDamp * Integral(du, w)
           + internal
           - BoundaryIntegral(traction, w).over(rightBC)
           - aMass * Integral(uPred, w)
           - c * Integral(rhsDamp, w)
           + DirichletBC(du, zero).on(leftBC);

    CG linearSolver(newton);
    NewtonSolver solver(linearSolver);
    solver.setMaxIterations(50)
          .setAbsoluteTolerance(1e-10)
          .setRelativeTolerance(1e-8);
    Math::Vector<Real> x;
    x.resize(u.getSize());
    x.setZero();
    solver.solve(x);
    u.setData(x);
    std::cout << "Step " << step << ", time " << t << std::endl;

    // ---- Newmark correctors -----------------------------------------------
    aNew = (1.0 / (beta * dt * dt)) * (u - uPred);
    vel  = vPred + (gamma * dt) * aNew;
    acc  = aNew;

    xdmf.write(t).flush();
  }

  xdmf.close();
  return 0;
}
