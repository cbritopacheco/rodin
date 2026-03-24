/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file CantileverBeam.cpp
 * @brief Quasi-static hyperelastic cantilever beam under time-varying traction.
 *
 * A rectangular 2D beam (4 × 1) is clamped on the left edge and subject to
 * a downward traction on the right edge.  The traction ramps up linearly
 * over the first half of the time steps and is then suddenly removed, so the
 * NeoHookean material springs back to its undeformed configuration.
 *
 * Output is written to XDMF for visualization in ParaView (apply "Warp by
 * Vector" to see the deformed shape).
 */
#include <cstddef>

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
  // UniformGrid({nx, ny}) produces nx × ny vertices on [0, nx-1] × [0, ny-1].
  // After scaling by 1/(ny-1) the domain becomes [0, Lx] × [0, 1].
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

  // ---- solution -----------------------------------------------------------
  GridFunction u(Vh);
  u.setName("Displacement");
  u = VectorFunction{ Zero(), Zero() };

  // ---- XDMF output --------------------------------------------------------
  IO::XDMF xdmf("CantileverBeam");
  auto grid = xdmf.grid();
  grid.setMesh(mesh);
  grid.add(u);

  // ---- quasi-static time stepping -----------------------------------------
  // The traction ramps up linearly then drops to zero (sudden release).
  constexpr size_t nSteps = 20;
  constexpr Real maxTraction = -5.0;  // downward

  TrialFunction du(Vh);
  TestFunction  v(Vh);
  auto zero = VectorFunction{ Zero(), Zero() };

  for (size_t step = 0; step <= nSteps; ++step)
  {
    Real ty = 0;
    if (step <= nSteps / 2)
      ty = maxTraction * static_cast<Real>(step) / static_cast<Real>(nSteps / 2);

    auto traction = VectorFunction{ Zero(), RealFunction(ty) };

    Solid::MaterialTangent tangent(law, du, v);
    tangent.setLinearizationPoint(u.getData());

    Solid::InternalForce residual(law, v);
    residual.setLinearizationPoint(u.getData());

    // Newton linearization:  K δu = -F_int(u) + F_ext
    Problem newton(du, v);
    newton = tangent
           + residual
           - BoundaryIntegral(traction, v).over(rightBC)
           + DirichletBC(du, zero).on(leftBC);

    SparseLU linearSolver(newton);
    NewtonSolver solver(linearSolver);
    solver.setMaxIterations(50)
          .setAbsoluteTolerance(1e-10)
          .setRelativeTolerance(1e-8);

    solver.solve(u.getData());

    xdmf.write(static_cast<Real>(step));
  }

  return 0;
}
