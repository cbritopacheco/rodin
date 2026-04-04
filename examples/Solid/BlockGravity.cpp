/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file BlockGravity.cpp
 * @brief Hyperelastic block under time-varying gravity.
 *
 * A unit square block is fixed on its bottom edge and subject to a
 * downward body force (gravity).  The load is ramped up, held steady,
 * and then gradually removed so the NeoHookean material returns to its
 * undeformed configuration, demonstrating the internal-force–driven
 * elastic spring-back.
 *
 * Output is written to XDMF for visualization in ParaView (apply "Warp
 * by Vector" to see the deformed shape).
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
  constexpr size_t nc = 17;  // 17 × 17 vertices → 16 × 16 cells
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { nc, nc });
  mesh.scale(1.0 / static_cast<Real>(nc - 1));  // unit square [0,1]²
  mesh.getConnectivity().compute(1, 2);

  // ---- label bottom boundary ----------------------------------------------
  constexpr Attribute bottomBC = 1;
  constexpr Real eps = 1e-10;

  for (auto it = mesh.getBoundary(); !it.end(); ++it)
  {
    const auto& verts = it->getVertices();
    const size_t nv = verts.size();

    Real ySum = 0;
    for (size_t i = 0; i < nv; ++i)
      ySum += mesh.getVertexCoordinates(verts[i])(1);
    const Real yMid = ySum / static_cast<Real>(nv);

    if (yMid < eps)
      mesh.setAttribute({ 1, it->getIndex() }, bottomBC);
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
  IO::XDMF xdmf("BlockGravity");
  auto grid = xdmf.grid();
  grid.setMesh(mesh);
  grid.add(u);

  // ---- quasi-static time stepping -----------------------------------------
  // Schedule:  steps 0–5   ramp gravity up
  //            steps 5–10  hold at maximum
  //            steps 10–15 ramp gravity down
  constexpr size_t nSteps = 15;
  constexpr Real maxGravity = -10.0;  // downward body force density

  TrialFunction du(Vh);
  TestFunction  v(Vh);
  auto zero = VectorFunction{ Zero(), Zero() };

  for (size_t step = 0; step <= nSteps; ++step)
  {
    Real gy = 0;
    if (step <= 5)
      gy = maxGravity * static_cast<Real>(step) / 5.0;
    else if (step <= 10)
      gy = maxGravity;
    else
      gy = maxGravity * static_cast<Real>(nSteps - step) / 5.0;

    auto bodyForce = VectorFunction{ Zero(), RealFunction(gy) };

    Solid::MaterialTangent tangent(law, du, v);
    tangent.setDisplacement(u);

    Solid::InternalForce residual(law, v);
    residual.setDisplacement(u);

    // Newton linearization:  K δu = -F_int(u) + F_body
    Problem newton(du, v);
    newton = tangent
           + residual
           - Integral(bodyForce, v)
           + DirichletBC(du, zero).on(bottomBC);

    SparseLU linearSolver(newton);
    NewtonSolver solver(linearSolver);
    solver.setMaxIterations(50)
          .setAbsoluteTolerance(1e-10)
          .setRelativeTolerance(1e-8);

    solver.solve(u);

    xdmf.write(static_cast<Real>(step));
  }

  return 0;
}
