/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ActiveContraction.cpp
 * @brief Minimal monolithic active-contraction demo without local variable elimination.
 *
 * Unknowns at each time step are solved in one coupled system:
 *   - displacement u      (H1 vector)
 *   - active extension ec (P0g scalar, one dof over the mesh)
 *   - gamma               (P0g scalar)
 *   - beta                (P0g scalar)
 *
 * This keeps a monolithic block structure and avoids Schur/local elimination.
 * The model below is intentionally minimal and pedagogical (Hill–Maxwell-inspired),
 * designed to show the architecture direction requested by the user.
 */
#include <cstddef>

#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <Rodin/IO/XDMF.h>
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
  constexpr Real Lx = static_cast<Real>(nx - 1) / static_cast<Real>(ny - 1);

  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { nx, ny });
  mesh.scale(1.0 / static_cast<Real>(ny - 1));
  mesh.getConnectivity().compute(1, 2);

  // Clamp left boundary to remove rigid modes.
  constexpr Attribute leftBC = 1;
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
  }

  const size_t dim = mesh.getSpaceDimension();

  // ---- FE spaces ----------------------------------------------------------
  // u  : vector H1
  // ec, gamma, beta : P0g (single scalar dof over whole mesh)
  P1 Vh(mesh, dim);
  P0g Qh(mesh);

  // ---- state at time n ----------------------------------------------------
  GridFunction u_n(Vh);
  u_n.setName("Displacement");
  u_n = VectorFunction{ Zero(), Zero() };

  GridFunction ec_n(Qh);
  ec_n.setName("ec");
  ec_n = RealFunction(0.0);

  GridFunction gamma_n(Qh);
  gamma_n.setName("gamma");
  gamma_n = RealFunction(1.0);

  GridFunction beta_n(Qh);
  beta_n.setName("beta");
  beta_n = RealFunction(0.0);

  // ---- output -------------------------------------------------------------
  IO::XDMF xdmf("ActiveContraction");
  auto grid = xdmf.grid();
  grid.setMesh(mesh);
  grid.add(u_n);
  grid.add(ec_n);
  grid.add(gamma_n);
  grid.add(beta_n);

  // ---- model parameters (minimal Hill–Maxwell-inspired) -------------------
  const Real E  = 120.0;
  const Real nu = 0.35;
  const Real lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  const Real mu     = E / (2.0 * (1.0 + nu));

  const Real Es        = 40.0;
  const Real muParallel = 0.5;
  const Real n0         = 0.8;
  const Real k0         = 1.0;
  const Real sigma0     = 1.0;
  const Real alphaD     = 0.2;

  constexpr size_t nSteps = 20;
  const Real T = 1.0;
  const Real dt = T / static_cast<Real>(nSteps);
  const Real pi = Math::Constants::pi();

  // ---- time stepping ------------------------------------------------------
  for (size_t step = 1; step <= nSteps; ++step)
  {
    const Real time = dt * static_cast<Real>(step);

    // Unknowns at t^{n+1}
    TrialFunction u(Vh);
    TrialFunction ec(Qh);
    TrialFunction gamma(Qh);
    TrialFunction beta(Qh);

    // Tests
    TestFunction v(Vh);
    TestFunction w(Qh);
    TestFunction z(Qh);
    TestFunction r(Qh);

    // Smooth prescribed activation in [0, 1] with mild x-modulation.
    const auto activation =
      0.5 * (1.0 + Sin(2.0 * pi * time)) * (0.6 + 0.4 * F::x / Lx);

    // kinematic proxy for fiber strain (minimal demo): div(u)
    const auto e1D = Div(u);

    // Active strain rate (backward Euler)
    const auto edot = (ec - ec_n) / dt;

    // Positive regularized damping factor (minimal smooth variant)
    const auto damp = activation + alphaD * edot * edot;

    // Residual block 1: mechanics (linear elastic + active contribution)
    const auto Ru =
        Integral(mu * Jacobian(u), Jacobian(v))
      + Integral(lambda * Div(u), Div(v))
      + Integral(Es * (e1D - ec) * Div(v));

    // Residual block 2: active branch constraint (minimal HM-inspired form)
    const auto Rec =
      Integral((gamma * beta + muParallel * edot) - Es * (e1D - ec), w);

    // Residual block 3: gamma evolution (no elimination)
    const auto Rg =
      Integral((gamma - gamma_n) / dt - (n0 * k0 * activation - damp * gamma), z);

    // Residual block 4: beta evolution (no elimination)
    const auto Rb =
      Integral((beta - beta_n) / dt - (n0 * sigma0 * activation - damp * beta + edot * gamma), r);

    Problem monolithic(u, ec, gamma, beta, v, w, z, r);
    monolithic =
        Ru + Rec + Rg + Rb
      + DirichletBC(u, VectorFunction{ Zero(), Zero() }).on(leftBC);

    monolithic.assemble();
    SparseLU(monolithic).solve();

    // Commit n+1 -> n
    u_n = u.getSolution();
    ec_n = ec.getSolution();
    gamma_n = gamma.getSolution();
    beta_n = beta.getSolution();

    xdmf.write(time);
  }

  return 0;
}
