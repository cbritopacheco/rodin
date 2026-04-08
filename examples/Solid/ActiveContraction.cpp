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
 * Unknowns at each time step are solved in a coupled system:
 *   - displacement u      (H1 vector)
 *   - active extension ec (P0g scalar, one DOF over the mesh)
 *   - gamma               (P0g scalar)
 *   - beta                (P0g scalar)
 *
 * This keeps the full block structure and avoids Schur/local elimination.
 * The model below is intentionally minimal and pedagogical (Hill–Maxwell-inspired),
 * designed to show the architecture direction.
 *
 * Because ec, gamma, beta are P0g (one DOF each, i.e. scalar numbers), and
 * the coupling terms include nonlinear products of unknowns (gamma * beta,
 * damp * gamma, edot * gamma), the system is solved with a staggered
 * fixed-point iteration at each time step:
 *
 *   1. Solve the linear mechanics for u given current ec
 *      (linear elastic + active coupling).
 *   2. Update ec, gamma, beta from their backward-Euler ODE equations
 *      given the current displacement field u.
 *   3. Repeat until convergence.
 *
 * This models the same physical problem as a monolithic Newton solve
 * but uses Rodin's linear form language for the mechanics block and
 * explicit scalar algebra for the internal variables.
 */
#include <cstddef>
#include <cmath>
#include <iostream>

#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
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
  // ec, gamma, beta : P0g (single scalar DOF over whole mesh)
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

    // Smooth prescribed activation in [0, 1] (spatially uniform here for simplicity)
    const Real activation = 0.5 * (1.0 + std::sin(2.0 * pi * time));

    // Previous step values for the scalar P0g variables
    const Real ecPrev    = ec_n.getData()(0);
    const Real gammaPrev = gamma_n.getData()(0);
    const Real betaPrev  = beta_n.getData()(0);

    // Current iterate values for the scalar variables
    Real ecCur    = ecPrev;
    Real gammaCur = gammaPrev;
    Real betaCur  = betaPrev;

    // Staggered fixed-point iteration
    constexpr size_t maxIter = 20;
    constexpr Real tol = 1e-10;
    for (size_t iter = 0; iter < maxIter; ++iter)
    {
      // ---- Block 1: Solve mechanics for u given current ec ----------------
      // Residual:
      //   ∫ μ ∇u:∇v + λ div(u) div(v) + Es div(u) div(v) = Es ec div(v)
      // where ec is treated as known scalar.
      {
        TrialFunction u(Vh);
        TestFunction  v(Vh);

        // Active coupling: Es * ec * div(v) is a body-like linear form
        // (ec is scalar, div(v) integrates over domain)
        auto ecFunc = RealFunction(ecCur);

        Problem mechanics(u, v);
        mechanics =
            Integral(mu * Jacobian(u), Jacobian(v))
          + Integral(lambda * Div(u), Div(v))
          + Integral(Es * Div(u), Div(v))
          - Integral(Es * ecFunc, Div(v))
          + DirichletBC(u, VectorFunction{ Zero(), Zero() }).on(leftBC);

        mechanics.assemble();
        SparseLU(mechanics).solve();

        u_n = u.getSolution();
      }

      // Compute volume-averaged div(u) as kinematic proxy for fiber strain.
      // For a P0g variable (one global DOF), the relevant coupling quantity
      // is ∫ div(u) dΩ / |Ω|.
      Real divU = 0.0;
      Real vol  = 0.0;
      {
        const auto& conn = mesh.getConnectivity();
        const size_t d = mesh.getDimension();
        const size_t nCells = conn.getCount(d);
        for (size_t c = 0; c < nCells; ++c)
        {
          Polytope cell(d, c, mesh);
          const Real cellVol = cell.getMeasure();

          // Evaluate div(u) at cell centroid using the FE basis
          const auto& fes = u_n.getFiniteElementSpace();
          const auto& fe = fes.getFiniteElement(d, c);
          const size_t ndof = fe.getCount();
          const size_t vdim = fes.getVectorDimension();

          // Get the centroid quadrature point
          const auto& qf = QF::PolytopeQuadratureFormula::get(0, cell.getGeometry());
          const auto& quadrature = cell.getQuadrature(qf);
          const auto& pt = quadrature.getPoint(0);
          const auto& rc = qf.getPoint(0);
          const auto& JacInv = pt.getJacobianInverse();

          Real cellDiv = 0.0;
          for (size_t dof = 0; dof < ndof; ++dof)
          {
            auto refJac = fe.getBasis(dof).getJacobian()(rc);
            auto physJac = refJac * JacInv;
            const Real uDof = u_n.getData()(fes.getGlobalIndex({d, c}, dof));
            for (size_t k = 0; k < std::min(vdim, d); ++k)
              cellDiv += uDof * physJac(k, k);
          }

          divU += cellDiv * cellVol;
          vol  += cellVol;
        }
      }
      divU /= vol;

      // ---- Blocks 2–4: Update scalar internal variables -------------------
      // Active strain rate (backward Euler)
      const Real edot = (ecCur - ecPrev) / dt;

      // Regularized damping
      const Real damp = activation + alphaD * edot * edot;

      // Block 2: active branch constraint
      //   gamma * beta + muParallel * edot - Es * (divU - ec) = 0
      //   => ec = divU - (gamma * beta + muParallel * edot) / Es
      const Real ecNew = divU - (gammaCur * betaCur + muParallel * edot) / Es;

      // Block 3: gamma evolution
      //   (gamma - gammaPrev) / dt = n0 * k0 * activation - damp * gamma
      //   => gamma * (1 + dt * damp) = gammaPrev + dt * n0 * k0 * activation
      const Real gammaNew = (gammaPrev + dt * n0 * k0 * activation)
                          / (1.0 + dt * damp);

      // Block 4: beta evolution
      //   (beta - betaPrev) / dt = n0 * sigma0 * activation - damp * beta + edot * gamma
      //   => beta * (1 + dt * damp) = betaPrev + dt * (n0 * sigma0 * activation + edot * gamma)
      const Real betaNew = (betaPrev + dt * (n0 * sigma0 * activation + edot * gammaNew))
                         / (1.0 + dt * damp);

      // Check convergence
      const Real diffEc    = std::abs(ecNew - ecCur);
      const Real diffGamma = std::abs(gammaNew - gammaCur);
      const Real diffBeta  = std::abs(betaNew - betaCur);
      const Real change = std::max({diffEc, diffGamma, diffBeta});

      ecCur    = ecNew;
      gammaCur = gammaNew;
      betaCur  = betaNew;

      if (change < tol)
        break;
    }

    // Commit n+1 -> n
    ec_n.getData()(0)    = ecCur;
    gamma_n.getData()(0) = gammaCur;
    beta_n.getData()(0)  = betaCur;

    std::cout << "Step " << step
              << "  t=" << time
              << "  ec=" << ecCur
              << "  gamma=" << gammaCur
              << "  beta=" << betaCur
              << std::endl;

    xdmf.write(time);
  }

  return 0;
}
