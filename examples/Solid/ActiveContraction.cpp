/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ActiveContraction.cpp
 * @brief Hill–Maxwell-inspired active contraction with explicit local internal updates.
 *
 * This example keeps active variables as FE fields (P0) and avoids Schur/local
 * elimination in the global mechanics solve. Internal variable updates are done
 * per-cell through a local nonlinear fixed-point/Newton loop.
 */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

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
  constexpr size_t nx = 65;
  constexpr size_t ny = 17;

  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { nx, ny });
  mesh.scale(1.0 / static_cast<Real>(ny - 1));
  mesh.getConnectivity().compute(1, 2);

  constexpr Attribute leftBC = 1;
  constexpr Real eps = 1e-10;
  for (auto it = mesh.getBoundary(); !it.end(); ++it)
  {
    const auto& verts = it->getVertices();
    Real xSum = 0;
    for (size_t i = 0; i < verts.size(); ++i)
      xSum += mesh.getVertexCoordinates(verts[i])(0);
    if (xSum / static_cast<Real>(verts.size()) < eps)
      mesh.setAttribute({ 1, it->getIndex() }, leftBC);
  }

  const size_t dim = mesh.getSpaceDimension();
  const size_t d = mesh.getDimension();

  // Mechanics on vector P1, internal variables on cellwise P0.
  P1 Vh(mesh, dim);
  P0 Qh(mesh);

  GridFunction u_n(Vh);
  u_n.setName("Displacement");
  u_n = VectorFunction{ Zero(), Zero() };

  GridFunction ec_n(Qh);    ec_n.setName("ec");    ec_n = RealFunction(0.0);
  GridFunction gamma_n(Qh); gamma_n.setName("gamma"); gamma_n = RealFunction(1.0);
  GridFunction beta_n(Qh);  beta_n.setName("beta");  beta_n = RealFunction(0.0);
  GridFunction activation_n(Qh); activation_n.setName("activation"); activation_n = RealFunction(0.0);

  IO::XDMF xdmf("ActiveContraction");
  auto grid = xdmf.grid();
  grid.setMesh(mesh);
  grid.add(u_n);
  grid.add(ec_n);
  grid.add(gamma_n);
  grid.add(beta_n);
  grid.add(activation_n);

  // Material/model parameters.
  const Real E  = 120.0;
  const Real nu = 0.35;
  const Real lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  const Real mu     = E / (2.0 * (1.0 + nu));

  const Real Es         = 40.0;
  const Real muParallel = 0.5;
  const Real alphaD     = 0.2;
  const Real k0         = 1.0;
  const Real sigma0     = 1.0;

  constexpr size_t nSteps = 20;
  const Real T = 1.0;
  const Real dt = T / static_cast<Real>(nSteps);
  const Real pi = Math::Constants::pi();

  GridFunction ecPrev(Qh), gammaPrev(Qh), betaPrev(Qh);
  GridFunction ecCur(Qh), gammaCur(Qh), betaCur(Qh);

  // Frank-Starling reduction factor (piecewise linear; clamped at 0).
  const auto starling = [](Real fib0) -> Real
  {
    const Real x1 = -0.4, y1 = 0.0;
    const Real x2 =  0.3, y2 = 0.38;
    const Real x3 =  0.73, y3 = 0.74;
    const Real x4 =  1.0, y4 = 1.0;
    const Real x5 =  1.3, y5 = 1.0;
    const Real x6 =  2.4, y6 = 0.0;

    Real n0 = 0.0;
    if (fib0 < x2) n0 = ((y2-y1)/(x2-x1)) * (fib0 - x2) + y2;
    else if (fib0 < x3) n0 = ((y3-y2)/(x3-x2)) * (fib0 - x3) + y3;
    else if (fib0 < x4) n0 = ((y4-y3)/(x4-x3)) * (fib0 - x4) + y4;
    else if (fib0 < x5) n0 = y4;
    else if (fib0 < x6) n0 = ((y6-y5)/(x6-x5)) * (fib0 - x6) + y6;
    return std::max<Real>(0.0, n0);
  };

  const size_t nCells = mesh.getConnectivity().getCount(d);
  std::vector<Real> e1D(nCells, 0.0);
  const Real invSqrt2 = 1.0 / std::sqrt(2.0);

  for (size_t step = 1; step <= nSteps; ++step)
  {
    const Real time = dt * static_cast<Real>(step);
    const Real activation = 0.5 * (1.0 + std::sin(2.0 * pi * time));
    const Real activationPlus = std::max<Real>(0.0, activation);
    for (size_t c = 0; c < nCells; ++c)
      activation_n.getData()(c) = activation;

    ecPrev = ec_n; gammaPrev = gamma_n; betaPrev = beta_n;
    ecCur = ecPrev; gammaCur = gammaPrev; betaCur = betaPrev;

    constexpr size_t maxIter = 25;
    constexpr Real tol = 1e-10;

    for (size_t iter = 0; iter < maxIter; ++iter)
    {
      // 1) Mechanics solve for u given current ec.
      {
        TrialFunction u(Vh);
        TestFunction  v(Vh);

        Problem mechanics(u, v);
        mechanics =
            Integral(mu * Jacobian(u), Jacobian(v))
          + Integral(lambda * Div(u), Div(v))
          + Integral(Es * Div(u), Div(v))
          - Integral(Es * ecCur, Div(v))
          + DirichletBC(u, VectorFunction{ Zero(), Zero() }).on(leftBC);

        mechanics.assemble();
        SparseLU(mechanics).solve();
        u_n = u.getSolution();
      }

      // 2) Cell-wise e1D proxy from centroid div(u) + moved-fiber stretch.
      //    We transport a reference fiber a0=(1,1)/sqrt(2) through a simple
      //    deformation proxy F = diag(1+div(u), 1) and set e1D = |F a0| - 1.
      {
        const auto& fes = u_n.getFiniteElementSpace();
        for (size_t c = 0; c < nCells; ++c)
        {
          Polytope cell(d, c, mesh);
          const auto& fe = fes.getFiniteElement(d, c);
          const size_t ndof = fe.getCount();
          const size_t vdim = fes.getVectorDimension();

          const auto& qf = QF::PolytopeQuadratureFormula::get(0, cell.getGeometry());
          const auto& quadrature = cell.getQuadrature(qf);
          const auto& pt = quadrature.getPoint(0);
          const auto& rc = qf.getPoint(0);
          const auto& JacInv = pt.getJacobianInverse();

          Real divCell = 0.0;
          for (size_t dof = 0; dof < ndof; ++dof)
          {
            const auto refJac  = fe.getBasis(dof).getJacobian()(rc);
            const auto physJac = refJac * JacInv;
            const Real uDof = u_n.getData()(fes.getGlobalIndex({d, c}, dof));
            for (size_t k = 0; k < std::min(vdim, d); ++k)
              divCell += uDof * physJac(k, k);
          }
          const Real fax = (1.0 + divCell) * invSqrt2;
          const Real fay = invSqrt2;
          const Real lambdaF = std::sqrt(fax * fax + fay * fay);
          e1D[c] = lambdaF - 1.0;
        }
      }

      // 3) Local update of (ec, gamma, beta) per cell.
      auto& ecData = ecCur.getData();
      auto& gammaData = gammaCur.getData();
      auto& betaData = betaCur.getData();
      const auto& ec0Data = ecPrev.getData();
      const auto& gamma0Data = gammaPrev.getData();
      const auto& beta0Data = betaPrev.getData();

      Real maxDiff = 0.0;
      for (size_t c = 0; c < nCells; ++c)
      {
        Real ec = ecData(c);
        Real gamma = gammaData(c);
        Real beta = betaData(c);

        const Real ec0 = ec0Data(c);
        const Real gamma0 = gamma0Data(c);
        const Real beta0 = beta0Data(c);
        const Real n0 = starling(ec0);

        for (size_t lit = 0; lit < 8; ++lit)
        {
          const Real edot = (ec - ec0) / dt;

          // Positivity-style gamma/beta updates.
          const Real Dg = 1.0 + dt * (std::abs(activation) + alphaD * std::abs(edot));
          const Real gammaSq = std::max<Real>(1e-16, (gamma0 * gamma0 + dt * n0 * k0 * activationPlus) / Dg);
          gamma = std::sqrt(gammaSq);

          const Real Db = 1.0
            + 0.5 * dt * n0 * k0 * activationPlus / gammaSq
            + 0.5 * dt * (std::abs(activation) + alphaD * std::abs(edot));

          beta = (beta0 + gamma * (ec - ec0) + dt * n0 * sigma0 * activationPlus / gamma) / Db;

          // Solve RA(ec)=0 with tau = gamma*beta (frozen during local Newton).
          const Real tau = gamma * beta;
          for (size_t nit = 0; nit < 6; ++nit)
          {
            const Real onep2ec = 1.0 + 2.0 * ec;
            const Real onep2e1D = 1.0 + 2.0 * e1D[c];

            const Real R =
              (tau + muParallel * (ec - ec0) / dt) * std::pow(onep2ec, 3)
              - Es * (e1D[c] - ec) * onep2e1D;

            const Real dR =
              (muParallel / dt) * std::pow(onep2ec, 3)
              + (tau + muParallel * (ec - ec0) / dt) * 6.0 * std::pow(onep2ec, 2)
              + Es * onep2e1D;

            const Real dec = -R / dR;
            ec += dec;
            ec = std::max<Real>(ec, -0.49);
            if (std::abs(dec) < 1e-12)
              break;
          }
        }

        maxDiff = std::max(maxDiff, std::abs(ec - ecData(c)));
        maxDiff = std::max(maxDiff, std::abs(gamma - gammaData(c)));
        maxDiff = std::max(maxDiff, std::abs(beta - betaData(c)));

        ecData(c) = ec;
        gammaData(c) = gamma;
        betaData(c) = beta;
      }

      if (maxDiff < tol)
        break;
    }

    ec_n = ecCur;
    gamma_n = gammaCur;
    beta_n = betaCur;

    std::cout << "Step " << step << " t=" << time << std::endl;
    xdmf.write(time).flush();
  }

  xdmf.close();
  return 0;
}
