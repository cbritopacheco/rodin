/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @example Transient 2D Navier–Stokes / Oseen-like problem with PETSc mixed assembly
 *
 * Mesh:
 *   ../resources/examples/PDEs/NavierStokes.medit.mesh
 *
 * Suggested PETSc options:
 *   -ksp_type gmres -pc_type lu
 *   -ksp_type fgmres -pc_type fieldsplit
 */

#include <Rodin/PETSc.h>

#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

#include <petscksp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace Rodin;
using namespace Rodin::Math;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

static Real InletPulse(const Real t, const Real freq)
{
  const Real period = 1.0 / freq;
  const Real halfPeriod = 0.5 * period;
  const Real tm = std::fmod(t, period);
  const Real pi = Constants::pi();
  if (tm < halfPeriod)
    return 200.0 * std::sin(2.0 * pi * freq * t);
  return 0.0;
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  constexpr Attribute outlet = 0;
  constexpr Attribute inlet = 1;
  constexpr Attribute wall = 3;

  const Real R = 1.5;
  const Real r = 5.0;
  const Real rho = 1.0;
  const Real mu = 0.035;
  const Real T = 12.0;
  const Index Nt = 200;
  const Real dt = T / Nt;
  const Real pd = 8.0 * 13332.2;
  const std::array<Real, 5> RdValues{{0.0, 100.0, 200.0, 300.0, 800.0}};
  const Index idxRd = 4;
  const Real Rd = RdValues[idxRd];
  const Real freq = 50.0 / 60.0;

  Mesh mesh;
  mesh.load("../resources/examples/PDEs/NavierStokes.medit.mesh");
  mesh.getConnectivity().compute(1, 2);
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);
  mesh.getConnectivity().compute(0, 1);
  mesh.save("NavierStokes.mesh");

  const size_t dim = mesh.getSpaceDimension();
  P1 uh(mesh, dim);
  P1 ph(mesh);
  P0 kh(mesh);

  PETSc::Variational::TrialFunction u(uh); u.setName("u");
  PETSc::Variational::TrialFunction p(ph); p.setName("p");
  PETSc::Variational::TestFunction  v(uh);
  PETSc::Variational::TestFunction  q(ph);

  PETSc::Variational::GridFunction u_old(uh);
  PETSc::Variational::GridFunction p_old(ph);
  PETSc::Variational::GridFunction tauK(kh);
  PETSc::Variational::GridFunction beta(kh);
  PETSc::Variational::GridFunction inletProfile(uh);

  u_old = Math::Vector<Real>{{0.0, 0.0}};
  p_old = 0.0;
  tauK = 0.0;
  beta = 0.0;
  inletProfile = Math::Vector<Real>{{0.0, 0.0}};

  std::ofstream fluxFile("flux.txt");
  std::ofstream pressureFile("pressure.txt");
  if (!fluxFile || !pressureFile)
  {
    std::cerr << "Failed to open diagnostic output files flux.txt / pressure.txt.\n";
    PetscFinalize();
    return 1;
  }

  Real pout = pd;
  Real t = 0.0;

  PETSc::Variational::GridFunction one(ph);
  one = 1.0;
  PETSc::Variational::TestFunction qFlux(ph);
  LinearForm fluxForm(qFlux);

  for (Index n = 0; n < Nt; n++)
  {
    t += dt;

    const auto uoldNorm2 = Dot(u_old, u_old);
    tauK = [&](const Point& x)
    {
      const Real umag2 = uoldNorm2(x);
      // TODO: replace hK = sqrt(|K|) by a directional cell diameter
      // (e.g. max edge length) once a mesh-quality diameter utility is exposed.
      const Real hK = std::sqrt(x.getPolytope().getMeasure());
      const Real hK2 = hK * hK;
      const Real denom = std::sqrt(4.0 * umag2 / hK2 + 16.0 * mu * mu / (hK2 * hK2));
      return (denom > 0.0) ? 0.1 / denom : 0.0;
    };

    const auto uyOld = u_old.y();
    beta = [&](const Point& x)
    {
      return std::max(-uyOld(x), Real(0));
    };

    const Real gt = InletPulse(t, freq);
    inletProfile = [&](const Point& x)
    {
      const Real uy = gt * (R * R - x.x() * x.x()) / (R * R);
      return Math::Vector<Real>{{0.0, uy}};
    };

    const auto conv_u = Mult(Jacobian(u), u_old);
    const auto conv_v = Mult(Jacobian(v), u_old);
    const auto grad_p = Grad(p);
    const auto grad_q = Grad(q);

    Problem flow(u, p, v, q);
    flow =
        (1.0 / dt) * Integral(u, v)
      - (1.0 / dt) * Integral(u_old, v)
      + Integral(Dot(conv_u, v))
      + mu * Integral(Jacobian(u), Jacobian(v))
      - Integral(p, Div(v))
      + Integral(Div(u), q)
      + 0.5 * Integral(Div(u_old) * Dot(u, v))
      - BoundaryIntegral(pout * v.y()).over(outlet)
      + BoundaryIntegral(0.5 * rho * beta * Dot(u, v)).over(outlet)

        // SUPG / PSPG expanded
      + Integral(tauK * Dot(conv_u, conv_v))
      + Integral(tauK * Dot(conv_u, grad_q))
      + Integral(tauK * Dot(grad_p, conv_v))
      + Integral(tauK * Dot(grad_p, grad_q))

      + DirichletBC(u, Zero(dim)).on(wall)
      + DirichletBC(u, inletProfile).on(inlet);

    flow.assemble();
    flow.setFieldSplits();

    Solver::KSP(flow).solve();

    // Use setData to copy PETSc vector entries without re-creating GridFunction storage.
    u_old.setData(u.getSolution().getData());
    p_old.setData(p.getSolution().getData());

    fluxForm = BoundaryIntegral(u_old.y(), qFlux).over(outlet);
    fluxForm.assemble();
    const Real qout = fluxForm(one);
    pout = pd - Rd * qout;

    fluxFile << t << " " << qout << "\n";
    pressureFile << t << " " << pout << "\n";

    std::ostringstream id;
    id << std::setfill('0') << std::setw(6) << n;
    u_old.save("NavierStokes_velocity_" + id.str() + ".gf");
    p_old.save("NavierStokes_pressure_" + id.str() + ".gf");
  }

  PetscFinalize();
  return 0;
}
