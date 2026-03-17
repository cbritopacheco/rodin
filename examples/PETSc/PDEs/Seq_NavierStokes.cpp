/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @example Transient 2D incompressible Navier–Stokes with Taylor–Hood elements and PETSc
 *
 * This example solves a transient incompressible flow problem in 2D using:
 *   - a mixed finite element discretization,
 *   - Taylor–Hood spaces (quadratic velocity, linear pressure),
 *   - a semi-implicit Oseen / Picard linearization of the convective term,
 *   - backward Euler time stepping,
 *   - PETSc for linear algebra and iterative solution.
 *
 * The weak problem solved at each time step is:
 *   Find (u^{n+1}, p^{n+1}) such that, for all test functions (v, q),
 *
 *     (rho / dt) (u^{n+1}, v)
 *   + rho ((u^n · ∇) u^{n+1}, v)
 *   + (rho / 2) ((∇ · u^n) u^{n+1}, v)
 *   + mu (∇u^{n+1}, ∇v)
 *   - (p^{n+1}, ∇ · v)
 *   + (∇ · u^{n+1}, q)
 *
 *   = (rho / dt) (u^n, v)
 *   - <p_out, v · n>_{Γ_out}
 *   + <(rho / 2) beta u^{n+1}, v>_{Γ_out},
 *
 * where:
 *   - u^n is the velocity from the previous time step,
 *   - the convective term is linearized by freezing the transport velocity at u^n,
 *   - the additional (∇·u^n) term yields the skew-symmetric Oseen form,
 *   - beta = max(-(u^n · n), 0) penalizes backflow at the outlet.
 *
 * Boundary conditions:
 *   - inlet: prescribed pulsatile parabolic velocity profile,
 *   - wall: homogeneous Dirichlet condition,
 *   - outlet: prescribed normal traction through a Windkessel-like resistance law
 *             p_out = p_d - R_d * Q_out.
 *
 * Notes:
 *   - This example intentionally uses a stable Taylor–Hood pair, so no PSPG/SUPG
 *     stabilization is required.
 *   - The nonlinear term is not fully iterated inside each time step: one Oseen
 *     solve is performed per time level using the previous velocity as the
 *     convecting field.
 *   - The outlet model here is a simple resistance closure driven by the current
 *     volumetric flux.
 *
 * Mesh:
 *   ../resources/examples/PDEs/NavierStokes.medit.mesh
 *
 * Output:
 *   - NavierStokes.mesh
 *   - NavierStokes_velocity_XXXXXX.gf
 *   - NavierStokes_pressure_XXXXXX.gf
 *   - flux.txt
 *   - pressure.txt
 *
 * Suggested command:
 * VECLIB_MAXIMUM_THREADS=1 OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 \
 * ./examples/PETSc/PDEs/PETSc_Seq_NavierStokes \
 * -ksp_type preonly -pc_type lu -ksp_rtol 1e-6 \
 * -ksp_monitor -ksp_converged_reason \
 * -pc_factor_shift_type nonzero  -pc_factor_shift_amount 1e-10
 */

#include "Rodin/IO/XDMF.h"
#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

#include <Rodin/PETSc.h>

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>

using namespace Rodin;
using namespace Rodin::Math;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

/**
 * @brief Returns a half-wave sinusoidal inlet pulse.
 *
 * Over each period, the inflow is sinusoidal during the first half-period
 * and zero during the second half-period.
 *
 * @param t Current physical time.
 * @param freq Pulse frequency.
 * @return Instantaneous pulse amplitude.
 */
static Real InletPulse(const Real t, const Real freq)
{
  const Real period = 1.0 / freq;
  const Real halfPeriod = 0.5 * period;
  const Real tm = std::fmod(t, period);
  const Real pi = Constants::pi();
  if (tm < halfPeriod)
    return 200 * std::sin(2.0 * pi * freq * t);
  return 0.0;
}

static const char* filename = "../resources/examples/PDEs/NavierStokes.medit.mesh";

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  // Boundary attribute tags as stored in the input mesh.
  constexpr Attribute outlet = 0;
  constexpr Attribute inlet  = 1;
  constexpr Attribute wall   = 3;

  // Geometric / physical parameters.
  const Real R   = 1.5;    // Inlet radius used in the parabolic profile.
  const Real rho = 1.0;    // Fluid density.
  const Real mu  = 0.035;  // Viscosity coefficient.

  // Time discretization parameters.
  const Real T   = 12.0;   // Final time.
  const Index Nt = 400;    // Number of time steps.
  const Real dt  = T / Nt; // Time step size.

  // Outlet pressure model:
  //   p_out = p_d - R_d * Q_out
  // where p_d is a distal pressure and R_d is a resistance parameter.
  const Real pd = 8.0 * 13332.2;
  const std::array<Real, 5> RdValues{{0.0, 100.0, 200.0, 300.0, 800.0}};
  const Index idxRd = 4;
  const Real Rd = RdValues[idxRd];

  // Pulse frequency for the inlet waveform.
  const Real freq = 30.0 / 60.0;

  // Load the mesh and compute boundary-to-cell connectivity.
  Mesh mesh;
  mesh.load(filename, IO::FileFormat::MEDIT);
  mesh.getConnectivity().compute(1, 2);

  // Save the mesh once for post-processing convenience.
  Alert::Info() << "Saving MEDIT mesh to NavierStokes.mesh..." << Alert::Raise;
  mesh.save("NavierStokes.mesh", IO::FileFormat::MEDIT);

  Alert::Info() << "Setting up NavierStokes.xdmf ..." << Alert::Raise;
  IO::XDMF xdmf("NavierStokes");
  xdmf.setMesh(mesh);

  const size_t dim = mesh.getSpaceDimension();

  // Taylor–Hood discretization:
  //   velocity in H1 of degree 2,
  //   pressure in H1 of degree 1.
  //
  // This is an inf-sup stable mixed pair for incompressible flow,
  // so no pressure stabilization is needed here.
  H1 uh(std::integral_constant<size_t, 2>{}, mesh, dim);
  H1 ph(std::integral_constant<size_t, 1>{}, mesh);

  {
    // Unknowns and test functions.
    PETSc::Variational::TrialFunction u(uh); u.setName("u");
    PETSc::Variational::TrialFunction p(ph); p.setName("p");
    PETSc::Variational::TestFunction  v(uh);
    PETSc::Variational::TestFunction  q(ph);

    xdmf.add("velocity", u.getSolution());
    xdmf.add("pressure", p.getSolution());

    // Solution fields from the previous time step.
    PETSc::Variational::GridFunction u_old(uh);
    PETSc::Variational::GridFunction p_old(ph);

    // Time-dependent inlet velocity profile.
    PETSc::Variational::GridFunction inletProfile(uh);

    // Initial conditions: fluid initially at rest and zero pressure.
    u_old = Math::Vector<Real>{{0.0, 0.0}};
    p_old = 0.0;
    inletProfile = Math::Vector<Real>{{0.0, 0.0}};

    // Diagnostic output files.
    std::ofstream fluxFile("flux.txt");
    std::ofstream pressureFile("pressure.txt");
    if (!fluxFile || !pressureFile)
    {
      std::cerr << "Failed to open diagnostic output files flux.txt / pressure.txt.\n";
      PetscFinalize();
      return 1;
    }

    // Current outlet pressure used in the traction boundary condition.
    Real pout = pd;

    // Current simulation time.
    Real t = 0.0;

    // Helper objects for computing the outlet volumetric flow rate:
    //   Q_out = ∫_{Γ_out} u · n
    PETSc::Variational::GridFunction one(ph);
    one = 1.0;
    PETSc::Variational::TestFunction qFlux(ph);
    LinearForm flux(qFlux);

    // Outward unit normal on boundary faces.
    auto n = BoundaryNormal(mesh);

    for (Index k = 0; k < Nt; k++)
    {
      t += dt;

      // Build the prescribed pulsatile inlet profile at the current time.
      //
      // Here the inflow is vertical and parabolic in x:
      //   u_in(x, t) = (0, g(t) * (R^2 - x^2) / R^2)
      //
      // This assumes the inlet cross section and coordinate convention of the
      // supplied mesh. If the geometry changes, this profile may need revision.
      const Real gt = InletPulse(t, freq);
      inletProfile = [&](const Point& x)
      {
        const Real uy = gt * (R * R - x.x() * x.x()) / (R * R);
        return Math::Vector<Real>{{0.0, uy}};
      };

      // Semi-implicit Oseen / Picard linearization:
      //
      // The nonlinear convective term (u · ∇)u is replaced by
      //   (u_old · ∇)u
      // at the new time step.
      //
      // In Rodin notation:
      //   Jacobian(u) * u_old = (u_old · ∇)u.
      const auto conv_u = Mult(Jacobian(u), u_old);

      // Divergence of the frozen transport velocity.
      //
      // This is used in the skew-symmetric correction
      //   0.5 * (div u_old) * u · v
      // to obtain an energy-friendlier Oseen form.
      const auto div_u_old = Div(u_old);

      // Backflow coefficient at the outlet:
      //   beta = max(-(u_old · n), 0)
      //
      // When the flow tries to re-enter through the outlet, beta becomes positive
      // and activates an additional boundary damping term.
      const auto beta = Max(-Dot(u_old, n), 0.0);

      // Mixed Oseen problem at the current time step.
      Problem flow(u, p, v, q);
      flow =
          // Backward Euler time derivative:
          //   (rho / dt) (u - u_old, v)
          (rho / dt) * Integral(u, v)
        - (rho / dt) * Integral(u_old, v)

          // Linearized convection:
          //   rho ((u_old · ∇)u, v)
        + rho * Integral(Dot(conv_u, v))

          // Skew-symmetric correction:
          //   (rho / 2) ((div u_old) u, v)
        + 0.5 * rho * Integral(div_u_old * Dot(u, v))

          // Viscous term:
          //   mu (∇u, ∇v)
        + mu * Integral(Jacobian(u), Jacobian(v))

          // Pressure-velocity coupling:
          //   -(p, div v)
        - Integral(p, Div(v))

          // Incompressibility constraint:
          //   (div u, q)
        + Integral(Div(u), q)

          // Outlet pressure traction:
          //   -<p_out, v · n>_{Γ_out}
        - BoundaryIntegral(pout * Dot(v, n)).over(outlet)

          // Backflow damping on the outlet:
          //   <(rho / 2) beta u, v>_{Γ_out}
        + BoundaryIntegral(0.5 * rho * beta * Dot(u, v)).over(outlet)

        // Tiny pressure-block diagonal filler used to
        // improve robustness of sparse direct factorization
        // for the mixed system.
        + 1e-6 * Integral(p, q)

          // No-slip wall condition.
        + DirichletBC(u, Zero(dim)).on(wall)

          // Prescribed inlet velocity.
        + DirichletBC(u, inletProfile).on(inlet);

      Alert::Info() << "Assembling linear system for time step " << k + 1 << " / " << Nt << "..."
        << Alert::Raise;

      // Assemble the linear system and define PETSc field splits for the mixed
      // velocity-pressure block structure.
      flow.assemble().setFieldSplits();

      Alert::Info() << "Solving linear system for time step " << k + 1 << " / " << Nt << "..."
        << Alert::Raise;

      // Solve the assembled linear system.
      Solver::KSP(flow).solve();

      // Advance the time history by copying the newly computed solution into
      // the "old" fields used by the next step.
      u_old = u.getSolution();
      p_old = p.getSolution();

      // Compute the outlet flow rate
      //   Q_out = ∫_{Γ_out} u · n
      flux = BoundaryIntegral(Dot(u_old, n), qFlux).over(outlet);
      flux.assemble();
      const Real qout = flux(one);

      // Update the outlet pressure through the resistance law:
      //   p_out = p_d - R_d * Q_out
      //
      // This is a simple lumped downstream model.
      pout = pd - Rd * qout;

      // Save scalar diagnostics for later plotting.
      fluxFile << t << " " << qout << "\n";
      pressureFile << t << " " << pout << "\n";

      // Save solution snapshots at every time step.
      xdmf.write(t).flush();
    }
  }

  Alert::Success() << "Simulation completed. Closing XDMF file..." << Alert::Raise;

  xdmf.close();

  PetscFinalize();

  return 0;
}
