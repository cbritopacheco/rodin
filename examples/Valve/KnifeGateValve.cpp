/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
// KnifeGateValve.cpp
//
// Hydrodynamic force on the gate (knife) of a knife gate valve at a fixed
// opening. Incompressible Navier-Stokes with a power-law viscosity in the
// fluid of a straight pipe with the knife removed (meshes from
// examples/Valve/make_meshes.py):
//
//   rho (u - u^n)/dt + rho (grad u) u - div(2 mu(gd) e(u)) + grad p = 0,
//   div u = 0,
//   mu(gd) = K max(gd, gd_min)^(n - 1), clipped to [mu_min, mu_max],
//   gd = sqrt(2 e(u):e(u)),
//
//   u = U w(r)/<w> e      on the inlet  (velocity inlet, default: developed
//                         profile, laminar power-law or 1/7 turbulent law,
//                         flux exactly U A_in), or
//   sigma n = -p_in n     on the inlet  (pressure inlet, u x n = 0 by penalty;
//                         the only option for the closed gate),
//   sigma n = -p_out n    on the outlet (free outflow, directional do-nothing:
//                         + rho/2 max(-u.n, 0) u against backflow; default
//                         p_out = 150 psi, the CWP line pressure),
//   u = 0                 on the pipe wall and on the knife.
//
// Discretisation: variable-step BDF2 (backward Euler on the first step,
// -kv_bdf 1 throughout); P1/P1 for (u, p) with split orthogonal subgrid
// scales (OSS, Codina 2000), the projections Pi_h onto P1 with the lumped
// mass matrix (-kv_lumped false: consistent L2, OrthogonalProjection.h):
//
//   + tau_M rho ((grad u) a - Pi[(grad u^n) a^n], (grad v) a)
//   + tau_M/rho (grad p - Pi[grad p^n], grad q)
//   + tau_C (div u - Pi[div u^n], div v),
//
//   tau_M = [(2/dt)^2 + (2 |u^k|/h)^2 + (4 nu^k/h^2)^2]^(-1/2),
//   tau_C = mu^k + rho |u^k| h / 2,   h = (6 sqrt(2) |T|)^(1/3),
//
// a the cell mean of u^k (P0; the Galerkin convection keeps (grad u) u^k),
// and the projections lagged at t^n (refreshed once per step: lagged at the
// iterate they make Picard contract only like tau/(tau + dt), see
// examples/Heat/ForcedConvection.cpp). The nonlinearity is resolved by Picard
// iterations within each step, starting from the linear extrapolation of
// (u^{n-1}, u^n): convection velocity u^k, viscosity mu(u^k), tau(u^k), all
// coefficients projected onto P0 (constant per cell, exact for P1
// velocities) so every integrand has a known polynomial order.
//
// Cost. Per Picard iteration: one assembly and one linear solve, GMRES
// preconditioned by MUMPS factors that are reused while GMRES converges in
// a few iterations (-kv_refactor_every, -kv_refactor_its): the factorisation
// was ~40 % of the step on an 8-core Mac (FPC8 25 %, 430 k cells), the
// lagged factors cut the solve 4-5x with the same forces to 1e-6. The
// assembly then dominates (~130-185 us per cell and assembly): each
// integrator walks the mesh on its own.
// The log prints the time per phase at every step and the load balance of
// the partition at start-up; without Scotch the fallback partitioner can
// leave the whole mesh on rank 0.
//
// Force on the gate. F = -int_knife sigma n dS (n the outward normal of the
// fluid), the force exerted by the fluid on the knife, global axes: x
// lateral, y along the stroke (+y opens), z along the flow. Two evaluations:
//
//  (1) Consistent reaction (primary). The discrete momentum residual with no
//      Dirichlet elimination, R_k(z) = a((u_h, p_h), z e_k) (every volume term
//      of the discrete equations, stabilisation included), equals
//      int_dOmega (sigma n)_k z for every P1 z, so the force is the sum of
//      the nodal reactions over the knife,
//          F_k = -sum_{i on the knife} (R_k(phi_i) + int_wall p n_k phi_i),
//      the second term removing, at the knife-pipe junction nodes, the
//      pressure of the pipe wall faces that share the node (exact for P1 p).
//      Only volume integrals of u_h, p_h appear: no velocity gradient is
//      evaluated on a boundary face. Closed FPC8 (coarse mesh): F_z equals
//      (p_in - p_out) times the wetted area of the knife to 1e-6.
//
// A uniform pressure p_out on the wetted knife adds p_out int_knife n dS to
// F: zero in z (both faces wetted) but not in y, where the edge of the knife
// is wetted and its sides (inside the slot of the sleeve) are not. That part
// depends on the pressure in the valve chest, which is not modelled; F_dyn =
// F - p_out int_knife n dS removes it and is the hydrodynamic force.
//
//  (2) Direct surface integral (check). int_knife (p n - mu (grad u +
//      grad u^T) n) dS with grad u taken from the tetrahedron adjacent to each
//      face. Rodin's Jacobian(u) evaluated at a boundary IntegrationPoint uses
//      the FACE element (tangential gradient only), which vanishes under the
//      no-slip trace: written as Integral(mu (J + J^T) n) the viscous force is
//      identically zero. Here the gradient is evaluated through a Point (in a
//      RealFunction), the path that resolves the face to its incident cell.
//
// (1) and (2) differ by O(h) in the viscous part (the cell gradient of P1 is
// first order at the wall); the pressure parts coincide to quadrature error.
//
// Run (from the build directory; the mesh is in metres):
//   water at 3 m/s in a 150 psi line, open gate:
//   mpirun -n 8 ./examples/Valve/KnifeGateValve -kv_mesh FPC8_open050.mesh -kv_umean 3
//   closed gate holding 150 psi:
//   mpirun -n 8 ./examples/Valve/KnifeGateValve -kv_mesh FPC8_open000.mesh \
//       -kv_inlet pressure -kv_pin 1.0342e6 -kv_pout 0 -kv_tend 0.1
//   slurry (power law):
//   ... -kv_rho 1400 -kv_K 0.3 -kv_n 0.6 -kv_umean 3
//
// Options (prefix -kv_): mesh, out, rho, K, n, gdmin, mumin, mumax,
// inlet (velocity|pressure), umean, profile (auto|laminar|turbulent|uniform), pin, pout,
// ramp, tangential (inlet penalty, kg/(m^2 s)), total_pressure (bool), dt,
// dtmin, dtmax, tend, adaptive_dt (bool), grow, shrink, easy, hard, maxit,
// tol, tau_dt (bool), bdf (1|2), extrapolate (bool), lumped (bool), refactor_every, refactor_its, picard_monitor (bool), tavg, steady_tol, steady_window, output_every, scale.
// Solvers: -kv_flow_ (MUMPS LU by default), -kv_proj_* (CG + Jacobi).
#include <algorithm>
#include <array>
#include <chrono>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

#include <mpi.h>
#include <petscksp.h>

#include <Rodin/Alert.h>
#include <Rodin/Configure.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/MPI.h>
#include <Rodin/PETSc.h>
#include <Rodin/Solver.h>
#include <Rodin/Types.h>
#include <Rodin/Variational.h>

#ifdef RODIN_USE_SCOTCH
#include <Rodin/Scotch/MeshPartitioner.h>
#endif

#include "OrthogonalProjection.h"

namespace Rodin::Examples::Valve
{
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  class KnifeGateValve
  {
    public:
      using MeshType = Mesh<Context::MPI>;
      using VectorFES = H1<1, Math::SpatialVector<Real>, MeshType>;
      using ScalarFES = H1<1, Real, MeshType>;
      using CellFES = P0<Real, MeshType>;
      using VectorGF = PETSc::Variational::GridFunction<VectorFES>;
      using ScalarGF = PETSc::Variational::GridFunction<ScalarFES>;
      using CellGF = PETSc::Variational::GridFunction<CellFES>;
      using VectorTrial = PETSc::Variational::TrialFunction<VectorGF, VectorFES>;
      using ScalarTrial = PETSc::Variational::TrialFunction<ScalarGF, ScalarFES>;
      using VectorTest = PETSc::Variational::TestFunction<VectorFES>;
      using ScalarTest = PETSc::Variational::TestFunction<ScalarFES>;
      using System = PETSc::Math::LinearSystem;
      using FlowProblem = Problem<System, VectorTrial, ScalarTrial, VectorTest, ScalarTest>;
      using ScalarLF = LinearForm<ScalarFES, ::Vec>;
      using VectorLF = LinearForm<VectorFES, ::Vec>;
      using Vec3 = std::array<Real, 3>;

      struct Labels
      {
        Attribute inlet = 1, outlet = 2, wall = 3, knife = 4;
      };

      struct Config
      {
        std::string meshPath = "FPC8_open050.mesh";
        std::string output = "KnifeGateValve";
        Labels labels;
        Real meshScale = 1.0;      ///< mesh coordinates -> m

        // Fluid: power law mu = K gd^(n-1). Default: water (n = 1, K = mu).
        Real rho = 998.2;          ///< kg/m^3 (water, 20 C)
        Real K = 1.002e-3;         ///< Pa s^n (water, 20 C: n = 1, K = mu)
        Real n = 1.0;              ///< flow index (n < 1 shear-thinning)
        Real gammaMin = 1.0e-2;    ///< 1/s, regularisation of gd -> 0
        Real muMin = 0.0;          ///< Pa s, clip of the apparent viscosity
        Real muMax = 1.0e3;        ///< Pa s

        // Boundary conditions.
        /// Inlet: "velocity" (developed profile of mean speed uMean; open
        /// gate) or "pressure" (sigma n = -p_in n; the only choice for the
        /// closed gate, where no flow rate can be imposed).
        std::string inlet = "velocity";
        Real uMean = 3.0;          ///< m/s, mean inlet velocity Q/A
        /// "auto" (laminar power-law profile if Re_MR < 2100, else 1/7 law),
        /// "laminar", "turbulent" (1/7 power law) or "uniform".
        std::string profile = "auto";
        Real pIn = 1.0342e6;       ///< Pa, inlet pressure (pressure inlet): 150 psi
        Real pOut = -1.0;          ///< Pa, outlet pressure (< 0: 150 psi with a velocity inlet, 0 with a pressure inlet)
        Real ramp = 0.05;          ///< s, smooth start of the inflow (or of p_in - p_out)
        Real tangentialPenalty = -1.0;  ///< kg/(m^2 s), u x n = 0 at the inlet (< 0: auto)
        bool totalPressure = false;     ///< p_in as total pressure (adds rho/2 |u.n| u at the inlet)

        // Time stepping and Picard.
        Real dt = 5.0e-4;          ///< s (Picard stops contracting at CFL ~ 5 in the gap jet)
        Real dtMin = 1.0e-6;
        Real dtMax = 5.0e-3;
        Real tEnd = 0.6;           ///< s (water at 3 m/s: the gate load settles in ~0.2-0.3 s)
        bool adaptiveDt = true;
        Real dtGrow = 1.25;
        Real dtShrink = 0.5;
        int picardEasy = 3;        ///< iterations at or below which dt grows
        int picardHard = 7;        ///< iterations above which dt shrinks
        int maxIterations = 12;
        Real tolerance = 1.0e-3;   ///< relative increment of u
        /// Include 2/dt in tau_M. Without it tau_M ~ min(h/2|u|, h^2/4 nu) is
        /// orders of magnitude above dt in water at low speed (h^2/4 nu ~ 1.5 s
        /// for h = 2.5 mm), and the lagged projections Pi[grad p^n],
        /// Pi[(grad u^n) u^n] penalise every change within the step with a
        /// weight tau/dt >> 1: the pressure lags the flow and Picard stalls
        /// (see examples/Heat/ForcedConvection.cpp). With it, tau_M <= dt/2.
        bool tauDt = true;
        int maxRetries = 4;
        bool picardMonitor = false;
        int bdf = 2;                  ///< 1: backward Euler; 2: BDF2 (variable step), BDF1 on the first step
        bool extrapolate = true;      ///< Picard starts from u^n + w (u^n - u^{n-1}), w = dt/dt_prev
        bool lumpedProjection = true; ///< OSS projections with the lumped mass matrix
        /// Lagged LU: refactor every refactorEvery solves (1: every solve,
        /// direct solver only) or when GMRES with the old factors needs more
        /// than refactorIterations iterations.
        int refactorEvery = 20;
        int refactorIterations = 15;

        // Output and stopping.
        Real tAverage = -1.0;      ///< s, start of the time average (< 0: tEnd/2)
        /// Early stop when the relative change of (F, Q) per step stays below
        /// this for steadyWindow steps. Off by default (0): the test is per
        /// step, so with small steps a slow drift passes it (FPC8 100 %:
        /// stopped at 2.7 D/U while F_y,dyn still went from 9 to 6 N).
        Real steadyTolerance = 0.0;
        int steadyWindow = 10;     ///< consecutive steps below steadyTolerance
        int outputEvery = 10;
      };

      KnifeGateValve(const Context::MPI& context, const Config& cfg)
        : m_cfg(cfg),
          m_mesh(makeMesh(context, m_cfg)),
          m_xdmf(context.getCommunicator(), m_cfg.output),
          m_vh(std::integral_constant<size_t, 1>{}, m_mesh, m_mesh.getSpaceDimension()),
          m_sh(std::integral_constant<size_t, 1>{}, m_mesh),
          m_ch(m_mesh),
          m_u(m_vh), m_p(m_sh), m_v(m_vh), m_q(m_sh),
          m_uOld(m_vh), m_uOld2(m_vh), m_uIt(m_vh), m_pIt(m_sh),
          m_muGF(m_ch), m_tauMGF(m_ch), m_tauMRhoGF(m_ch), m_tauCGF(m_ch),
          m_ax(m_ch), m_ay(m_ch), m_az(m_ch),
          m_load(m_vh), m_comp{ VectorGF(m_vh), VectorGF(m_vh), VectorGF(m_vh) },
          m_piConv(m_vh, "kv_proj_conv_", cfg.lumpedProjection),
          m_piGradP(m_vh, "kv_proj_gradp_", cfg.lumpedProjection),
          m_piDiv(m_sh, "kv_proj_div_", cfg.lumpedProjection),
          m_muOut(m_sh, "kv_proj_mu_", cfg.lumpedProjection),
          m_flow(m_u, m_p, m_v, m_q), m_flowKSP(m_flow),
          m_z(m_sh), m_one(m_sh),
          m_w(m_vh), m_flux(m_z), m_lfArea(m_w), m_lfKnife(m_w), m_lfRes(m_w), m_lfAux(m_w)
      {
        if (m_mesh.getSpaceDimension() != 3)
          throw std::runtime_error("KnifeGateValve expects a 3D tetrahedral mesh.");

        Math::SpatialVector<Real> zero(3);
        zero.setZero();
        m_uOld = zero;
        m_uOld2 = zero;
        m_uIt = zero;
        m_pIt = 0.0;
        m_one = Real(1);
        m_flowKSP.setPrefix("kv_flow_");

        const auto& L = m_cfg.labels;
        m_area[0] = boundaryIntegral(m_one, L.inlet);
        m_area[1] = boundaryIntegral(m_one, L.outlet);
        m_area[2] = boundaryIntegral(m_one, L.wall);
        m_area[3] = boundaryIntegral(m_one, L.knife);
        m_volume = volumeIntegral(m_one);
        if (m_area[0] <= 0.0 || m_area[1] <= 0.0 || m_area[3] <= 0.0)
          throw std::runtime_error("missing inlet (1), outlet (2) or knife (4) boundary in the mesh");

        // Lumped boundary areas: A_i (all boundaries; the nodal wall load
        // -R_i/A_i written to XDMF) and A_i^knife (> 0 marks the knife nodes).
        {
          const VectorFunction ones{ 1.0, 1.0, 1.0 };
          m_lfArea = BoundaryIntegral(Dot(ones, m_w)).over(
              FlatSet<Attribute>{ L.inlet, L.outlet, L.wall, L.knife });
          m_lfArea.assemble();
          m_lfKnife = BoundaryIntegral(Dot(ones, m_w)).over(L.knife);
          m_lfKnife.assemble();
          for (size_t k = 0; k < 3; ++k)
          {
            Math::SpatialVector<Real> ek(3);
            ek.setZero();
            ek(k) = 1.0;
            m_comp[k] = ek;  // 1 on the dofs of component k
          }
        }

        m_velocityInlet = (m_cfg.inlet == "velocity");
        if (!m_velocityInlet && m_cfg.inlet != "pressure")
          throw std::runtime_error("-kv_inlet must be velocity or pressure");
        if (m_cfg.pOut < 0.0)
          m_cfg.pOut = m_velocityInlet ? 1.0342e6 : 0.0;
        // Reference velocity: U (velocity inlet) or the Bernoulli velocity of p_in - p_out.
        m_uRef = m_velocityInlet ? std::abs(m_cfg.uMean)
                                 : std::sqrt(2.0 * std::abs(m_cfg.pIn - m_cfg.pOut) / m_cfg.rho);
        // int_knife n dS: a uniform pressure p_out on the wetted knife gives
        // p_out times this, which F_dyn removes (see the header).
        {
          const auto nb = BoundaryNormal(m_mesh);
          for (size_t k = 0; k < 3; ++k)
          {
            const VectorFunction e{ Real(k == 0), Real(k == 1), Real(k == 2) };
            m_knifeNormal[k] = boundaryIntegral(Dot(nb, e), L.knife);
          }
        }
        if (m_velocityInlet)
          setupInflow();
        if (!m_velocityInlet && m_cfg.tangentialPenalty < 0.0)
        {
          // Dominates rho |u| and mu/h at the inlet by ~10^3: a penalty, not a
          // stabilisation (MUMPS does not mind the conditioning).
          const Real R = std::sqrt(m_area[0] / M_PI);
          m_cfg.tangentialPenalty = 1.0e3 * std::max(m_cfg.rho * m_uRef, m_cfg.K / (0.05 * R));
        }

        // l2 norm of a dof vector with |u| = U_B at every node.
        m_velocityScale = m_uRef * std::sqrt(static_cast<Real>(m_vh.getSize()));

        updateCoefficients();
        updateStabilisation();

        m_xdmf.setMesh(m_mesh);
        m_xdmf.add("velocity", m_u.getSolution());
        m_xdmf.add("pressure", m_p.getSolution());
        m_xdmf.add("viscosity", m_muOut.get());
        m_xdmf.add("wall_load", m_load);  // force per area of the fluid on the walls (Pa)

        const auto cells = m_mesh.getCellCount(), vertices = m_mesh.getVertexCount();  // collective
        // Load balance: owned cells per rank. Without Scotch the fallback
        // partitioner can leave (almost) the whole mesh on rank 0, and then
        // every process but one idles: check this line first when a run is slow.
        PetscInt ownedCells = 0;
        {
          PetscErrorCode ierr = VecGetLocalSize(m_muGF.getData(), &ownedCells);
          assert(ierr == PETSC_SUCCESS);
          (void)ierr;
        }
        PetscInt minCells = 0, maxCells = 0;
        MPI_Allreduce(&ownedCells, &minCells, 1, MPIU_INT, MPI_MIN, PETSC_COMM_WORLD);
        MPI_Allreduce(&ownedCells, &maxCells, 1, MPIU_INT, MPI_MAX, PETSC_COMM_WORLD);
        const int ranks = m_mesh.getContext().getCommunicator().size();
        if (isRoot())
        {
          const Real D = 2.0 * std::sqrt(m_area[0] / M_PI);
          Alert::Info() << "[mesh] " << m_cfg.meshPath << "  cells=" << cells << " vertices=" << vertices
                        << "  dofs=" << 4 * m_sh.getSize()
                        << Alert::NewLine << "       D=" << D << " m  V=" << m_volume << " m^3"
                        << "  A_in=" << m_area[0] << " A_out=" << m_area[1] << " A_wall=" << m_area[2]
                        << " A_knife=" << m_area[3] << " m^2" << Alert::Raise;
          const Real imbalance = static_cast<Real>(maxCells) * ranks / static_cast<Real>(cells);
          Alert::Info() << "[partition] " << ranks << " ranks"
#ifdef RODIN_USE_SCOTCH
                        << " (Scotch)"
#else
                        << " (no Scotch: BalancedCompactPartitioner)"
#endif
                        << "  owned cells per rank: min " << minCells << " max " << maxCells
                        << "  max/mean = " << imbalance << Alert::Raise;
          if (ranks > 1 && imbalance > 1.5)
            Alert::Warning() << "[partition] badly balanced (max/mean = " << imbalance
                             << "): most of the work sits on few ranks. Build Rodin with "
                                "-DRODIN_USE_SCOTCH=ON." << Alert::Raise;
          Alert::Info() << "[fluid] rho=" << m_cfg.rho << " kg/m^3  K=" << m_cfg.K << " Pa s^n  n=" << m_cfg.n
                        << "  mu in [" << m_cfg.muMin << ", " << m_cfg.muMax << "] Pa s"
                        << Alert::NewLine << "[bc] inlet: "
                        << (m_velocityInlet
                              ? "velocity, U=" + std::to_string(m_cfg.uMean) + " m/s, " + m_profileName
                                + " profile, Re_MR=" + std::to_string(m_reynolds)
                              : "pressure, p_in=" + std::to_string(m_cfg.pIn) + " Pa ("
                                + (m_cfg.totalPressure ? "total" : "static") + ")")
                        << "  p_out=" << m_cfg.pOut << " Pa  ramp=" << m_cfg.ramp << " s" << Alert::Raise;
          m_csv.open(m_cfg.output + ".csv");
          m_csv << "t,dt,picard,Q_in,Q_out,p_in_mean,p_out_mean,"
                   "Fx,Fy,Fz,Fx_direct,Fy_direct,Fz_direct,Fx_pressure,Fy_pressure,Fz_pressure,"
                   "Fx_junction,Fy_junction,Fz_junction,Fx_dyn,Fy_dyn,Fz_dyn,zeta\n";
          m_csv << std::setprecision(10);
        }
      }

      int run()
      {
        const Real tAvg = (m_cfg.tAverage >= 0.0) ? m_cfg.tAverage : 0.5 * m_cfg.tEnd;
        int n = 0, quiet = 0;
        bool steady = false;
        while (m_t < m_cfg.tEnd - 1.0e-12 && !steady)
        {
          ++n;
          int iterations = 0;
          bool retried = false;
          for (int attempt = 0;; ++attempt)
          {
            m_dt = std::min(m_cfg.dt, m_cfg.tEnd - m_t);
            m_t += m_dt;
            bool threw = false;
            try
            {
              iterations = step();
            }
            catch (const std::runtime_error& e)
            {
              threw = true;
              if (isRoot())
                Alert::Warning() << e.what() << Alert::Raise;
            }
            if ((!threw && m_converged) || attempt >= m_cfg.maxRetries || 0.5 * m_cfg.dt < m_cfg.dtMin)
            {
              if (threw)
                throw std::runtime_error("step failed at t = " + std::to_string(m_t));
              if (!m_converged && isRoot())  // after all retries: accepted, flagged in the log
                Alert::Warning() << "  step accepted UNCONVERGED at t = " << m_t << Alert::Raise;
              break;
            }
            m_t -= m_dt;
            m_uIt.setData(m_uOld.getData());
            m_cfg.dt *= 0.5;
            if (isRoot())
              Alert::Info() << "  step rejected, retrying with dt = " << m_cfg.dt << Alert::Raise;
            retried = true;
          }

          // Forces at t^{n+1} need u^n: before the history update.
          Timer te(m_timing.forces);
          const auto r = evaluate();
          te.stop();
          Timer tp(m_timing.projections);
          m_uOld2.setData(m_uOld.getData());
          m_uOld.setData(m_u.getSolution().getData());
          m_dtPrev = m_dt;
          ++m_history;
          updateStabilisation();
          tp.stop();

          if (m_t >= tAvg - 1.0e-12)
          {
            m_avgTime += m_dt;
            for (size_t k = 0; k < 3; ++k)
            {
              m_avgF[k] += m_dt * r.F[k];
              m_avgFd[k] += m_dt * r.Fd[k];
            }
            m_avgQ += m_dt * r.Qin;
          }
          report(n, iterations, r);

          // Steady state: (F, Q) frozen for steadyWindow steps after the ramp.
          const Real scaleF = norm(r.F) + 1.0e-6 * 0.5 * m_cfg.rho * m_uRef * m_uRef * m_area[0];
          const Real scaleQ = std::abs(r.Qin) + 1.0e-9;
          const Real change = std::max(dist(r.F, m_lastF) / scaleF, std::abs(r.Qin - m_lastQ) / scaleQ);
          m_lastF = r.F;
          m_lastQ = r.Qin;
          quiet = (m_t > m_cfg.ramp && m_converged && change < m_cfg.steadyTolerance) ? quiet + 1 : 0;
          steady = quiet >= m_cfg.steadyWindow;

          const bool last = steady || !(m_t < m_cfg.tEnd - 1.0e-12);
          if (n % m_cfg.outputEvery == 0 || last)
            write();

          if (m_cfg.adaptiveDt)
          {
            // No growth during the ramp: the speed is still rising and a dt
            // grown on the easy first steps overshoots the CFL of the jet.
            Real dt = m_cfg.dt;
            if (!m_converged || iterations > m_cfg.picardHard)
              dt *= m_cfg.dtShrink;
            else if (iterations <= m_cfg.picardEasy && !retried && m_t >= m_cfg.ramp)
              dt *= m_cfg.dtGrow;
            m_cfg.dt = std::clamp(dt, m_cfg.dtMin, m_cfg.dtMax);
          }
        }
        m_xdmf.close();
        summary(steady, r_last());
        return 0;
      }

    private:
      /// Wall time per step by phase (root rank), reset after each report.
      struct Timing
      {
        Real setup = 0, coefficients = 0, assemble = 0, solve = 0, forces = 0, projections = 0, output = 0;
        int factorizations = 0, kspIterations = 0;
      };

      /// Adds the elapsed wall time to a Timing field on stop() or destruction.
      class Timer
      {
        public:
          explicit Timer(Real& acc) : m_acc(&acc), m_start(std::chrono::steady_clock::now()) {}
          ~Timer() { stop(); }
          void stop()
          {
            if (m_acc)
              *m_acc += std::chrono::duration<Real>(std::chrono::steady_clock::now() - m_start).count();
            m_acc = nullptr;
          }
        private:
          Real* m_acc;
          std::chrono::steady_clock::time_point m_start;
      };

      struct Result
      {
        Vec3 F{}, Fd{}, Fp{}, Fw{}, Fdyn{};
        Real Qin = 0.0, Qout = 0.0, pInMean = 0.0, pOutMean = 0.0, zeta = 0.0;
      };

      const Result& r_last() const { return m_last; }

      static Real norm(const Vec3& a) { return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]); }
      static Real dist(const Vec3& a, const Vec3& b)
      {
        return norm(Vec3{ a[0] - b[0], a[1] - b[1], a[2] - b[2] });
      }

      /// Smooth (C^1) ramp from 0 to 1 over [0, ramp].
      Real rampFactor() const
      {
        const Real s = (m_cfg.ramp > 0.0) ? std::min(m_t / m_cfg.ramp, 1.0) : 1.0;
        return 0.5 * (1.0 - std::cos(M_PI * s));
      }

      /// p_in(t): p_out + ramp (p_in - p_out).
      Real inletPressure() const
      {
        return m_cfg.pOut + rampFactor() * (m_cfg.pIn - m_cfg.pOut);
      }

      /// Velocity inlet: axis and centre of the inlet disc, profile shape and
      /// its normalisation so that the discrete flux (of the P1 interpolant
      /// the Dirichlet condition imposes) is exactly U A_in.
      void setupInflow()
      {
        const auto& L = m_cfg.labels;
        const auto nb = BoundaryNormal(m_mesh);
        for (size_t k = 0; k < 3; ++k)
        {
          const VectorFunction e{ Real(k == 0), Real(k == 1), Real(k == 2) };
          m_flowDir[k] = -boundaryIntegral(Dot(nb, e), L.inlet) / m_area[0];  // inward
          m_inletCentre[k] = boundaryIntegral(
              RealFunction([k](const Point& p) { return p.getCoordinates()(k); }), L.inlet) / m_area[0];
        }
        const Real d = norm(m_flowDir);
        for (auto& c : m_flowDir)
          c /= d;
        m_inletRadius = std::sqrt(m_area[0] / M_PI);

        // Metzner-Reed Reynolds number (= rho U D / mu for n = 1).
        const Real n = m_cfg.n, U = std::abs(m_cfg.uMean), D = 2.0 * m_inletRadius;
        m_reynolds = m_cfg.rho * std::pow(U, 2.0 - n) * std::pow(D, n)
                   / (m_cfg.K * std::pow(8.0, n - 1.0) * std::pow((3.0 * n + 1.0) / (4.0 * n), n));
        m_profileName = m_cfg.profile;
        if (m_profileName == "auto")
          m_profileName = (m_reynolds < 2100.0) ? "laminar" : "turbulent";
        if (m_profileName != "laminar" && m_profileName != "turbulent" && m_profileName != "uniform")
          throw std::runtime_error("-kv_profile must be auto, laminar, turbulent or uniform");

        ScalarGF shape(m_sh);
        shape.project(RealFunction([this](const Point& p) { return profileShape(p); }));
        m_profileScale = std::abs(m_cfg.uMean) * m_area[0] / boundaryIntegral(shape, L.inlet);
        if (m_cfg.uMean < 0.0)
          m_profileScale = -m_profileScale;
      }

      /// Profile shape w(r/R) at a point (unnormalised), zero at the wall.
      /// Laminar: fully developed power-law fluid, 1 - (r/R)^((n+1)/n)
      /// (Poiseuille for n = 1). Turbulent: (1 - r/R)^(1/7).
      Real profileShape(const Point& p) const
      {
        const auto x = p.getCoordinates();
        Real axial = 0.0, r2 = 0.0;
        for (size_t k = 0; k < 3; ++k)
          axial += (x(k) - m_inletCentre[k]) * m_flowDir[k];
        for (size_t k = 0; k < 3; ++k)
        {
          const Real dk = x(k) - m_inletCentre[k] - axial * m_flowDir[k];
          r2 += dk * dk;
        }
        const Real s = std::min(std::sqrt(r2) / m_inletRadius, 1.0);
        if (m_profileName == "laminar")
          return 1.0 - std::pow(s, (m_cfg.n + 1.0) / m_cfg.n);
        if (m_profileName == "turbulent")
          return std::pow(1.0 - s, 1.0 / 7.0);
        return s < 1.0 ? 1.0 : 0.0;
      }

      // ---- Picard iterations of one time step --------------------------------
      int step()
      {
        Timer timer(m_timing.setup);
        // Variable-step BDF2 (w = dt/dt_prev): (b0 u - b1 u^n + b2 u^{n-1})/dt,
        // b0 = (1 + 2w)/(1 + w), b1 = 1 + w, b2 = w^2/(1 + w); BDF1 without history.
        const bool history = m_history >= 1 && m_dtPrev > 0.0;
        const Real w = history ? m_dt / m_dtPrev : 0.0;
        if (m_cfg.bdf == 2 && history)
        {
          m_bdf[0] = (1.0 + 2.0 * w) / (1.0 + w);
          m_bdf[1] = 1.0 + w;
          m_bdf[2] = w * w / (1.0 + w);
        }
        else
        {
          m_bdf[0] = 1.0;
          m_bdf[1] = 1.0;
          m_bdf[2] = 0.0;
        }
        // Picard starts from the linear extrapolation of (u^{n-1}, u^n).
        m_uIt.setData(m_uOld.getData());
        if (m_cfg.extrapolate && history)
        {
          m_uIt.axpy(w, m_uOld);
          m_uIt.axpy(-w, m_uOld2);
        }
        setupFlow();  // holds 1/dt and p_in(t^{n+1})
        timer.stop();
        int k = 0;
        for (; k < m_cfg.maxIterations; ++k)
        {
          {
            Timer tc(m_timing.coefficients);
            updateCoefficients();  // mu, tau at u^k
          }
          solve();
          // Relative increment, with a velocity floor of 1e-3 U_B (U_B the
          // Bernoulli velocity of p_in - p_out): a closed gate has u ~ 0 and a
          // purely relative test would chase round-off.
          const auto [delta, scale] = increment(m_u.getSolution(), m_uIt);
          const Real dU = delta / std::max(scale, 1.0e-3 * m_velocityScale);
          m_uIt.setData(m_u.getSolution().getData());
          m_pIt.setData(m_p.getSolution().getData());
          if (m_cfg.picardMonitor && isRoot())
            Alert::Info() << "  Picard " << k + 1 << ": du/u = " << dU << Alert::Raise;
          if (dU < m_cfg.tolerance)
          {
            m_converged = true;
            return k + 1;
          }
        }
        m_converged = false;
        if (isRoot())
          Alert::Warning() << "Picard did not converge at t = " << m_t << " (dt = " << m_dt << ")"
                           << Alert::Raise;
        return k;
      }

      void setupFlow()
      {
        const auto& L = m_cfg.labels;
        const Real rho = m_cfg.rho, dt = m_dt;
        const auto n = BoundaryNormal(m_mesh);

        const auto convU = Mult(Jacobian(m_u), m_uIt);    // (grad u) u^k (Galerkin)
        const auto adv = VectorFunction(m_ax, m_ay, m_az);  // P0 mean of u^k (OSS)
        const auto convA = Mult(Jacobian(m_u), adv);
        const auto streamV = Mult(Jacobian(m_v), adv);
        const auto symU = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));
        const auto symV = 0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));
        const auto tangentialU = m_u - Dot(m_u, n) * n;
        const auto backflow = 0.5 * rho * Max(-Dot(m_uIt, n), 0.0);
        const Real pIn = inletPressure(), pOut = m_cfg.pOut;
        const Real totalPressure = m_cfg.totalPressure ? 1.0 : 0.0;

        const Real b0 = m_bdf[0], b1 = m_bdf[1], b2 = m_bdf[2];
        auto body = (rho * b0 / dt) * Integral(m_u, m_v)
               - (rho / dt) * Integral(b1 * m_uOld - b2 * m_uOld2, m_v)
               + rho * Integral(convU, m_v)
               + 2.0 * Integral(m_muGF * symU, symV)
               - Integral(m_p, Div(m_v))
               + Integral(Div(m_u), m_q)

               // Split OSS, projections at t^n.
               + rho * Integral(m_tauMGF * convA, streamV)
               - rho * Integral(m_tauMGF * m_piConv.get(), streamV)
               + Integral(m_tauMRhoGF * Grad(m_p), Grad(m_q))
               - Integral(m_tauMRhoGF * m_piGradP.get(), Grad(m_q))
               + Integral(m_tauCGF * Div(m_u), Div(m_v))
               - Integral(m_tauCGF * m_piDiv.get(), Div(m_v))

               // Free outlet: sigma n = -p_out n, directional do-nothing.
               + BoundaryIntegral(pOut * Dot(m_v, n)).over(L.outlet)
               + BoundaryIntegral(backflow * Dot(m_u, m_v)).over(L.outlet)

               + DirichletBC(m_u, Zero(3)).on(FlatSet<Attribute>{ L.wall, L.knife });

        if (m_velocityInlet)
        {
          // Developed profile, ramped in time, flux exactly ramp U A_in.
          const Real a = rampFactor() * m_profileScale;
          const auto inflow = VectorFunction(size_t{ 3 }, [this, a](const Point& p) {
            Math::SpatialVector<Real> u(3);
            const Real w = a * profileShape(p);
            for (size_t k = 0; k < 3; ++k)
              u(k) = w * m_flowDir[k];
            return u;
          });
          m_flow = body + DirichletBC(m_u, inflow).on(L.inlet);
        }
        else
        {
          // Pressure inlet: sigma n = -p_in n, u x n = 0 (penalty).
          m_flow = body
                 + BoundaryIntegral(pIn * Dot(m_v, n)).over(L.inlet)
                 + m_cfg.tangentialPenalty * BoundaryIntegral(Dot(tangentialU, m_v)).over(L.inlet)
                 + totalPressure * BoundaryIntegral(backflow * Dot(m_u, m_v)).over(L.inlet);
        }
      }

      void solve()
      {
        {
          Timer ta(m_timing.assemble);
          m_flow.assemble();
        }
        Timer ts(m_timing.solve);
        // Lagged LU: the MUMPS factors of an earlier matrix precondition
        // GMRES on the current one (Picard iterations and steps change the
        // matrix little). Refactor every refactorEvery solves, or as soon as
        // GMRES needs more than refactorIterations iterations or fails.
        ::KSP ksp = m_flowKSP.getHandle();
        const bool reuse = m_cfg.refactorEvery > 1 && m_factored
                         && m_solvesSinceFactor < m_cfg.refactorEvery;
        PetscErrorCode ierr = KSPSetReusePreconditioner(ksp, reuse ? PETSC_TRUE : PETSC_FALSE);
        assert(ierr == PETSC_SUCCESS);
        m_flow.solve(m_flowKSP);
        ::KSPConvergedReason reason;
        PetscInt its = 0;
        ierr = KSPGetConvergedReason(ksp, &reason);
        assert(ierr == PETSC_SUCCESS);
        ierr = KSPGetIterationNumber(ksp, &its);
        assert(ierr == PETSC_SUCCESS);
        if (reuse && (reason < 0 || its > m_cfg.refactorIterations))
        {
          ierr = KSPSetReusePreconditioner(ksp, PETSC_FALSE);
          assert(ierr == PETSC_SUCCESS);
          m_flow.solve(m_flowKSP);
          ierr = KSPGetConvergedReason(ksp, &reason);
          assert(ierr == PETSC_SUCCESS);
          ierr = KSPGetIterationNumber(ksp, &its);
          assert(ierr == PETSC_SUCCESS);
          m_solvesSinceFactor = 0;
          ++m_timing.factorizations;
        }
        else if (!reuse)
        {
          m_solvesSinceFactor = 0;
          ++m_timing.factorizations;
        }
        m_factored = true;
        ++m_solvesSinceFactor;
        m_timing.kspIterations += static_cast<int>(its);
        ts.stop();
        (void)ierr;
        if (reason < 0)
          throw std::runtime_error("KSP diverged at t = " + std::to_string(m_t));
      }

      // ---- Constitutive law and stabilisation ---------------------------------
      Real viscosity(Real gd) const
      {
        const Real mu = m_cfg.K * std::pow(std::max(gd, m_cfg.gammaMin), m_cfg.n - 1.0);
        return std::clamp(mu, m_cfg.muMin, m_cfg.muMax);
      }

      template <class Matrix>
      static Real shearRate(const Matrix& G)
      {
        Real s = 0.0;
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j)
          {
            const Real e = 0.5 * (G(i, j) + G(j, i));
            s += e * e;
          }
        return std::sqrt(2.0 * s);
      }

      /// Element size h = (6 sqrt(2) |T|)^(1/3): the edge of the regular
      /// tetrahedron of the same volume (the length the tau constants assume).
      static Real cellSize(const Point& p)
      {
        return std::cbrt(6.0 * std::sqrt(2.0) * p.getPolytope().getMeasure());
      }

      /// mu, tau_M, tau_M/rho, tau_C at the iterate u^k, constant per cell.
      void updateCoefficients()
      {
        const auto J = Jacobian(m_uIt);
        const Real rho = m_cfg.rho;
        const auto muAt = [this, &J](const Point& p) { return viscosity(shearRate(J.getValue(p))); };
        const auto tauAt = [this, &muAt, rho](const Point& p) {
          const Real h = cellSize(p);
          const auto u = m_uIt.getValue(p);
          const Real a = m_cfg.tauDt ? 2.0 / m_dt : 0.0;
          const Real b = 2.0 * std::sqrt(Math::dot(u, u)) / h;
          const Real c = 4.0 * muAt(p) / (rho * h * h);
          return 1.0 / std::sqrt(a * a + b * b + c * c);
        };
        updateAdvection();
        m_muGF.project(RealFunction(muAt));
        m_tauMGF.project(RealFunction(tauAt));
        m_tauMRhoGF.project(RealFunction([&tauAt, rho](const Point& p) { return tauAt(p) / rho; }));
        m_tauCGF.project(RealFunction([this, &muAt, rho](const Point& p) {
          const auto u = m_uIt.getValue(p);
          return muAt(p) + 0.5 * rho * std::sqrt(Math::dot(u, u)) * cellSize(p); }));
      }

      /// Advection velocity a of the OSS convective terms, constant per cell:
      /// the cell mean of u^k (its value at the centroid, exact for P1).
      /// tau_M rho ((grad u) a - Pi, (grad v) a) with a in P0 instead of u^k in
      /// P1: the same O(h) consistency, and the integrand no longer evaluates
      /// the P1 grid function u^k twice at every quadrature point and for
      /// every pair of basis functions, which made this one term ~30 % of the
      /// assembly. The Galerkin convection keeps u^k.
      void updateAdvection()
      {
        m_ax.project(RealFunction([this](const Point& p) { return m_uIt.getValue(p)(0); }));
        m_ay.project(RealFunction([this](const Point& p) { return m_uIt.getValue(p)(1); }));
        m_az.project(RealFunction([this](const Point& p) { return m_uIt.getValue(p)(2); }));
      }

      /// OSS projections Pi_h at the converged state of the step (t^n of the next).
      void updateStabilisation()
      {
        updateAdvection();
        const auto adv = VectorFunction(m_ax, m_ay, m_az);
        m_piConv.project(Mult(Jacobian(m_uIt), adv));
        m_piGradP.project(Grad(m_pIt));
        m_piDiv.project(Div(m_uIt));
      }

      // ---- Forces and diagnostics ---------------------------------------------
      Result evaluate()
      {
        const auto& L = m_cfg.labels;
        const auto& u = m_u.getSolution();
        const auto& p = m_p.getSolution();
        const auto n = BoundaryNormal(m_mesh);
        const Real rho = m_cfg.rho, dt = m_dt;
        Result r;

        // (1) Consistent reaction. R(w) = a((u_h, p_h), w): the momentum rows
        // of setupFlow() for every P1 vector test w, without boundary terms,
        // one assembly for the three components. Integrands: P1/P0 grid
        // functions and P1 test functions only.
        const auto J = Jacobian(u);
        const auto conv = Mult(J, u);
        const auto adv = VectorFunction(m_ax, m_ay, m_az);
        const auto convA = Mult(J, adv);
        const Real b0 = m_bdf[0], b1 = m_bdf[1], b2 = m_bdf[2];
        m_lfRes = Integral((rho / dt) * (b0 * u - b1 * m_uOld + b2 * m_uOld2), m_w)
                + Integral(rho * conv, m_w)
                + Integral(m_muGF * (J + Transpose(J)), Jacobian(m_w))  // 2 mu e(u) : e(w)
                + Integral(-1.0 * p, Div(m_w))
                + Integral(rho * m_tauMGF * (convA - m_piConv.get()), Mult(Jacobian(m_w), adv))
                + Integral(m_tauCGF * (Div(u) - m_piDiv.get()), Div(m_w));
        m_lfRes.assemble();
        // Pressure of the pipe wall on each node, int_wall p (n . w): the wall
        // share of the reaction at the knife-wall junction nodes.
        m_lfAux = BoundaryIntegral(p * Dot(n, m_w)).over(L.wall);
        m_lfAux.assemble();
        reaction(r);

        // (2) Direct surface integral, gradient from the incident cell (Point path).
        for (size_t k = 0; k < 3; ++k)
        {
          const VectorFunction e{ Real(k == 0), Real(k == 1), Real(k == 2) };
          r.Fd[k] = boundaryIntegral(RealFunction([this, k, &J, &n, &p](const Point& x) {
            const auto G = J.getValue(x);   // Point on a boundary face -> incident tetrahedron
            const auto nn = n.getValue(x);  // outward normal of the fluid
            const Real mu = viscosity(shearRate(G));
            Real sn = -p.getValue(x) * nn(k);
            for (size_t j = 0; j < 3; ++j)
              sn += mu * (G(k, j) + G(j, k)) * nn(j);
            return -sn;  // force on the knife = -(sigma n)
          }), L.knife);
          r.Fp[k] = boundaryIntegral(p * Dot(n, e), L.knife);
          r.Fdyn[k] = r.F[k] - m_cfg.pOut * m_knifeNormal[k];
        }

        r.Qin = -boundaryIntegral(Dot(u, n), L.inlet);
        r.Qout = boundaryIntegral(Dot(u, n), L.outlet);
        r.pInMean = boundaryIntegral(p, L.inlet) / m_area[0];
        r.pOutMean = boundaryIntegral(p, L.outlet) / m_area[1];
        const Real U = r.Qin / m_area[0];
        r.zeta = (std::abs(U) > 1.0e-6 * m_velocityScale / std::sqrt(static_cast<Real>(m_vh.getSize()))) ? (r.pInMean - r.pOutMean) / (0.5 * rho * U * U) : 0.0;
        m_last = r;
        return r;
      }

      /// From the assembled residual R (m_lfRes, one entry per velocity dof
      /// (node i, component k)) and the wall pressure P_k,i = int_wall p n_k
      /// phi_i (m_lfAux):
      ///
      ///   F_k = -sum_{i on the knife} (R_k,i + P_k,i).
      ///
      /// The sum runs over every node of the knife, junction with the pipe
      /// included. At a junction node the reaction R_k,i also holds the pipe
      /// wall faces around it; their pressure share, -P_k,i, is removed
      /// exactly (P1 p, same faces), so what is left of the wall is its shear
      /// on a strip one cell wide. (Weighting R_k,i by A_i^knife/A_i instead
      /// loses the knife part of the junction reaction: -1.7% of the thrust of
      /// the closed FPC8 on the coarse mesh.) The component of each dof comes
      /// from the interpolants of e_k (m_comp). Also writes the nodal wall load
      /// -R_k,i/A_i (lumped, Pa) into m_load, for visualisation only.
      void reaction(Result& out)
      {
        const ::Vec R = m_lfRes.getVector(), A = m_lfArea.getVector();
        const ::Vec Ak = m_lfKnife.getVector(), P = m_lfAux.getVector();
        ::Vec t;
        PetscErrorCode ierr = VecDuplicate(R, &t);
        assert(ierr == PETSC_SUCCESS);
        PetscInt nloc = 0;
        ierr = VecGetLocalSize(R, &nloc);
        assert(ierr == PETSC_SUCCESS);
        const PetscScalar *r, *a, *ak, *pw, *e[3];
        PetscScalar* tt;
        VecGetArrayRead(R, &r);
        VecGetArrayRead(A, &a);
        VecGetArrayRead(Ak, &ak);
        VecGetArrayRead(P, &pw);
        for (size_t k = 0; k < 3; ++k)
          VecGetArrayRead(m_comp[k].getData(), &e[k]);
        VecGetArray(t, &tt);
        Real local[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
        for (PetscInt i = 0; i < nloc; ++i)
        {
          tt[i] = (a[i] > 0.0) ? -r[i] / a[i] : 0.0;  // interior nodes: no load
          if (ak[i] > 0.0)                             // dof of a knife node
            for (size_t k = 0; k < 3; ++k)
            {
              local[k] -= e[k][i] * (r[i] + pw[i]);
              local[3 + k] += e[k][i] * pw[i];
            }
        }
        VecRestoreArray(t, &tt);
        for (size_t k = 0; k < 3; ++k)
          VecRestoreArrayRead(m_comp[k].getData(), &e[k]);
        VecRestoreArrayRead(P, &pw);
        VecRestoreArrayRead(Ak, &ak);
        VecRestoreArrayRead(A, &a);
        VecRestoreArrayRead(R, &r);
        Real global[6];
        MPI_Allreduce(local, global, 6, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        for (size_t k = 0; k < 3; ++k)
        {
          out.F[k] = global[k];
          out.Fw[k] = global[3 + k];
        }
        m_load.setData(t);
        ierr = VecDestroy(&t);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
      }

      void report(int step, int iterations, const Result& r)
      {
        if (!isRoot())
          return;
        m_csv << m_t << ',' << m_dt << ',' << iterations << ',' << r.Qin << ',' << r.Qout << ','
              << r.pInMean << ',' << r.pOutMean;
        for (const auto* f : { &r.F, &r.Fd, &r.Fp, &r.Fw, &r.Fdyn })
          for (size_t k = 0; k < 3; ++k)
            m_csv << ',' << (*f)[k];
        m_csv << ',' << r.zeta << '\n';
        m_csv.flush();
        Alert::Info() << "step " << step << "  t=" << m_t << " s  dt=" << m_dt << "  picard=" << iterations
                      << (m_converged ? "" : " (NOT converged)")
                      << "  Q=" << r.Qin << " m^3/s (out " << r.Qout << ")  p_in=" << r.pInMean
                      << " Pa" << Alert::NewLine
                      << "   F_knife = (" << r.F[0] << ", " << r.F[1] << ", " << r.F[2] << ") N"
                      << "   direct = (" << r.Fd[0] << ", " << r.Fd[1] << ", " << r.Fd[2] << ") N"
                      << "   pressure = (" << r.Fp[0] << ", " << r.Fp[1] << ", " << r.Fp[2] << ") N"
                      << Alert::NewLine
                      << "   F_dyn (F - p_out int n) = (" << r.Fdyn[0] << ", " << r.Fdyn[1] << ", "
                      << r.Fdyn[2] << ") N" << Alert::NewLine
                      << "   time [s]: setup " << m_timing.setup << "  coefficients " << m_timing.coefficients
                      << "  assemble " << m_timing.assemble << "  solve " << m_timing.solve
                      << " (" << m_timing.factorizations << " LU, " << m_timing.kspIterations << " GMRES its)"
                      << "  forces " << m_timing.forces << "  projections " << m_timing.projections
                      << "  output " << m_timing.output << Alert::Raise;
        m_timing = Timing{};
      }

      void summary(bool steady, const Result& r)
      {
        if (!isRoot())
          return;
        Alert::Info() << (steady ? "Steady state reached" : "End of the time interval") << " at t=" << m_t
                      << " s" << Alert::NewLine
                      << "  last:     F_knife = (" << r.F[0] << ", " << r.F[1] << ", " << r.F[2] << ") N"
                      << "  Q = " << r.Qin << " m^3/s  zeta = " << r.zeta << Alert::Raise;
        if (m_avgTime > 0.0)
          Alert::Info() << "  average over the last " << m_avgTime << " s:  F_knife = ("
                        << m_avgF[0] / m_avgTime << ", " << m_avgF[1] / m_avgTime << ", "
                        << m_avgF[2] / m_avgTime << ") N   direct = (" << m_avgFd[0] / m_avgTime << ", "
                        << m_avgFd[1] / m_avgTime << ", " << m_avgFd[2] / m_avgTime << ") N   Q = "
                        << m_avgQ / m_avgTime << " m^3/s" << Alert::NewLine
                        << "  F_dyn = (" << m_avgF[0] / m_avgTime - m_cfg.pOut * m_knifeNormal[0] << ", "
                        << m_avgF[1] / m_avgTime - m_cfg.pOut * m_knifeNormal[1] << ", "
                        << m_avgF[2] / m_avgTime - m_cfg.pOut * m_knifeNormal[2] << ") N" << Alert::Raise;
      }

      void write()
      {
        Timer tw(m_timing.output);
        m_muOut.project(m_muGF);
        m_xdmf.write(m_t).flush();
      }

      // ---- Utilities ------------------------------------------------------------
      /// (||current - previous||, ||current||), l2 norms of the dof vectors.
      template <class GF>
      static std::pair<Real, Real> increment(const GF& current, const GF& previous)
      {
        GF diff(current.getFiniteElementSpace());
        diff.setData(current.getData());
        PetscReal scale = 0.0, delta = 0.0;
        PetscErrorCode ierr = VecAXPY(diff.getData(), -1.0, previous.getData());
        assert(ierr == PETSC_SUCCESS);
        ierr = VecNorm(current.getData(), NORM_2, &scale);
        assert(ierr == PETSC_SUCCESS);
        ierr = VecNorm(diff.getData(), NORM_2, &delta);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
        return { delta, scale };
      }

      template <class Expression>
      Real boundaryIntegral(const Expression& f, Attribute a)
      {
        m_flux = BoundaryIntegral(f, m_z).over(a);
        m_flux.assemble();
        return m_flux(m_one);
      }

      template <class Expression>
      Real volumeIntegral(const Expression& f)
      {
        m_flux = Integral(f, m_z);
        m_flux.assemble();
        return m_flux(m_one);
      }

      template <class M>
      static void connect(M& mesh)
      {
        const size_t D = mesh.getDimension();
        mesh.getConnectivity().compute(D, D);
        mesh.getConnectivity().compute(D, 0);
        mesh.getConnectivity().compute(D, D - 1);
        mesh.getConnectivity().compute(D - 1, D);
        mesh.getConnectivity().compute(D - 1, 0);
        mesh.getConnectivity().compute(D - 1, 1);
        mesh.getConnectivity().compute(1, 0);
      }

      static MeshType makeMesh(const Context::MPI& context, const Config& cfg)
      {
        const auto& comm = context.getCommunicator();
        Rodin::MPI::Sharder sharder(context);
        if (comm.rank() == 0)
        {
          Mesh<Context::Local> mesh;
          mesh.load(cfg.meshPath, IO::FileFormat::MEDIT);
          if (mesh.getSpaceDimension() != 3 || mesh.getDimension() != 3)
            throw std::runtime_error("KnifeGateValve expects a tetrahedral mesh.");
          connect(mesh);
#ifdef RODIN_USE_SCOTCH
          Scotch::Partitioner partitioner(mesh);
#else
          BalancedCompactPartitioner partitioner(mesh);
#endif
          partitioner.partition(static_cast<size_t>(comm.size()));
          sharder.shard(partitioner);
          sharder.scatter(0);
        }
        MeshType mesh = sharder.gather(0);
        if (cfg.meshScale != 1.0)
          mesh.scale(cfg.meshScale);
        connect(mesh);
        mesh.reconcile(2);
        mesh.reconcile(1);
        return mesh;
      }

      bool isRoot() const { return m_mesh.getContext().getCommunicator().rank() == 0; }

      Config m_cfg;
      MeshType m_mesh;
      IO::XDMF m_xdmf;
      VectorFES m_vh;
      ScalarFES m_sh;
      CellFES m_ch;

      VectorTrial m_u;
      ScalarTrial m_p;
      VectorTest m_v;
      ScalarTest m_q;

      VectorGF m_uOld, m_uOld2, m_uIt;  ///< time levels n, n-1, Picard iterate k
      ScalarGF m_pIt;
      CellGF m_muGF, m_tauMGF, m_tauMRhoGF, m_tauCGF;
      CellGF m_ax, m_ay, m_az;  ///< advection velocity of the OSS terms: cell mean of u^k
      VectorGF m_load;     ///< wall load -sigma n (Pa), consistent reaction (lumped)
      VectorGF m_comp[3];  ///< interpolants of e_k: component of each velocity dof

      OrthogonalProjection<VectorFES> m_piConv, m_piGradP;
      OrthogonalProjection<ScalarFES> m_piDiv, m_muOut;

      FlowProblem m_flow;
      Solver::KSP m_flowKSP;

      ScalarTest m_z;
      ScalarGF m_one;
      VectorTest m_w;  ///< test of the reaction forms
      ScalarLF m_flux;
      VectorLF m_lfArea, m_lfKnife, m_lfRes, m_lfAux;

      Real m_t = 0.0, m_dt = 0.0;
      Real m_dtPrev = 0.0;                  ///< dt of the last accepted step
      bool m_factored = false;              ///< MUMPS factors available for reuse
      int m_solvesSinceFactor = 0;
      int m_history = 0;                    ///< accepted steps (BDF2 and extrapolation need 1)
      Real m_bdf[3] = { 1.0, 1.0, 0.0 };    ///< (b0 u - b1 u^n + b2 u^{n-1})/dt of this step
      bool m_converged = true;
      Real m_area[4] = { 0.0, 0.0, 0.0, 0.0 };
      Real m_volume = 0.0;
      Real m_velocityScale = 1.0;
      bool m_velocityInlet = true;
      Real m_uRef = 1.0;
      Vec3 m_knifeNormal{};
      Vec3 m_flowDir{}, m_inletCentre{};
      Real m_inletRadius = 0.0, m_profileScale = 0.0, m_reynolds = 0.0;
      std::string m_profileName = "none";
      Vec3 m_lastF{}, m_avgF{}, m_avgFd{};
      Real m_lastQ = 0.0, m_avgQ = 0.0, m_avgTime = 0.0;
      Result m_last;
      Timing m_timing;
      std::ofstream m_csv;
  };
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  const auto setDefault = [](const std::string& name, const char* value) {
    PetscBool set = PETSC_FALSE;
    PetscErrorCode ierr = PetscOptionsHasName(PETSC_NULLPTR, PETSC_NULLPTR, name.c_str(), &set);
    if (ierr == PETSC_SUCCESS && !set)
      ierr = PetscOptionsSetValue(PETSC_NULLPTR, name.c_str(), value);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
  };
  // Flow: direct (MUMPS). This MUMPS build faults in its distributed RHS
  // scatter: keep the right-hand side and the solution centralised.
  // GMRES preconditioned by (possibly lagged) MUMPS factors; with fresh
  // factors it converges in one iteration. -kv_refactor_every 1 together
  // with -kv_flow_ksp_type preonly is the plain direct solve.
  setDefault("-kv_flow_ksp_type", "gmres");
  setDefault("-kv_flow_ksp_rtol", "1e-9");
  setDefault("-kv_flow_ksp_atol", "1e-14");
  setDefault("-kv_flow_ksp_max_it", "40");
  setDefault("-kv_flow_ksp_gmres_restart", "40");
  setDefault("-kv_flow_pc_type", "lu");
  setDefault("-kv_flow_pc_factor_mat_solver_type", "mumps");
  setDefault("-kv_flow_mat_mumps_icntl_20", "0");
  setDefault("-kv_flow_mat_mumps_icntl_21", "0");
  setDefault("-kv_flow_mat_mumps_icntl_14", "50");
  for (const char* prefix : { "-kv_proj_conv_", "-kv_proj_gradp_", "-kv_proj_div_", "-kv_proj_mu_" })
  {
    setDefault(std::string(prefix) + "ksp_type", "cg");
    setDefault(std::string(prefix) + "pc_type", "jacobi");
    setDefault(std::string(prefix) + "ksp_rtol", "1e-10");
  }

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world(PETSC_COMM_WORLD, boost::mpi::comm_attach);
  Rodin::Context::MPI context(env, world);

  int status = 0;
  try
  {
    using Rodin::Examples::Valve::KnifeGateValve;
    KnifeGateValve::Config cfg;

    char buffer[PETSC_MAX_PATH_LEN];
    PetscBool got = PETSC_FALSE;
    PetscOptionsGetString(PETSC_NULLPTR, PETSC_NULLPTR, "-kv_mesh", buffer, sizeof(buffer), &got);
    if (got)
      cfg.meshPath = buffer;
    PetscOptionsGetString(PETSC_NULLPTR, PETSC_NULLPTR, "-kv_out", buffer, sizeof(buffer), &got);
    if (got)
      cfg.output = buffer;

    const auto real = [](const char* name, Rodin::Real& value) {
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, name, &value, PETSC_NULLPTR);
    };
    const auto integer = [](const char* name, int& value) {
      PetscInt v = value;
      PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, name, &v, PETSC_NULLPTR);
      value = static_cast<int>(v);
    };
    const auto boolean = [](const char* name, bool& value) {
      PetscBool v = value ? PETSC_TRUE : PETSC_FALSE;
      PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, name, &v, PETSC_NULLPTR);
      value = (v == PETSC_TRUE);
    };
    real("-kv_scale", cfg.meshScale);
    real("-kv_rho", cfg.rho);
    real("-kv_K", cfg.K);
    real("-kv_n", cfg.n);
    real("-kv_gdmin", cfg.gammaMin);
    real("-kv_mumin", cfg.muMin);
    real("-kv_mumax", cfg.muMax);
    PetscOptionsGetString(PETSC_NULLPTR, PETSC_NULLPTR, "-kv_inlet", buffer, sizeof(buffer), &got);
    if (got)
      cfg.inlet = buffer;
    PetscOptionsGetString(PETSC_NULLPTR, PETSC_NULLPTR, "-kv_profile", buffer, sizeof(buffer), &got);
    if (got)
      cfg.profile = buffer;
    real("-kv_umean", cfg.uMean);
    real("-kv_pin", cfg.pIn);
    real("-kv_pout", cfg.pOut);
    real("-kv_ramp", cfg.ramp);
    real("-kv_tangential", cfg.tangentialPenalty);
    boolean("-kv_total_pressure", cfg.totalPressure);
    real("-kv_dt", cfg.dt);
    real("-kv_dtmin", cfg.dtMin);
    real("-kv_dtmax", cfg.dtMax);
    real("-kv_tend", cfg.tEnd);
    boolean("-kv_adaptive_dt", cfg.adaptiveDt);
    real("-kv_grow", cfg.dtGrow);
    real("-kv_shrink", cfg.dtShrink);
    integer("-kv_easy", cfg.picardEasy);
    integer("-kv_hard", cfg.picardHard);
    integer("-kv_maxit", cfg.maxIterations);
    real("-kv_tol", cfg.tolerance);
    boolean("-kv_tau_dt", cfg.tauDt);
    integer("-kv_retries", cfg.maxRetries);
    boolean("-kv_picard_monitor", cfg.picardMonitor);
    integer("-kv_bdf", cfg.bdf);
    boolean("-kv_extrapolate", cfg.extrapolate);
    boolean("-kv_lumped", cfg.lumpedProjection);
    integer("-kv_refactor_every", cfg.refactorEvery);
    integer("-kv_refactor_its", cfg.refactorIterations);
    real("-kv_tavg", cfg.tAverage);
    real("-kv_steady_tol", cfg.steadyTolerance);
    integer("-kv_steady_window", cfg.steadyWindow);
    integer("-kv_output_every", cfg.outputEvery);
    cfg.outputEvery = std::max(cfg.outputEvery, 1);

    KnifeGateValve simulation(context, cfg);
    status = simulation.run();
  }
  catch (const std::exception& e)
  {
    std::cerr << "Fatal error: " << e.what() << "\n";
    status = 1;
  }
  PetscFinalize();
  return status;
}
