// ForcedConvection.cpp
//
// Flow boiling of water in a 3D channel with the Lee phase-change model.
// One-fluid formulation with vapour fraction phi in [0, 1]:
//
//   d phi/dt + div(phi u) = gamma div(eps grad phi - phi (1 - phi) n) + mdot / rho_v
//   div u                                         = mdot (1/rho_v - 1/rho_l)
//   rho (u - u^n)/dt + rho (grad u) u - div(2 mu e(u)) + grad p
//                                                 = rho g + sigma kappa grad phi
//   rho cp (dT/dt + u . grad T) - div(k grad T)   = - mdot h_lv
//
// with mixture properties phi-weighted (arithmetic for rho, mu, rho cp;
// harmonic for k) and the Lee source
//
//   mdot = beta(phi, T) (T - T_sat)/T_sat,
//   beta = r_l (1 - phi) rho_l   (T > T_sat, evaporation)
//        = r_v phi rho_v         (T < T_sat, condensation),
//
// so mdot is continuous at T_sat. The switch is smoothed with tanh and the
// source is linearised in (phi, T) about the previous Picard iterate: both
// resulting reaction terms are positive on the left-hand side.
//
// The conservative Allen-Cahn term (Chiu & Lin 2011) keeps the interface at a
// fixed thickness eps (tanh profile) and phi in [0, 1]; n = Pi_h[grad phi /
// |grad phi|] is shared with the curvature. The compression is semi-implicit:
// gamma ((1 - phi^k) n^k phi, grad a). gamma is a velocity scale (default: the
// peak inlet velocity) and eps defaults to 1.5 mean cell sizes.
//
// Discretisation: backward Euler. (phi, T) are solved monolithically in P2/P2,
// (u, p) monolithically in equal-order P1/P1. Both blocks are stabilised with
// split orthogonal subgrid scales (OSS): tau (A(u_h) - Pi_h[A(u_h^n)], B(v_h))
// with Pi_h the L2 projection from OrthogonalProjection.h, lagged at t^n. The
// two blocks are coupled by Picard iterations within each time step
// (properties, mdot, normal and curvature refreshed every iteration).
// Curvature is kappa = -div Pi_h[grad phi / |grad phi|].
//
// Run from the build directory:
//   mpirun -n 4 ./examples/Heat/ForcedConvection -fc_dt 1e-5 -fc_tend 0.02 -fc_tin 370
//
// Options: -fc_mesh, -fc_dt, -fc_tend, -fc_G, -fc_q, -fc_tin, -fc_r, -fc_g,
// -fc_eps, -fc_gamma (0 disables Allen-Cahn), -fc_maxit, -fc_tol, -fc_output_every,
// -fc_print_every (progress line on screen, default every step), -fc_picard_monitor,
// -fc_neps (normal regularisation, 1/m; default 0.05/eps); solvers under the prefixes -fc_flow_
// (MUMPS), -fc_phase_ (GMRES) and -fc_proj_* (CG + Jacobi mass solves).
#include "Rodin/Variational/ForwardDecls.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

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

namespace Rodin::Examples::Heat
{
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  class ForcedConvection
  {
    public:
      using MeshType = Mesh<Context::MPI>;
      using VectorFES = H1<1, Math::SpatialVector<Real>, MeshType>;
      using ScalarFES = H1<1, Real, MeshType>;
      using ScalarFES2 = H1<2, Real, MeshType>;
      using VectorGF = PETSc::Variational::GridFunction<VectorFES>;
      using ScalarGF = PETSc::Variational::GridFunction<ScalarFES>;
      using ScalarGF2 = PETSc::Variational::GridFunction<ScalarFES2>;
      using VectorTrial = PETSc::Variational::TrialFunction<VectorGF, VectorFES>;
      using ScalarTrial = PETSc::Variational::TrialFunction<ScalarGF, ScalarFES>;
      using ScalarTrial2 = PETSc::Variational::TrialFunction<ScalarGF2, ScalarFES2>;
      using VectorTest = PETSc::Variational::TestFunction<VectorFES>;
      using ScalarTest = PETSc::Variational::TestFunction<ScalarFES>;
      using ScalarTest2 = PETSc::Variational::TestFunction<ScalarFES2>;
      using System = PETSc::Math::LinearSystem;
      using FlowProblem = Problem<System, VectorTrial, ScalarTrial, VectorTest, ScalarTest>;
      using PhaseProblem = Problem<System, ScalarTrial2, ScalarTrial2, ScalarTest2, ScalarTest2>;
      using Coefficient = RealFunction<std::function<Real(const Point&)>>;

      struct Labels
      {
        Attribute inlet = 1, outlet = 2, bottom = 3, top = 4, sides = 5;
      };

      struct Config
      {
        std::string meshPath = "../resources/examples/Heat/DivergingMicrochannel3D.mesh";
        std::string output = "ForcedConvection";
        Labels labels;

        // Saturated water at p = 101.325 kPa (steam tables).
        Real rhoL = 958.4;     ///< kg/m^3
        Real rhoV = 0.5978;    ///< kg/m^3
        Real muL = 2.82e-4;    ///< Pa s
        Real muV = 1.23e-5;    ///< Pa s
        Real kL = 0.679;       ///< W/(m K)
        Real kV = 0.0251;      ///< W/(m K)
        Real cpL = 4216.0;     ///< J/(kg K)
        Real cpV = 2080.0;     ///< J/(kg K)
        Real hLV = 2.257e6;    ///< J/kg
        Real sigma = 0.0589;   ///< N/m
        Real Tsat = 373.15;    ///< K

        // Lee model. r is not physical: calibrate on the 1D Stefan problem.
        Real rL = 100.0;       ///< 1/s, evaporation
        Real rV = 100.0;       ///< 1/s, condensation
        Real dTsmooth = 0.05;  ///< K, half-width of the tanh switch at T_sat
        Real gravity = 0.0;    ///< m/s^2 along -z; 0 keeps the microchannel horizontal
        /// 1/m, regularisation of n = grad phi / sqrt(|grad phi|^2 + neps^2).
        /// Must be comparable to interface gradients (~1/(4 eps)), not ~0: with
        /// neps -> 0, n is a unit vector even where phi ~ 1e-5, and the
        /// compression gamma div(phi n) turns into an anti-diffusive reaction of
        /// rate ~gamma |div n| ~ gamma/h. Negative = automatic: 0.05/eps.
        Real normalEps = -1.0;

        // Conservative Allen-Cahn interface regularisation. Negative = automatic.
        Real epsilon = -1.0;   ///< m, interface half-thickness (auto: 1.5 h_mean)
        Real gamma = -1.0;     ///< m/s, mobility (auto: peak inlet velocity); 0 disables

        Real massFlux = 240.0;          ///< G (kg/m^2 s); mean inlet velocity G/rho_l
        Real inletWidth = 2.0e-4;       ///< m, inlet centred on y = 0
        Real depth = 4.0e-4;            ///< m, z in [0, depth]
        Real inletTemperature = 300.0;  ///< K
        Real heatFlux = 1.5e5;          ///< W/m^2 on the bottom
        Real outletPressure = 0.0;      ///< Pa (gauge)

        Real dt = 1.0e-5;
        Real tEnd = 0.02;
        Real ramp = 5.0e-3;  ///< s, smooth start of the inflow
        int maxIterations = 10;   ///< Picard iterations per step
        Real tolerance = 1.0e-4;  ///< relative increment of (u, phi, T)
        int outputEvery = 20;
        int printEvery = 1;       ///< steps between progress lines on screen
        bool picardMonitor = false; ///< print (dPhi, dT, dU) at every Picard iteration
      };

      ForcedConvection(const Context::MPI& context, const Config& cfg)
        : m_cfg(cfg),
          m_mesh(makeMesh(context, m_cfg.meshPath)),
          m_xdmf(context.getCommunicator(), m_cfg.output),
          m_vh(std::integral_constant<size_t, 1>{}, m_mesh, m_mesh.getSpaceDimension()),
          m_sh(std::integral_constant<size_t, 1>{}, m_mesh),
          m_sh2(std::integral_constant<size_t, 2>{}, m_mesh),
          m_u(m_vh), m_p(m_sh), m_v(m_vh), m_q(m_sh),
          m_phi(m_sh2), m_T(m_sh2), m_a(m_sh2), m_w(m_sh2),
          m_uOld(m_vh), m_uIt(m_vh), m_pIt(m_sh),
          m_phiOld(m_sh2), m_TOld(m_sh2), m_phiIt(m_sh2), m_TIt(m_sh2),
          m_phiOut(m_sh), m_TOut(m_sh),
          m_rho([this](const Point& p) { return rhoAt(p); }),
          m_mu([this](const Point& p) { return muAt(p); }),
          m_k([this](const Point& p) { return kAt(p); }),
          m_rhoCp([this](const Point& p) { return rhoCpAt(p); }),
          m_liquidFraction([this](const Point& p) { return 1.0 - phiAt(p); }),
          m_cT([this](const Point& p) { return betaAt(p) / m_cfg.Tsat; }),
          m_cPhi([this](const Point& p) { return dbetaAt(p) * thetaAt(p); }),
          m_c0([this](const Point& p) {
            return -betaAt(p) - dbetaAt(p) * thetaAt(p) * phiAt(p); }),
          m_divSource([this](const Point& p) {
            return betaAt(p) * thetaAt(p) * (1.0 / m_cfg.rhoV - 1.0 / m_cfg.rhoL); }),
          m_tauM([this](const Point& p) { return tau(p, muAt(p) / rhoAt(p), 0.0); }),
          m_tauT([this](const Point& p) {
            const Real C = rhoCpAt(p);
            return tau(p, kAt(p) / C, m_cfg.hLV * betaAt(p) / (m_cfg.Tsat * C)); }),
          m_tauPhi([this](const Point& p) {
            return tau(p, m_gamma * m_eps, -dbetaAt(p) * thetaAt(p) / m_cfg.rhoV); }),
          m_tauC([this](const Point& p) {
            return muAt(p) + 0.5 * rhoAt(p) * speed(p) * cellSize(p); }),
          m_piConv(m_vh, "fc_proj_conv_"), m_piGradP(m_vh, "fc_proj_gradp_"),
          m_piNormal(m_vh, "fc_proj_normal_"), m_piDiv(m_sh, "fc_proj_div_"),
          m_piKappa(m_sh, "fc_proj_kappa_"),
          m_piConvPhi(m_sh2, "fc_proj_phi_"), m_piConvT(m_sh2, "fc_proj_T_"),
          m_flow(m_u, m_p, m_v, m_q), m_flowKSP(m_flow),
          m_phase(m_phi, m_T, m_a, m_w), m_phaseKSP(m_phase),
          m_z(m_sh), m_one(m_sh), m_flux(m_z)
      {
        const size_t dim = m_mesh.getSpaceDimension();
        Math::SpatialVector<Real> zero(dim);
        for (size_t i = 0; i < dim; ++i)
          zero(i) = 0.0;
        m_uOld = zero;
        m_uIt = zero;
        m_pIt = 0.0;
        m_TOld = m_cfg.inletTemperature;
        m_TIt = m_cfg.inletTemperature;
        m_phiOld = 0.0;  // all liquid
        m_phiIt = 0.0;
        m_one = Real(1);

        m_flowKSP.setPrefix("fc_flow_");
        m_phaseKSP.setPrefix("fc_phase_");

        const auto& L = m_cfg.labels;
        m_bottomArea = boundaryIntegral(m_one, L.bottom);
        m_inletArea = boundaryIntegral(m_one, L.inlet);
        m_outletArea = boundaryIntegral(m_one, L.outlet);

        // Allen-Cahn parameters from the mesh and the inflow when not given.
        const Real hMean = std::cbrt(volumeIntegral(m_one) / static_cast<Real>(m_mesh.getCellCount()));
        m_eps = (m_cfg.epsilon > 0.0) ? m_cfg.epsilon : 1.5 * hMean;
        m_normalEps = (m_cfg.normalEps > 0.0) ? m_cfg.normalEps : 0.05 / m_eps;
        m_gamma = (m_cfg.gamma >= 0.0) ? m_cfg.gamma : 2.25 * m_cfg.massFlux / m_cfg.rhoL;

        // Projections must hold consistent data before the forms are assembled.
        updateStabilisation();
        updateInterface();

        setupFlow();
        setupPhase();

        m_xdmf.setMesh(m_mesh);
        m_xdmf.add("velocity", m_u.getSolution());
        m_xdmf.add("pressure", m_p.getSolution());
        m_xdmf.add("temperature", m_TOut);  // P2 fields interpolated to P1 for output
        m_xdmf.add("fraction", m_phiOut);
        m_xdmf.add("curvature", m_piKappa.get());

        const auto cells = m_mesh.getCellCount(), vertices = m_mesh.getVertexCount();  // collective
        if (isRoot())
        {
          m_csv.open(m_cfg.output + ".csv");
          m_csv << "t,iterations,maxU,maxT,minPhi,maxPhi,vaporVolume,bulkTout,dp\n";
          Alert::Info() << "[mesh] cells=" << cells << " vertices=" << vertices
                        << "  bottom area=" << m_bottomArea << " m^2"
                        << "  Allen-Cahn eps=" << m_eps << " m gamma=" << m_gamma << " m/s"
                        << " neps=" << m_normalEps << " 1/m"
                        << Alert::Raise;
        }
      }

      int run()
      {
        const int steps = static_cast<int>(m_cfg.tEnd / m_cfg.dt + 0.5);
        for (int n = 1; n <= steps; ++n)
        {
          m_t = n * m_cfg.dt;
          const int iterations = step();

          m_uOld.setData(m_u.getSolution().getData());
          m_phiOld.setData(m_phi.getSolution().getData());
          m_TOld.setData(m_T.getSolution().getData());
          updateStabilisation();  // Pi_h at t^n for the next step

          report(n, iterations);
          if (n % m_cfg.outputEvery == 0 || n == steps)
          {
            m_phiOut = m_phi.getSolution();
            m_TOut = m_T.getSolution();
            m_xdmf.write(m_t).flush();
          }
        }
        m_xdmf.close();
        return 0;
      }

    private:
      /// One time step: Picard iterations phase -> flow until the relative
      /// increments of (phi, T, u) fall below the tolerance.
      int step()
      {
        int k = 0;
        for (; k < m_cfg.maxIterations; ++k)
        {
          // (phi, T) monolithic, with properties, mdot and Pi_h at the iterate.
          solve(m_phase, m_phaseKSP);
          const Real dPhi = increment(m_phi.getSolution(), m_phiIt);
          const Real dT = increment(m_T.getSolution(), m_TIt);
          m_phiIt.setData(m_phi.getSolution().getData());
          m_TIt.setData(m_T.getSolution().getData());

          // (u, p) with the fresh fraction: density, viscosity, capillary force
          // and continuity source see (phi, T)^{k+1}. The OSS projections stay
          // at t^n: lagged at the iterate they make Picard contract only like
          // tau_M/(tau_M + dt) on smooth pressure modes (~0.37 measured).
          updateInterface();
          solve(m_flow, m_flowKSP);
          const Real dU = increment(m_u.getSolution(), m_uIt);
          m_uIt.setData(m_u.getSolution().getData());
          m_pIt.setData(m_p.getSolution().getData());

          if (m_cfg.picardMonitor && isRoot())
            Alert::Info() << "  Picard " << k + 1 << ": dPhi=" << dPhi << " dT=" << dT
                          << " dU=" << dU << Alert::Raise;

          if (std::max({ dPhi, dT, dU }) < m_cfg.tolerance)
            return k + 1;
        }
        if (isRoot())
          Alert::Warning() << "Picard did not converge at t = " << m_t << Alert::Raise;
        return k;
      }

      /// Relative L2 increment ||new - old|| / ||new|| (collective).
      template <class GF>
      static Real increment(const GF& current, const GF& previous)
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
        return delta / (scale > 0.0 ? scale : 1.0);
      }

      /// OSS projections Pi_h of the residuals. Called once per time step, at
      /// the converged state, so every Picard iteration of the next step sees
      /// Pi_h at t^n.
      void updateStabilisation()
      {
        m_piConv.project(Mult(Jacobian(m_uIt), m_uIt));
        m_piGradP.project(Grad(m_pIt));
        m_piDiv.project(Div(m_uIt));
        m_piConvPhi.project(Dot(m_uIt, Grad(m_phiIt)));
        m_piConvT.project(Dot(m_uIt, Grad(m_TIt)));
      }

      /// Interface normal and curvature at the current iterate phi^k (physics,
      /// not stabilisation: refreshed at every Picard iteration).
      void updateInterface()
      {
        const auto g = Grad(m_phiIt);
        const Real eps = m_normalEps;
        const auto normal = (1.0 / Sqrt(Dot(g, g) + eps * eps)) * g;

        m_piNormal.project(normal);
        m_piKappa.project(-1.0 * Div(m_piNormal.get()));
      }

      static MeshType makeMesh(const Context::MPI& context, const std::string& path)
      {
        const auto& comm = context.getCommunicator();
        Rodin::MPI::Sharder sharder(context);
        if (comm.rank() == 0)
        {
          Mesh<Context::Local> mesh;
          mesh.load(path, IO::FileFormat::MEDIT);
          if (mesh.getSpaceDimension() != 3 || mesh.getDimension() != 3)
            throw std::runtime_error("ForcedConvection expects a tetrahedral mesh.");
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
        connect(mesh);
        mesh.reconcile(2);
        mesh.reconcile(1);
        return mesh;
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

      static Real cellSize(const Point& p)
      {
        return std::pow(p.getPolytope().getMeasure(), 1.0 / p.getPolytope().getDimension());
      }

      Real speed(const Point& p) const
      {
        const auto u = m_uIt.getValue(p);
        return std::sqrt(Math::dot(u, u));
      }

      // ---- Mixture properties and Lee coefficients at the last iterate ----
      // Evaluated at quadrature points from phi^k, T^k (P2), never from nodal
      // averages. phi is clipped only for property evaluation.

      Real phiAt(const Point& p) const
      {
        return std::min(std::max(m_phiIt.getValue(p), Real(0)), Real(1));
      }

      Real rhoAt(const Point& p) const
      {
        const Real f = phiAt(p);
        return f * m_cfg.rhoV + (1.0 - f) * m_cfg.rhoL;
      }

      Real muAt(const Point& p) const
      {
        const Real f = phiAt(p);
        return f * m_cfg.muV + (1.0 - f) * m_cfg.muL;
      }

      /// Harmonic mean: correct normal heat flux across the interface.
      Real kAt(const Point& p) const
      {
        const Real f = phiAt(p);
        return 1.0 / (f / m_cfg.kV + (1.0 - f) / m_cfg.kL);
      }

      Real rhoCpAt(const Point& p) const
      {
        const Real f = phiAt(p);
        return f * m_cfg.rhoV * m_cfg.cpV + (1.0 - f) * m_cfg.rhoL * m_cfg.cpL;
      }

      /// Smooth evaporation switch s in [0, 1]: 1 above T_sat, 0 below.
      Real switchAt(const Point& p) const
      {
        return 0.5 * (1.0 + std::tanh((m_TIt.getValue(p) - m_cfg.Tsat) / m_cfg.dTsmooth));
      }

      /// (T^k - T_sat)/T_sat
      Real thetaAt(const Point& p) const
      {
        return (m_TIt.getValue(p) - m_cfg.Tsat) / m_cfg.Tsat;
      }

      /// beta(phi^k, T^k): mdot = beta (T - T_sat)/T_sat.
      Real betaAt(const Point& p) const
      {
        const Real f = phiAt(p), s = switchAt(p);
        return s * m_cfg.rL * (1.0 - f) * m_cfg.rhoL + (1.0 - s) * m_cfg.rV * f * m_cfg.rhoV;
      }

      /// d beta / d phi at the iterate.
      Real dbetaAt(const Point& p) const
      {
        const Real s = switchAt(p);
        return -s * m_cfg.rL * m_cfg.rhoL + (1.0 - s) * m_cfg.rV * m_cfg.rhoV;
      }

      /// tau = [(2/dt)^2 + (2|u^k|/h)^2 + (4 D/h^2)^2 + s^2]^{-1/2}
      /// with D a diffusivity and s a (nonnegative) reaction rate.
      Real tau(const Point& p, Real diffusivity, Real reaction) const
      {
        const Real h = cellSize(p);
        const Real a = 2.0 / m_cfg.dt, b = 2.0 * speed(p) / h, c = 4.0 * diffusivity / (h * h);
        return 1.0 / std::sqrt(a * a + b * b + c * c + reaction * reaction);
      }

      void setupFlow()
      {
        const size_t dim = m_mesh.getSpaceDimension();
        const Real dt = m_cfg.dt;
        const auto& L = m_cfg.labels;
        const auto n = BoundaryNormal(m_mesh);

        const auto convU = Mult(Jacobian(m_u), m_uIt);    // (grad u) u^k
        const auto streamV = Mult(Jacobian(m_v), m_uIt);  // (grad v) u^k
        const auto symU = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));
        const auto symV = 0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));
        const auto backflow = 0.5 * m_rho * Max(-Dot(m_uIt, n), 0.0);
        const RealFunction pOut = m_cfg.outletPressure;

        const auto gravity = VectorFunction(dim, [this, dim](const Point&) {
          Math::SpatialVector<Real> g(dim);
          for (size_t i = 0; i < dim; ++i)
            g(i) = 0.0;
          g(dim - 1) = -m_cfg.gravity;
          return g;
        });

        // CSF: sigma kappa^k grad phi^k, kappa^k = -div Pi_h[n^k].
        const auto capillary = m_cfg.sigma * m_piKappa.get() * Grad(m_phiIt);

        // Product of parabolas in (y, z), mean G/rho_l, smooth ramp in time.
        const auto inflow = VectorFunction(dim, [this, dim](const Point& p) {
          const Real W = m_cfg.inletWidth, d = m_cfg.depth;
          const Real s = std::min(m_t / m_cfg.ramp, 1.0);
          const Real U = 2.25 * m_cfg.massFlux / m_cfg.rhoL * 0.5 * (1.0 - std::cos(M_PI * s));
          const Real eta = 2.0 * p.y() / W, zeta = 2.0 * p.z() / d - 1.0;
          Math::SpatialVector<Real> u(dim);
          for (size_t i = 0; i < dim; ++i)
            u(i) = 0.0;
          u(0) = std::max(U * (1.0 - eta * eta) * (1.0 - zeta * zeta), 0.0);
          return u;
        });

        // Galerkin with rho^k, mu^k; the convective form is not
        // skew-symmetrised because div u = mdot (1/rho_v - 1/rho_l) != 0.
        // Split OSS: tau_M rho ((grad u) u^k - Pi[(grad u^k) u^k], (grad v) u^k)
        //          + tau_M/rho (grad p - Pi[grad p^k], grad q)
        //          + tau_C (div u - Pi[div u^k], div v).
        m_flow = (1.0 / dt) * Integral(m_rho * m_u, m_v) - (1.0 / dt) * Integral(m_rho * m_uOld, m_v)
               + Integral(m_rho * convU, m_v)
               + 2.0 * Integral(m_mu * symU, symV)
               - Integral(m_p, Div(m_v))
               + Integral(Div(m_u), m_q) - Integral(m_divSource, m_q)
               - Integral(m_rho * gravity, m_v)
               - Integral(capillary, m_v)

               + Integral(m_tauM * m_rho * convU, streamV)
               - Integral(m_tauM * m_rho * m_piConv.get(), streamV)

               + Integral((m_tauM / m_rho) * Grad(m_p), Grad(m_q))
               - Integral((m_tauM / m_rho) * m_piGradP.get(), Grad(m_q))

               + Integral(m_tauC * Div(m_u), Div(m_v))
               - Integral(m_tauC * m_piDiv.get(), Div(m_v))

               + BoundaryIntegral(pOut * Dot(m_v, n)).over(L.outlet)
               + BoundaryIntegral(backflow * Dot(m_u, m_v)).over(L.outlet)

               + DirichletBC(m_u, Zero(dim)).on(FlatSet<Attribute>{ L.bottom, L.top, L.sides })
               + DirichletBC(m_u, inflow).on(L.inlet);
      }

      void setupPhase()
      {
        const Real dt = m_cfg.dt, rhoV = m_cfg.rhoV, hLV = m_cfg.hLV;
        const auto& L = m_cfg.labels;
        const RealFunction q = m_cfg.heatFlux;
        const RealFunction Tin = m_cfg.inletTemperature;
        const RealFunction liquid = 0.0;
        const auto streamA = Dot(m_uIt, Grad(m_a));
        const auto streamW = Dot(m_uIt, Grad(m_w));
        const Real gammaEps = m_gamma * m_eps;
        const auto compression = m_gamma * (m_liquidFraction * m_piNormal.get());  // (1 - phi^k) n^k

        // Linearised Lee source about (phi^k, T^k):
        //   mdot ~= cT T + cPhi phi + c0,
        //   cT = beta/T_sat >= 0,  cPhi = beta' theta <= 0,
        //   c0 = -beta - beta' theta phi^k.
        // phi-equation: conservative transport div(phi u^k), source mdot/rho_v.
        // T-equation:   rho cp (T - T^n)/dt + rho cp u^k.grad T - div(k grad T)
        //               = -mdot h_lv.
        // Allen-Cahn (weak): gamma eps (grad phi, grad a) - gamma ((1 - phi^k) n^k phi, grad a).
        // Split OSS on the convective operators only.
        m_phase = (1.0 / dt) * Integral(m_phi, m_a) - (1.0 / dt) * Integral(m_phiOld, m_a)
                + Integral(Dot(m_uIt, Grad(m_phi)), m_a)
                + Integral(Div(m_uIt) * m_phi, m_a)
                - (1.0 / rhoV) * Integral(m_cT * m_T, m_a)
                - (1.0 / rhoV) * Integral(m_cPhi * m_phi, m_a)
                - (1.0 / rhoV) * Integral(m_c0, m_a)
                + gammaEps * Integral(Grad(m_phi), Grad(m_a))
                - Integral(compression * m_phi, Grad(m_a))

                + Integral(m_tauPhi * Dot(m_uIt, Grad(m_phi)), streamA)
                - Integral(m_tauPhi * m_piConvPhi.get(), streamA)

                + (1.0 / dt) * Integral(m_rhoCp * m_T, m_w) - (1.0 / dt) * Integral(m_rhoCp * m_TOld, m_w)
                + Integral(m_rhoCp * Dot(m_uIt, Grad(m_T)), m_w)
                + Integral(m_k * Grad(m_T), Grad(m_w))
                - BoundaryIntegral(q * m_w).over(L.bottom)
                + hLV * Integral(m_cT * m_T, m_w)
                + hLV * Integral(m_cPhi * m_phi, m_w)
                + hLV * Integral(m_c0, m_w)

                + Integral(m_tauT * m_rhoCp * Dot(m_uIt, Grad(m_T)), streamW)
                - Integral(m_tauT * m_rhoCp * m_piConvT.get(), streamW)

                + DirichletBC(m_T, Tin).on(L.inlet)
                + DirichletBC(m_phi, liquid).on(L.inlet);
      }

      template <class ProblemType>
      void solve(ProblemType& problem, Solver::KSP& ksp)
      {
        problem.assemble();
        problem.solve(ksp);
        ::KSPConvergedReason reason;
        PetscErrorCode ierr = KSPGetConvergedReason(ksp.getHandle(), &reason);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
        if (reason < 0)
          throw std::runtime_error("KSP diverged at t = " + std::to_string(m_t));
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

      /// Outlet bulk temperature, inlet-outlet mean pressure drop, vapour
      /// volume and bounds of phi.
      void report(int step, int iterations)
      {
        const auto& L = m_cfg.labels;
        const auto n = BoundaryNormal(m_mesh);
        const auto& u = m_u.getSolution();
        const auto& T = m_T.getSolution();
        const auto& p = m_p.getSolution();
        const auto& phi = m_phi.getSolution();

        const Real flow = boundaryIntegral(Dot(u, n), L.outlet);
        const Real enthalpy = boundaryIntegral(T * Dot(u, n), L.outlet);
        const Real dp = boundaryIntegral(p, L.inlet) / m_inletArea
                      - boundaryIntegral(p, L.outlet) / m_outletArea;
        const Real bulk = (flow > 0.0) ? enthalpy / flow : m_cfg.inletTemperature;
        const Real vapour = volumeIntegral(phi);
        const Real umax = std::max(std::abs(u.max()), std::abs(u.min()));
        const Real Tmax = T.max();
        const Real phiMin = phi.min(), phiMax = phi.max();

        if (!isRoot())
          return;
        m_csv << m_t << ',' << iterations << ',' << umax << ',' << Tmax << ',' << phiMin << ','
              << phiMax << ',' << vapour << ',' << bulk << ',' << dp << '\n';
        m_csv.flush();
        if (step % m_cfg.printEvery == 0)
          Alert::Info() << "t=" << m_t << " s  it=" << iterations << "  max|u|=" << umax
                        << " m/s  maxT=" << Tmax << " K  phi in [" << phiMin << ", " << phiMax
                        << "]  Vv=" << vapour << " m^3  Tb,out=" << bulk << " K  dp=" << dp
                        << " Pa" << Alert::Raise;
      }

      bool isRoot() const { return m_mesh.getContext().getCommunicator().rank() == 0; }

      Config m_cfg;
      MeshType m_mesh;
      IO::XDMF m_xdmf;
      VectorFES m_vh;
      ScalarFES m_sh;
      ScalarFES2 m_sh2;

      // Flow unknowns (P1/P1) and phase unknowns (P2/P2).
      VectorTrial m_u;
      ScalarTrial m_p;
      VectorTest m_v;
      ScalarTest m_q;
      ScalarTrial2 m_phi;
      ScalarTrial2 m_T;
      ScalarTest2 m_a;
      ScalarTest2 m_w;

      // Time level n and Picard iterate k.
      VectorGF m_uOld, m_uIt;
      ScalarGF m_pIt;
      ScalarGF2 m_phiOld, m_TOld, m_phiIt, m_TIt;
      ScalarGF m_phiOut, m_TOut;  ///< P1 copies for XDMF

      // Mixture properties, Lee coefficients and stabilisation parameters.
      Coefficient m_rho, m_mu, m_k, m_rhoCp, m_liquidFraction;
      Coefficient m_cT, m_cPhi, m_c0, m_divSource;
      Coefficient m_tauM, m_tauT, m_tauPhi, m_tauC;

      // Orthogonal projections Pi_h at the iterate.
      OrthogonalProjection<VectorFES> m_piConv, m_piGradP, m_piNormal;
      OrthogonalProjection<ScalarFES> m_piDiv, m_piKappa;
      OrthogonalProjection<ScalarFES2> m_piConvPhi, m_piConvT;

      FlowProblem m_flow;
      Solver::KSP m_flowKSP;
      PhaseProblem m_phase;
      Solver::KSP m_phaseKSP;

      ScalarTest m_z;
      ScalarGF m_one;
      LinearForm<ScalarFES, ::Vec> m_flux;

      Real m_t = 0.0;
      Real m_eps = 0.0, m_gamma = 0.0, m_normalEps = 0.0;
      Real m_bottomArea = 0.0, m_inletArea = 0.0, m_outletArea = 0.0;
      std::ofstream m_csv;
  };
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  const auto setDefault = [](const char* name, const char* value) {
    PetscBool set = PETSC_FALSE;
    PetscErrorCode ierr = PetscOptionsHasName(PETSC_NULLPTR, PETSC_NULLPTR, name, &set);
    if (ierr == PETSC_SUCCESS && !set)
      ierr = PetscOptionsSetValue(PETSC_NULLPTR, name, value);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
  };
  // Flow: direct (MUMPS). Phase block: GMRES + block Jacobi/ILU.
  // Projections: CG + Jacobi on the mass matrices.
  setDefault("-fc_flow_ksp_type", "preonly");
  setDefault("-fc_flow_pc_type", "lu");
  setDefault("-fc_flow_pc_factor_mat_solver_type", "mumps");
  // This MUMPS build faults in its distributed RHS scatter: keep the
  // right-hand side and solution centralized (as in the Heart examples).
  for (const char* prefix : { "-fc_flow_", "-fc_phase_" })
  {
    setDefault((std::string(prefix) + "mat_mumps_icntl_20").c_str(), "0");
    setDefault((std::string(prefix) + "mat_mumps_icntl_21").c_str(), "0");
  }
  setDefault("-fc_phase_ksp_type", "gmres");
  setDefault("-fc_phase_pc_type", "bjacobi");
  setDefault("-fc_phase_ksp_rtol", "1e-8");
  setDefault("-fc_phase_ksp_gmres_restart", "100");
  for (const char* prefix : { "-fc_proj_conv_", "-fc_proj_gradp_", "-fc_proj_normal_",
                              "-fc_proj_div_", "-fc_proj_kappa_", "-fc_proj_phi_", "-fc_proj_T_" })
  {
    setDefault((std::string(prefix) + "ksp_type").c_str(), "cg");
    setDefault((std::string(prefix) + "pc_type").c_str(), "jacobi");
    setDefault((std::string(prefix) + "ksp_rtol").c_str(), "1e-10");
  }

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world(PETSC_COMM_WORLD, boost::mpi::comm_attach);
  Rodin::Context::MPI context(env, world);

  int status = 0;
  try
  {
    using Rodin::Examples::Heat::ForcedConvection;
    ForcedConvection::Config cfg;

    char buffer[512];
    PetscBool got = PETSC_FALSE;
    PetscOptionsGetString(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_mesh", buffer, sizeof(buffer), &got);
    if (got)
      cfg.meshPath = buffer;
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_dt", &cfg.dt, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_tend", &cfg.tEnd, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_G", &cfg.massFlux, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_q", &cfg.heatFlux, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_tin", &cfg.inletTemperature, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_g", &cfg.gravity, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_tol", &cfg.tolerance, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_eps", &cfg.epsilon, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_gamma", &cfg.gamma, PETSC_NULLPTR);
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_neps", &cfg.normalEps, PETSC_NULLPTR);
    Rodin::Real r = cfg.rL;
    PetscBool gotR = PETSC_FALSE;
    PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_r", &r, &gotR);
    if (gotR)
      cfg.rL = cfg.rV = r;
    PetscInt every = cfg.outputEvery, maxit = cfg.maxIterations;
    PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_output_every", &every, PETSC_NULLPTR);
    PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_maxit", &maxit, PETSC_NULLPTR);
    PetscInt print = cfg.printEvery;
    PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_print_every", &print, PETSC_NULLPTR);
    cfg.printEvery = std::max(static_cast<int>(print), 1);
    PetscBool monitor = PETSC_FALSE;
    PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_picard_monitor", &monitor, PETSC_NULLPTR);
    cfg.picardMonitor = (monitor == PETSC_TRUE);
    cfg.outputEvery = static_cast<int>(every);
    cfg.maxIterations = static_cast<int>(maxit);

    ForcedConvection simulation(context, cfg);
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
