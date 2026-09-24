// ForcedConvectionNewton.cpp
//
// Same problem as ForcedConvection, solved monolithically in (u, p, T) with
// Newton (PETSc SNES) at every backward-Euler step, so the flow and the energy
// equation are iterated together and u^{n+1} convects T^{n+1}.
//
// Stabilisation: split (term-by-term) orthogonal subgrid scales (OSGS),
//   rho tau_1 ( (grad u) u - Pi[(grad u) u], (grad v) u )   convection
//   tau_1/rho ( grad p     - Pi[grad p],     grad q )       pressure
//   tau_2     ( div u      - Pi[div u],      div v )        divergence
//   rho cp tau_T ( u.grad T - Pi[u.grad T],  u.grad w )     energy
// with Pi the L2 projection onto P1 (OrthogonalProjection.h) of the fields at
// t^n (lagged once per step, as tau), and
//   tau_1 = [4 nu/h^2 + 2|u^n|/h]^{-1},  tau_2 = rho h^2/(4 tau_1),
//   tau_T = [4 alpha/h^2 + 2|u^n|/h]^{-1}.
// With Pi and tau lagged the Jacobian is exact, including the derivative of
// the test operators (grad v) u and u.grad w. Pi evaluated at the iterate
// instead would leave -tau Pi' out of the Jacobian, and on smooth modes
// Pi ~ I, so Newton would degrade to a Picard iteration with contraction
// close to one.
//
// Boundary labels and data as in ForcedConvection.cpp.
//
// Run from the build directory:
//   mpirun -n 4 ./examples/Heat/ForcedConvectionNewton -fc_dt 2.5e-4 -fc_tend 0.2 \
//     -snes_monitor -snes_converged_reason
//
// Options: -fc_mesh, -fc_dt, -fc_tend, -fc_G, -fc_q, -fc_tin, -fc_output_every;
// Newton via -snes_* (defaults rtol 1e-8, atol 1e-10, 20 iterations), its
// linear solver via -ksp_*/-pc_* (MUMPS), projections under -fcn_proj_ (CG).
#include <cassert>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

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

  class ForcedConvectionNewton
  {
    public:
      using MeshType = Mesh<Context::MPI>;
      using VectorFES = H1<1, Math::SpatialVector<Real>, MeshType>;
      using ScalarFES = H1<1, Real, MeshType>;
      using VectorGF = PETSc::Variational::GridFunction<VectorFES>;
      using ScalarGF = PETSc::Variational::GridFunction<ScalarFES>;
      using VectorTrial = PETSc::Variational::TrialFunction<VectorGF, VectorFES>;
      using ScalarTrial = PETSc::Variational::TrialFunction<ScalarGF, ScalarFES>;
      using VectorTest = PETSc::Variational::TestFunction<VectorFES>;
      using ScalarTest = PETSc::Variational::TestFunction<ScalarFES>;
      using System = PETSc::Math::LinearSystem;
      using CoupledProblem = Problem<System, VectorTrial, ScalarTrial, ScalarTrial,
                                     VectorTest, ScalarTest, ScalarTest>;
      using CellFES = P0<Real, MeshType>;
      using CellGF = PETSc::Variational::GridFunction<CellFES>;
      using Clock = std::chrono::steady_clock;

      struct Labels
      {
        Attribute inlet = 1, outlet = 2, bottom = 3, top = 4, sides = 5;
      };

      struct Config
      {
        std::string meshPath = "../resources/examples/Heat/DivergingMicrochannel3D.mesh";
        std::string output = "ForcedConvectionNewton";
        Labels labels;

        Real rho = 996.5;    ///< kg/m^3
        Real mu = 2.82e-4;   ///< Pa s
        Real k = 0.606;      ///< W/(m K)
        Real cp = 4179.0;    ///< J/(kg K)

        Real massFlux = 240.0;          ///< G (kg/m^2 s); mean inlet velocity G/rho
        Real inletWidth = 2.0e-4;       ///< m, inlet centred on y = 0
        Real depth = 4.0e-4;            ///< m, z in [0, depth]
        Real inletTemperature = 300.0;  ///< K
        Real heatFlux = 1.5e5;          ///< W/m^2 on the bottom
        Real outletPressure = 0.0;      ///< Pa (gauge)

        Real dt = 2.5e-4;
        Real tEnd = 0.2;
        Real ramp = 5.0e-3;  ///< s, smooth start of the inflow
        int outputEvery = 20;
      };

      ForcedConvectionNewton(const Context::MPI& context, const Config& cfg)
        : m_cfg(cfg),
          m_mesh(makeMesh(context, m_cfg.meshPath)),
          m_xdmf(context.getCommunicator(), m_cfg.output),
          m_vh(std::integral_constant<size_t, 1>{}, m_mesh, m_mesh.getSpaceDimension()),
          m_sh(std::integral_constant<size_t, 1>{}, m_mesh),
          m_ch(m_mesh),
          m_u(m_vh), m_p(m_sh), m_T(m_sh), m_v(m_vh), m_q(m_sh), m_w(m_sh),
          m_uOld(m_vh), m_TOld(m_sh),
          m_piConv(m_vh, "fcn_proj_"), m_piGradP(m_vh, "fcn_proj_"),
          m_piDiv(m_sh, "fcn_proj_"), m_piAdv(m_sh, "fcn_proj_"),
          m_tau1(m_ch), m_tau2(m_ch), m_tauT(m_ch),
          m_problem(m_u, m_p, m_T, m_v, m_q, m_w), m_ksp(m_problem), m_snes(m_ksp),
          m_z(m_sh), m_one(m_sh), m_flux(m_z)
      {
        const size_t dim = m_mesh.getSpaceDimension();
        Math::SpatialVector<Real> zero(dim);
        for (size_t i = 0; i < dim; ++i)
          zero(i) = 0.0;
        for (auto* g : { &m_uOld, &m_u.getSolution(), &m_piConv.get(), &m_piGradP.get() })
          *g = zero;
        for (auto* g : { &m_p.getSolution(), &m_piDiv.get(), &m_piAdv.get() })
          *g = Real(0);
        m_TOld = m_cfg.inletTemperature;
        m_T.getSolution() = m_cfg.inletTemperature;
        m_one = Real(1);

        const auto& L = m_cfg.labels;
        m_bottomArea = boundaryIntegral(m_one, L.bottom);
        m_inletArea = boundaryIntegral(m_one, L.inlet);
        m_outletArea = boundaryIntegral(m_one, L.outlet);

        setupProblem();

        // The SNES iterate x = (u, p, T) lives in the trial solutions.
        const size_t nu = m_vh.getSize(), np = m_sh.getSize();
        m_snes.setTolerances(1e-10, 1e-8, 1e-12, 20, 10000)
          .setStateUpdate([this, nu, np](const PETSc::Math::Vector& x) {
            m_u.getSolution().setData(x, 0);
            m_p.getSolution().setData(x, nu);
            m_T.getSolution().setData(x, nu + np);
          });

        m_xdmf.setMesh(m_mesh);
        m_xdmf.add("velocity", m_u.getSolution());
        m_xdmf.add("pressure", m_p.getSolution());
        m_xdmf.add("temperature", m_T.getSolution());

        const auto cells = m_mesh.getCellCount(), vertices = m_mesh.getVertexCount();  // collective
        if (isRoot())
        {
          m_csv.open(m_cfg.output + ".csv");
          m_csv << "t,newtonIterations,maxU,maxT,bulkTout,dp,energyBalance\n";
          Alert::Info() << "[mesh] cells=" << cells << " vertices=" << vertices
                        << "  bottom area=" << m_bottomArea << " m^2" << Alert::Raise;
        }
      }

      int run()
      {
        const int steps = static_cast<int>(m_cfg.tEnd / m_cfg.dt + 0.5);
        m_problem.assemble();  // seeds the SNES iterate from the trial solutions
        for (int n = 1; n <= steps; ++n)
        {
          m_t = n * m_cfg.dt;

          const auto t0 = Clock::now();
          updateStabilization();  // tau and Pi at (u^n, p^n, T^n)
          m_snes.solve();
          const Real elapsed = std::chrono::duration<Real>(Clock::now() - t0).count();
          if (!m_snes.converged())
            throw std::runtime_error("Newton did not converge at t = " + std::to_string(m_t));

          m_uOld.setData(m_u.getSolution().getData());
          m_TOld.setData(m_T.getSolution().getData());

          report(n, steps, elapsed);
          if (n % m_cfg.outputEvery == 0 || n == steps)
            m_xdmf.write(m_t).flush();
        }
        m_xdmf.close();
        return 0;
      }

    private:
      static MeshType makeMesh(const Context::MPI& context, const std::string& path)
      {
        const auto& comm = context.getCommunicator();
        Rodin::MPI::Sharder sharder(context);
        if (comm.rank() == 0)
        {
          Mesh<Context::Local> mesh;
          mesh.load(path, IO::FileFormat::MEDIT);
          if (mesh.getSpaceDimension() != 3 || mesh.getDimension() != 3)
            throw std::runtime_error("ForcedConvectionNewton expects a tetrahedral mesh.");
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

      /// tau = [4 D/h^2 + 2|u^n|/h]^{-1}, D = nu or alpha; lagged in time, so
      /// the Newton Jacobian needs no derivative of it.
      Real tau(const Point& p, Real diffusivity) const
      {
        const auto u = m_uOld.getValue(p);
        const Real h = cellSize(p);
        return 1.0 / (4.0 * diffusivity / (h * h) + 2.0 * std::sqrt(Math::dot(u, u)) / h);
      }

      /// Cellwise tau (P0) and the projections, both at t^n. P0 matters: a
      /// coefficient of known order lets the residual be integrated with the
      /// same rule as the Jacobian, which keeps Newton quadratic.
      void updateStabilization()
      {
        const Real nu = m_cfg.mu / m_cfg.rho, alpha = m_cfg.k / (m_cfg.rho * m_cfg.cp);
        m_tau1.project(RealFunction([this, nu](const Point& p) { return tau(p, nu); }));
        m_tau2.project(RealFunction([this, nu](const Point& p) {
          const Real h = cellSize(p);
          return m_cfg.rho * h * h / (4.0 * tau(p, nu)); }));
        m_tauT.project(RealFunction([this, alpha](const Point& p) { return tau(p, alpha); }));

        const auto& U = m_u.getSolution();
        m_piConv.project(Mult(Jacobian(U), U));
        m_piGradP.project(Grad(m_p.getSolution()));
        m_piDiv.project(Div(U));
        m_piAdv.project(Dot(U, Grad(m_T.getSolution())));
      }

      void setupProblem()
      {
        const size_t dim = m_mesh.getSpaceDimension();
        const Real rho = m_cfg.rho, dt = m_cfg.dt, C = m_cfg.rho * m_cfg.cp;
        const auto& L = m_cfg.labels;
        const auto n = BoundaryNormal(m_mesh);

        // State (current Newton iterate) and corrections (trial functions).
        const auto& U = m_u.getSolution();
        const auto& P = m_p.getSolution();
        const auto& Th = m_T.getSolution();
        const auto& du = m_u;
        const auto& dp = m_p;
        const auto& dT = m_T;

        const auto convU = Mult(Jacobian(U), U);                                // (grad U) U
        const auto dConvU = Mult(Jacobian(du), U) + Mult(Jacobian(U), du);      // its derivative
        const auto streamV = Mult(Jacobian(m_v), U);                            // (grad v) U
        const auto streamW = Dot(U, Grad(m_w));                                 // U.grad w
        const auto advT = Dot(U, Grad(Th));                                     // U.grad T
        const auto dAdvT = Dot(U, Grad(dT));                                    // its derivative,
        const auto dAdvU = Dot(du, Grad(Th));                                   // split by field
        const auto backflow = 0.5 * rho * Max(-Dot(m_uOld, n), 0.0);

        const RealFunction pOut = m_cfg.outletPressure;
        const RealFunction q = m_cfg.heatFlux;

        // Dirichlet data for the corrections: g - state.
        const auto inflowCorrection = VectorFunction(dim, [this, dim](const Point& p) {
          const Real W = m_cfg.inletWidth, d = m_cfg.depth;
          const Real s = std::min(m_t / m_cfg.ramp, 1.0);
          const Real Umax = 2.25 * m_cfg.massFlux / m_cfg.rho * 0.5 * (1.0 - std::cos(M_PI * s));
          const Real eta = 2.0 * p.y() / W, zeta = 2.0 * p.z() / d - 1.0;
          Math::SpatialVector<Real> g = m_u.getSolution().getValue(p);
          for (size_t i = 0; i < dim; ++i)
            g(i) = -g(i);
          g(0) += std::max(Umax * (1.0 - eta * eta) * (1.0 - zeta * zeta), 0.0);
          return g;
        });
        const RealFunction temperatureCorrection = [this](const Point& p) {
          return m_cfg.inletTemperature - m_T.getSolution().getValue(p);
        };

        m_problem =
            // ---- Jacobian: momentum --------------------------------------
              Integral(Dot((rho / dt) * du + 0.5 * rho * (Div(du) * U) + 0.5 * rho * (Div(U) * du), m_v))
            + Integral(dConvU, rho * m_v + rho * (m_tau1 * streamV))
            + rho * Integral(du, Mult(Transpose(Jacobian(m_v)), m_tau1 * (convU - m_piConv.get())))
            + m_cfg.mu * Integral(Jacobian(du), Jacobian(m_v))
            + Integral(m_tau2 * Div(du), Div(m_v))
            - Integral(dp, Div(m_v))
            + BoundaryIntegral(backflow * Dot(du, m_v)).over(L.outlet)

            // ---- Jacobian: continuity ------------------------------------
            + Integral(Div(du), m_q)
            + (1.0 / rho) * Integral(m_tau1 * Grad(dp), Grad(m_q))

            // ---- Jacobian: energy ----------------------------------------
            + Integral((C / dt) * dT + C * dAdvT, m_w)
            + C * Integral(m_tauT * dAdvT, streamW)
            + m_cfg.k * Integral(Grad(dT), Grad(m_w))
            + Integral(dAdvU, C * m_w + C * (m_tauT * streamW))
            + C * Integral(m_tauT * (advT - m_piAdv.get()) * du, Grad(m_w))

            // ---- Residual: momentum --------------------------------------
            + Integral(Dot((rho / dt) * (U - m_uOld) + rho * convU + 0.5 * rho * (Div(U) * U), m_v))
            + rho * Integral(m_tau1 * (convU - m_piConv.get()), streamV)
            + m_cfg.mu * Integral(Jacobian(U), Jacobian(m_v))
            + Integral(m_tau2 * (Div(U) - m_piDiv.get()) - P, Div(m_v))
            + BoundaryIntegral(Dot(backflow * U + pOut * n, m_v)).over(L.outlet)

            // ---- Residual: continuity ------------------------------------
            + Integral(Div(U), m_q)
            + (1.0 / rho) * Integral(m_tau1 * (Grad(P) - m_piGradP.get()), Grad(m_q))

            // ---- Residual: energy ----------------------------------------
            + Integral((C / dt) * (Th - m_TOld) + C * advT, m_w)
            + C * Integral(m_tauT * (advT - m_piAdv.get()), streamW)
            + m_cfg.k * Integral(Grad(Th), Grad(m_w))
            - BoundaryIntegral(q * m_w).over(L.bottom)

            // ---- Dirichlet conditions on the corrections: g - state ------
            + DirichletBC(du, -U).on(FlatSet<Attribute>{ L.bottom, L.top, L.sides })
            + DirichletBC(du, inflowCorrection).on(L.inlet)
            + DirichletBC(dT, temperatureCorrection).on(L.inlet);
      }

      template <class Expression>
      Real boundaryIntegral(const Expression& f, Attribute a)
      {
        m_flux = BoundaryIntegral(f, m_z).over(a);
        m_flux.assemble();
        return m_flux(m_one);
      }

      /// Outlet bulk temperature, inlet-outlet mean pressure drop and
      /// rho cp int_out (T - T_in) u.n / (q'' A_bottom), which tends to 1.
      void report(int step, int steps, Real elapsed)
      {
        const auto& L = m_cfg.labels;
        const auto n = BoundaryNormal(m_mesh);
        const auto& u = m_u.getSolution();
        const auto& T = m_T.getSolution();
        const auto& p = m_p.getSolution();

        const Real flow = boundaryIntegral(Dot(u, n), L.outlet);
        const Real enthalpy = boundaryIntegral(T * Dot(u, n), L.outlet);
        const Real dp = boundaryIntegral(p, L.inlet) / m_inletArea
                      - boundaryIntegral(p, L.outlet) / m_outletArea;
        const Real bulk = (flow > 0.0) ? enthalpy / flow : m_cfg.inletTemperature;
        const Real balance = m_cfg.rho * m_cfg.cp *
          (enthalpy - m_cfg.inletTemperature * flow) / (m_cfg.heatFlux * m_bottomArea);
        const Real umax = std::max(std::abs(u.max()), std::abs(u.min()));
        const Real Tmax = T.max();
        const int its = static_cast<int>(m_snes.getIterationNumber());

        if (!isRoot())
          return;
        m_csv << m_t << ',' << its << ',' << umax << ',' << Tmax << ',' << bulk << ','
              << dp << ',' << balance << '\n';
        m_csv.flush();
        Alert::Info() << "step " << step << "/" << steps << "  t=" << m_t
                      << " s  newton=" << its << "  max|u|=" << umax << " m/s  maxT=" << Tmax
                      << " K  Tb,out=" << bulk << " K  dp=" << dp
                      << " Pa  balance=" << balance << "  [" << elapsed << " s]"
                      << Alert::Raise;
        std::cout.flush();
      }

      bool isRoot() const { return m_mesh.getContext().getCommunicator().rank() == 0; }

      Config m_cfg;
      MeshType m_mesh;
      IO::XDMF m_xdmf;
      VectorFES m_vh;
      ScalarFES m_sh;
      CellFES m_ch;

      // Trial functions are Newton corrections; their solutions hold the state.
      VectorTrial m_u;
      ScalarTrial m_p;
      ScalarTrial m_T;
      VectorTest m_v;
      ScalarTest m_q;
      ScalarTest m_w;
      VectorGF m_uOld;
      ScalarGF m_TOld;

      OrthogonalProjection<VectorFES> m_piConv, m_piGradP;
      OrthogonalProjection<ScalarFES> m_piDiv, m_piAdv;

      CellGF m_tau1, m_tau2, m_tauT;

      CoupledProblem m_problem;
      Solver::KSP m_ksp;
      Solver::SNES m_snes;

      ScalarTest m_z;
      ScalarGF m_one;
      LinearForm<ScalarFES, ::Vec> m_flux;

      Real m_t = 0.0;
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
  // Newton step: direct (MUMPS). Projections: CG + Jacobi on the mass matrix.
  // SNES configures its KSP from the unprefixed options, as in CoronaryArtery.
  setDefault("-ksp_type", "preonly");
  setDefault("-pc_type", "lu");
  setDefault("-pc_factor_mat_solver_type", "mumps");
  setDefault("-mat_mumps_icntl_20", "0");  // centralized RHS and solution,
  setDefault("-mat_mumps_icntl_21", "0");  // as in the Heart examples
  setDefault("-fcn_proj_ksp_type", "cg");
  setDefault("-fcn_proj_pc_type", "jacobi");
  setDefault("-fcn_proj_ksp_rtol", "1e-10");

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world(PETSC_COMM_WORLD, boost::mpi::comm_attach);
  Rodin::Context::MPI context(env, world);

  int status = 0;
  try
  {
    using Rodin::Examples::Heat::ForcedConvectionNewton;
    ForcedConvectionNewton::Config cfg;

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
    PetscInt every = cfg.outputEvery;
    PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-fc_output_every", &every, PETSC_NULLPTR);
    cfg.outputEvery = static_cast<int>(every);

    ForcedConvectionNewton simulation(context, cfg);
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
