#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <type_traits>

#include <boost/mpi/collectives.hpp>

#include <petscvec.h>

#include <Rodin/Alert.h>
#include <Rodin/Configure.h>
#include <Rodin/Solid.h>
#include <Rodin/Solver.h>

#ifdef RODIN_USE_SCOTCH
#include <Rodin/Scotch/MeshPartitioner.h>
#endif

#include <Rodin/Alert/Exception.h>
#include <Rodin/Alert/NewLine.h>
#include <Rodin/Math/RootFinding/NewtonRaphson.h>
#include <Rodin/Math/RungeKutta/RK4.h>
#include <Rodin/Variational/ForwardDecls.h>

#include "CoronaryArteryAlerts.h"
#include "CoupledLV0DCoronary3D.h"

namespace Rodin::Examples::Heart
{
  using namespace Rodin;
  using namespace Rodin::Math;
  using namespace Rodin::Solver;
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  namespace
  {
    constexpr int RootRank = 0;

    const char *flowModeName(CoupledLV0DCoronary3D::FlowMode mode)
    {
      switch (mode)
      {
        case CoupledLV0DCoronary3D::FlowMode::Newton:
          return "newton";
        case CoupledLV0DCoronary3D::FlowMode::Oseen:
          return "oseen";
      }
      return "unknown";
    }

    class VecSnapshot
    {
      public:
        explicit VecSnapshot(::Vec source)
        {
          PetscErrorCode ierr = VecDuplicate(source, &m_data);
          assert(ierr == PETSC_SUCCESS);
          ierr = VecCopy(source, m_data);
          assert(ierr == PETSC_SUCCESS);
          (void)ierr;
        }

        VecSnapshot(const VecSnapshot &) = delete;
        VecSnapshot &operator=(const VecSnapshot &) = delete;

        ~VecSnapshot()
        {
          if (m_data)
          {
            PetscErrorCode ierr = VecDestroy(&m_data);
            assert(ierr == PETSC_SUCCESS);
            (void)ierr;
          }
        }

      void restore(::Vec target) const
      {
          PetscErrorCode ierr = VecCopy(m_data, target);
          assert(ierr == PETSC_SUCCESS);
          (void)ierr;
      }

    private:
      ::Vec m_data = PETSC_NULLPTR;
    };
  } // namespace

  CoupledLV0DCoronary3D::CoupledLV0DCoronary3D(const Context::MPI &context)
      : CoupledLV0DCoronary3D(context, Config{})
  {}

  CoupledLV0DCoronary3D::CoupledLV0DCoronary3D(
      const Context::MPI &context, const Config &cfg)
    : m_cfg(cfg),
      m_input(makeInput(m_cfg)),
      m_model(m_input),
      m_mesh(makeMesh(context, m_cfg)),
      m_xdmf(context.getCommunicator(), m_cfg.xdmfBasename),

      /*
       * Function spaces
       */
      m_uh(std::integral_constant<size_t, 2>{},
           m_mesh,
           m_mesh.getSpaceDimension()),

      m_ph(std::integral_constant<size_t, 1>{}, m_mesh),

      m_dh(std::integral_constant<size_t, 1>{},
           m_mesh,
           m_mesh.getSpaceDimension()),

      m_uph(std::integral_constant<size_t, 2>{},
            m_mesh,
            m_mesh.getSpaceDimension()),

      m_tauh(m_mesh),

      /*
       * Trial functions
       */
      m_u(m_uh),
      m_p(m_ph),
      m_d(m_dh),

      m_mu(m_ph),

      /*
       * Test functions
       */
      m_v(m_uh),
      m_q(m_ph),
      m_w(m_dh),

      m_r(m_ph),

      /*
       * Old-time / auxiliary grid functions
       */
      m_uOld(m_uh),
      m_pOld(m_ph),
      m_dOld(m_dh),
      m_dOldOld(m_dh),

      m_one(m_ph),
      m_qFlux(m_ph),

      /*
       * VMS / projection fields
       */
      m_sub(m_uph),
      m_subOld(m_uph),
      m_up(m_uph),
      m_vp(m_uph),

      m_tau(m_tauh),
      m_t(m_tauh),
      m_tauOld(m_tauh),

      /*
       * Forms and solvers
       */
      m_flux(m_qFlux),

      m_flow(m_u, m_p, m_d, m_v, m_q, m_w),
      m_flowKSP(m_flow),
      m_flowSolver(m_flowKSP),

      m_viscosityProjection(m_mu, m_r),
      m_viscosityProjectionKSP(m_viscosityProjection),

      m_l2ConvU(m_up, m_vp),
      m_l2ConvUSolver(m_l2ConvU),

      m_subProjection(m_sub, m_vp),
      m_subProjectionSolver(m_subProjection),

      m_tauProjection(m_tau, m_t),
      m_tauProjectionSolver(m_tauProjection)
  {
    auto cellCount = m_mesh.getCellCount();
    auto vertexCount = m_mesh.getVertexCount();
    auto spaceDim = m_mesh.getSpaceDimension();
    auto velocityDOFs = m_uh.getSize();
    auto pressureDOFs = m_ph.getSize();

    if (isRoot())
    {
      Alert::Info() << "---- Mesh ----" << Alert::NewLine
                    << "Number of elements in mesh: " << cellCount
                    << Alert::NewLine
                    << "Number of vertices in mesh: " << vertexCount
                    << Alert::NewLine << "Mesh space dimension: " << spaceDim
                    << Alert::Raise;

      Alert::Info() << "---- Function spaces ----" << Alert::NewLine
                    << "Velocity space has " << velocityDOFs << " DOFs."
                    << Alert::NewLine << "Pressure space has " << pressureDOFs
                    << " DOFs." << Alert::Raise;
    }
  }

  CoupledLV0DCoronary3D::~CoupledLV0DCoronary3D() = default;

  CoupledLV0DCoronary3D::MeshType CoupledLV0DCoronary3D::makeMesh(
      const Context::MPI &context, const Config &cfg)
  {
    const auto &comm = context.getCommunicator();

    Rodin::MPI::Sharder sharder(context);
    if (comm.rank() == RootRank)
    {
      Geometry::Mesh<Context::Local> mesh;
      mesh.load(cfg.meshPath, IO::FileFormat::MEDIT);

      if (mesh.getSpaceDimension() != 3)
        throw std::runtime_error("Expected a 3D coronary mesh.");

      Alert::Info() << "Computing connectivity for " << cfg.meshPath << " ..."
                    << Alert::Raise;

      const size_t D = mesh.getDimension();
      mesh.getConnectivity().compute(D, D);
      mesh.getConnectivity().compute(D, 0);
      mesh.getConnectivity().compute(D, D - 1);
      mesh.getConnectivity().compute(D - 1, D);
      mesh.getConnectivity().compute(D - 1, 0);
      mesh.getConnectivity().compute(D - 1, 1);
      mesh.getConnectivity().compute(1, 0);
      mesh.save("CoronaryArtery.medit.mesh", IO::FileFormat::MEDIT);

      Alert::Info() << "Partitioning coronary mesh over " << comm.size()
                    << " MPI ranks ..." << Alert::Raise;

#ifdef RODIN_USE_SCOTCH
      Alert::Info() << "Using SCOTCH mesh partitioner." << Alert::Raise;
      Scotch::Partitioner partitioner(mesh);
#else
      Alert::Info() << "Using balanced compact mesh partitioner." << Alert::Raise;
      Geometry::BalancedCompactPartitioner partitioner(mesh);
#endif
      partitioner.partition(static_cast<size_t>(comm.size()));
      sharder.shard(partitioner);
      sharder.scatter(RootRank);
    }

    MeshType mesh = sharder.gather(RootRank);
    mesh.scale(cfg.meshScale);

    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);
    mesh.getConnectivity().compute(D, 0);
    mesh.getConnectivity().compute(D, D - 1);
    mesh.getConnectivity().compute(D - 1, D);
    mesh.getConnectivity().compute(D - 1, 0);
    mesh.getConnectivity().compute(D - 1, 1);
    mesh.getConnectivity().compute(1, 0);
    mesh.reconcile(2);
    mesh.reconcile(1);

    mesh.save(cfg.xdmfBasename + "_partitioned" +
                  std::to_string(mesh.getContext().getCommunicator().rank()) +
                  ".mesh",
              IO::FileFormat::MEDIT);

    return mesh;
  }

  CoupledLV0DCoronary3D::Model::Input
  CoupledLV0DCoronary3D::makeInput(const Config &cfg)
  {
    Model::Input input;

    input.rho = cfg.lv.rho;
    input.R0 = cfg.lv.R0;
    input.d0 = cfg.lv.d0;

    input.Es = cfg.lv.Es;
    input.mu = cfg.lv.mu;
    input.eta = cfg.lv.eta;
    input.alpha = cfg.lv.alpha;
    input.alphaR = cfg.lv.alphaR;
    input.k0 = cfg.lv.k0;
    input.sigma0 = cfg.lv.sigma0;

    input.Rp = cfg.lv.Rp;
    input.Cp = cfg.lv.Cp;
    input.Rd = cfg.lv.Rd;
    input.Cd = cfg.lv.Cd;

    input.proximalRadius = cfg.lv.proximalRadius;
    input.proximalLength = cfg.lv.proximalLength;
    input.distalRadius = cfg.lv.distalRadius;
    input.distalLength = cfg.lv.distalLength;

    input.mu_0 = cfg.lv.mu_0;
    input.mu_Inf = cfg.lv.mu_Inf;
    input.lambda = cfg.lv.lambda;
    input.n = cfg.lv.n;
    input.yasuda = cfg.lv.yasuda;

    input.Kat = cfg.lv.Kat;
    input.Kp = cfg.lv.Kp;
    input.Kar = cfg.lv.Kar;

    input.cavityCapacity = cfg.lv.cavityCapacity;

    input.localTolerance = cfg.lv.localTolerance;
    input.localMaxIterations = cfg.lv.localMaxIterations;
    input.localDamping = cfg.lv.localDamping;
    input.absRegularization = cfg.lv.absRegularization;

    input.initFibDef = cfg.lv.initFibDef;
    input.initActiveStiffness = cfg.lv.initActiveStiffness;
    input.initActiveStress = cfg.lv.initActiveStress;

    input.pSv = [p = cfg.lv.systemicVenousPressure](Real) { return p; };
    input.pAt = [p = cfg.atrialPressure](Real t)
    {
      return atrial_pressure(p, t);
    };
    input.u = [a = cfg.activation](Real t) { return periodic_activation(a, t); };
    input.m0 =
      [ low = cfg.lv.relaxationM0Low, high = cfg.lv.relaxationM0High,
        lowEc = cfg.lv.relaxationM0LowEc,
        highEc = cfg.lv.relaxationM0HighEc](Real ec)
    {
      if (highEc <= lowEc)
        return high;
      if (ec <= lowEc)
        return low;
      if (ec >= highEc)
        return high;
      const Real s = (ec - lowEc) / (highEc - lowEc);
      return (1.0 - s) * low + s * high;
    };

    input.dm0 = [low = cfg.lv.relaxationM0Low, high = cfg.lv.relaxationM0High,
                 lowEc = cfg.lv.relaxationM0LowEc,
                 highEc = cfg.lv.relaxationM0HighEc](Real ec)
    {
      if (highEc <= lowEc || ec <= lowEc || ec >= highEc)
        return 0.0;
      return (high - low) / (highEc - lowEc);
    };

    {
      using PassiveEnergy = std::decay_t<decltype(input.passiveEnergy)>;

      typename PassiveEnergy::Parameters hp;
      hp.mu1 = cfg.lv.passiveMu1;
      hp.mu2 = cfg.lv.passiveMu2;
      hp.C0 = cfg.lv.passiveC0;
      hp.C1 = cfg.lv.passiveC1;
      hp.C2 = cfg.lv.passiveC2;
      hp.C3 = cfg.lv.passiveC3;

      input.passiveEnergy = PassiveEnergy(hp);
    }

    return input;
  }

  void CoupledLV0DCoronary3D::updateRCR(
      const Model &model, RCR &bc, Real Q, Real dt)
  {
    const Real a = bc.C / dt;

    const auto &s = model.getState();
    const auto &h = model.getHistory();

    const Real dPim = (s.pv - h.nm1.pv) / dt;
    const Real Qim = bc.C * 0.3 * dPim;

    bc.pc = (a * bc.pc + Q + Qim + s.pv / bc.Rd) / (a + 1.0 / bc.Rd);

    bc.qd = (bc.pc - s.pv) / bc.Rd;
    bc.pout = bc.pc + bc.Rp * Q;
  }

  void CoupledLV0DCoronary3D::updateRCRNonNew(
      const Config &cfg, const Attribute &tag,
      const Model &model, RCR &bc, Real Q, Real dt)
  {
    const auto &s = model.getState();
    const auto &h = model.getHistory();

    const Real cap = bc.C / dt;
    const Real pcOld = bc.pc;

    const auto &law = cfg.outletFlowLaw;
    const auto &cy = cfg.viscosity;

    const Real radiusP = law.geometricParam.at(tag).Rp;
    const Real lengthP = law.geometricParam.at(tag).Lp;

    const Real radiusD = law.geometricParam.at(tag).Rd;
    const Real lengthD = law.geometricParam.at(tag).Ld;

    auto flowLaw = [&](Real dp, Real L, Real radius) -> std::pair<Real, Real>
    {
      const Real mu0 = cy.mu0;
      const Real muInf = cy.muInf;
      const Real lambda = cy.lambda;
      const Real n = cy.n;
      const Real yasuda = cy.yasuda;
      const Real delta = mu0 - muInf;

      const Real sgn = (dp >= 0.0) ? 1.0 : -1.0;
      const Real adp = std::abs(dp);

      const Real R0 =
          8.0 * mu0 * L / (std::numbers::pi_v<Real> * std::pow(radius, 4.0));

      if (adp < law.pressureDropTolerance)
        return {dp / R0, 1.0 / R0};

      const Real tauW = radius * adp / (2.0 * L);

      auto mu = [&](Real g) -> Real
      {
        return muInf + delta * std::pow(1.0 + std::pow(lambda * g, yasuda),
                                        (n - 1.0) / yasuda);
      };

      auto dmu = [&](Real g) -> Real
      {
        const Real base = 1.0 + std::pow(lambda * g, yasuda);

        return delta * (n - 1.0) * std::pow(base, (n - 1.0 - yasuda) / yasuda) *
               std::pow(lambda, yasuda) * std::pow(g, yasuda - 1.0);
      };

      auto tauMinusTauW = [&](Real g) -> std::pair<Real, Real>
      {
        const Real m = mu(g);
        const Real dm = dmu(g);
        return {g * m - tauW, m + g * dm};
      };

      Math::RootFinding::NewtonRaphson<Real> rootFinder(
          law.shearAbsoluteTolerance, law.shearRelativeTolerance,
          law.shearStepTolerance, law.shearMaxIterations);

      Real gHi = std::max<Real>(tauW / muInf, law.minShearRate);

      for (int k = 0;
           k < law.maxBracketIterations && tauMinusTauW(gHi).first < 0.0; ++k)
        gHi *= 2.0;

      if (tauMinusTauW(gHi).first < 0.0) {
        std::cerr << "Warning: failed to bracket wall shear rate. "
                  << "Using Poiseuille fallback.\n";
        return {dp / R0, 1.0 / R0};
      }

      const auto gammaRoot =
          rootFinder.solve(tauMinusTauW, 0.5 * gHi, law.shearStepTolerance, gHi);

      if (!gammaRoot)
      {
        std::cerr << "Warning: failed to solve wall shear rate. "
                  << "Using Poiseuille fallback.\n";
        return {dp / R0, 1.0 / R0};
      }

      const Real gammaW = *gammaRoot;

      auto integrand = [&](Real g) -> Real
      {
        if (g <= 0.0)
          return 0.0;

        const Real m = mu(g);
        const Real dm = dmu(g);
        const Real dtau = m + g * dm;

        return std::pow(g, 3.0) * m * m * dtau;
      };

      Math::RungeKutta::RK4 integrator;

      const int steps = law.integralSteps;
      const Real h = gammaW / static_cast<Real>(steps);

      Real I = 0.0;

      auto rhs = [&](Real g, Real y) -> Real
      {
        (void) y;
        return integrand(g);
      };

      for (int i = 0; i < steps; ++i)
      {
        const Real g = static_cast<Real>(i) * h;
        integrator.step(I, g, h, I, rhs);
      }

      if (I <= 0.0 || !std::isfinite(I))
      {
        std::cerr << "Warning: invalid WRMS integral. "
                  << "Using Poiseuille fallback.\n";
        return {dp / R0, 1.0 / R0};
      }

      const Real qAbs = std::numbers::pi_v<Real> * std::pow(radius, 3.0) * I /
                        std::pow(tauW, 3.0);

      const Real dqAbs =
          (std::numbers::pi_v<Real> * std::pow(radius, 3.0) * gammaW -
           3.0 * qAbs) /
          adp;

      if (!std::isfinite(qAbs) || !std::isfinite(dqAbs) || dqAbs <= 0.0)
      {
        std::cerr << "Warning: invalid WRMS flow derivative. "
                  << "Using Poiseuille fallback.\n";
        return {dp / R0, 1.0 / R0};
      }

      return {sgn * qAbs, dqAbs};
    };

    auto solvePressureDropForFlow = [&](Real targetQ, Real L, Real radius,
                                        Real guess) -> Real
    {
      if (std::abs(targetQ) < law.zeroFlowTolerance)
        return 0.0;

      const Real sgn = (targetQ >= 0.0) ? 1.0 : -1.0;
      const Real qAbs = std::abs(targetQ);

      auto F = [&](Real x) -> std::pair<Real, Real>
      {
        const auto [q, dq] = flowLaw(sgn * x, L, radius);
        return {sgn * q - qAbs, dq};
      };

      Real hi = std::max<Real>(std::abs(guess), law.pressureDropBracketMin);

      for (int k = 0; k < law.maxBracketIterations && F(hi).first < 0.0; ++k)
        hi *= 2.0;

      if (F(hi).first < 0.0)
      {
        std::cerr << "Warning: failed to bracket pressure drop for targetQ = "
                  << targetQ << ". Returning last upper bound.\n";
        return sgn * hi;
      }

      Math::RootFinding::NewtonRaphson<Real> solver(
          law.flowAbsoluteTolerance, law.flowRelativeTolerance,
          law.flowStepTolerance, law.flowMaxIterations);

      const auto root = solver.solve(F, std::min(std::abs(guess), hi), 0.0, hi);

      if (!root)
      {
        std::cerr << "Warning: failed to invert flow law for targetQ = "
                  << targetQ << ". Returning bracket upper bound.\n";
        return sgn * hi;
      }

      return sgn * (*root);
    };

    const Real alpha = 0.7;
    const Real dPim = alpha * (s.pv - h.nm1.pv) / dt;
    const Real Qim = bc.C * dPim;

    auto distalResidual = [&](Real pc) -> std::pair<Real, Real>
    {
      const Real pim = alpha * s.pv;
      const Real x = std::max(pc - pim, Real(0.0));
      const auto [qd, dqd] = flowLaw(x, lengthD, radiusD);

      const Real f = cap * (pc - pcOld) - Qim + qd - Q;
      const Real df = cap + (pc > pim ? dqd : Real(0.0));
      return {f, df};
    };

    Math::RootFinding::NewtonRaphson<Real> solver(
        law.flowAbsoluteTolerance, law.flowRelativeTolerance,
        law.flowStepTolerance, law.flowMaxIterations);

    Real span = std::max<Real>(std::abs(Q) / cap + law.distalPressureBracketPad,
                               law.distalPressureBracketPad);

    Real lo = std::min(pcOld, s.pv) - span;
    Real hi = std::max(pcOld, s.pv) + span;

    for (int k = 0; k < law.maxBracketIterations && distalResidual(lo).first * distalResidual(hi).first > 0.0; ++k)
    {
      span *= 2.0;
      lo = std::min(pcOld, s.pv) - span;
      hi = std::max(pcOld, s.pv) + span;
    }

    if (distalResidual(lo).first * distalResidual(hi).first > 0.0)
    {
      std::cerr << "Warning: failed to bracket distal capacitor pressure. "
                << "Keeping previous pc.\n";
      bc.pc = pcOld;
    }
    else
    {
      const auto pcNew = solver.solve(distalResidual, pcOld, lo, hi);

      if (!pcNew)
      {
        std::cerr << "Warning: failed to solve distal capacitor equation. "
                  << "Keeping previous pc.\n";
        bc.pc = pcOld;
      }
      else
      {
        bc.pc = *pcNew;
      }
    }

    const Real pim_f = alpha * s.pv;
    const auto [qd, dqd_f] = flowLaw(std::max(bc.pc - pim_f, Real(0.0)), lengthD, radiusD);
    (void)dqd_f;
    bc.qd = qd;

    const Real oldGuess = bc.pout - bc.pc;
    const Real dpP = solvePressureDropForFlow(Q, lengthP, radiusP, oldGuess);

    bc.pout = bc.pc + dpP;
  }

  CoupledLV0DCoronary3D::Real
  CoupledLV0DCoronary3D::periodic_activation(const Activation &cfg, Real t)
  {
    const Real T = cfg.period;
    const Real tau = t - T * std::floor(t / T);

    if (tau < cfg.tRampStart)
      return 0.0;
    if (tau < cfg.tRampEnd)
      return cfg.positiveValue *
             ((tau - cfg.tRampStart) / (cfg.tRampEnd - cfg.tRampStart));
    if (tau < cfg.tPlateauEnd)
      return cfg.positiveValue;
    if (tau < cfg.tRelaxEnd)
      return cfg.positiveValue +
             (cfg.negativeValue - cfg.positiveValue) *
                 ((tau - cfg.tPlateauEnd) / (cfg.tRelaxEnd - cfg.tPlateauEnd));
    if (tau < cfg.tNegativeEnd)
      return cfg.negativeValue;

    return 0.0;
  }

  CoupledLV0DCoronary3D::Real
  CoupledLV0DCoronary3D::atrial_pressure(const AtrialPressure &cfg, Real t)
  {
    const Real T = cfg.period;
    const Real tau = t - T * std::floor(t / T);

    Real alpha = 0.0;
    Real value = cfg.minValue;

    if (tau < cfg.t1)
    {
      alpha = -(tau - cfg.t1) / cfg.t1;
      value = alpha * cfg.minValue + (1.0 - alpha) * cfg.maxValue;
    }
    else if (tau < cfg.t2)
    {
      value = cfg.maxValue;
    }
    else if (tau < cfg.t3)
    {
      alpha = -(tau - cfg.t3) / (cfg.t3 - cfg.t2);
      value = alpha * cfg.maxValue + (1.0 - alpha) * cfg.minValue;
    }
    else if (tau < cfg.t4)
    {
      alpha = -(tau - cfg.t4) / (cfg.t4 - cfg.t3);
      value = alpha * cfg.minValue + (1.0 - alpha) * cfg.secondThreshold;
    }
    else if (tau < cfg.t5)
    {
      value = cfg.secondThreshold;
    }
    else if (tau < cfg.t6)
    {
      alpha = -(tau - cfg.t6) / (cfg.t6 - cfg.t5);
      value = alpha * cfg.secondThreshold + (1.0 - alpha) * cfg.minValue;
    }

    return value;
  }

  CoupledLV0DCoronary3D &CoupledLV0DCoronary3D::initialize()
  {
    setupModel();
    setupMeshAndSpaces();
    setupDiagnostics();
    printInitialState();

    m_initialized = true;
    return *this;
  }

  void CoupledLV0DCoronary3D::setupModel()
  {
    m_model.setMaxIterations(m_cfg.lv.maxIterations)
        .setAbsoluteTolerance(m_cfg.lv.absoluteTolerance)
        .setRelativeTolerance(m_cfg.lv.relativeTolerance)
        .setStepTolerance(m_cfg.lv.stepTolerance)
        .setDampingFactor(m_cfg.lv.dampingFactor);

    Model::State s0;
    s0.t = 0.0;
    s0.y = m_cfg.lv.initialY;
    s0.v = m_cfg.lv.initialV;
    s0.pv = m_input.pAt(0.0) + m_cfg.lv.initialPvOffset;
    s0.par = m_cfg.lv.initialPar;
    s0.pd = m_cfg.lv.initialPd;
    s0.ec = m_input.initFibDef;
    s0.gamma = std::sqrt(std::max<Real>(m_input.initActiveStiffness, 0.0));
    s0.beta = (s0.gamma > 0.0) ? (m_input.initActiveStress / s0.gamma) : 0.0;
    s0.kc = s0.gamma * s0.gamma;
    s0.tauc = s0.gamma * s0.beta;
    s0.w = m_input.m0(s0.ec);

    m_model.initialize(s0);
  }

  void CoupledLV0DCoronary3D::setupMeshAndSpaces()
  {
    if (isRoot())
      Alert::Info() << "Setting up " << m_cfg.xdmfBasename << ".xdmf ..."
                    << Alert::Raise;

    // The ALE mesh moves every step, so the geometry must be re-exported at
    // every snapshot.
    m_xdmf.setMesh(m_mesh, IO::XDMF::MeshPolicy::Transient);

    // Capture the undeformed configuration used as the ALE reference.
    saveReferenceVertices();

    int rank = m_mesh.getContext().getCommunicator().rank();
    m_mesh.save("CoronaryArtery." + std::to_string(rank) + ".medit.mesh", IO::FileFormat::MEDIT);

    m_u.setName("FluidVelocity");
    m_p.setName("FluidPressure");
    m_d.setName("Displacement");
    m_mu.setName("BloodViscosity");
    m_up.setName("ProjectedConvection");
    m_sub.setName("SubgridVelocity");
    m_tau.setName("StabilizationTau");

    m_uOld = Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};
    m_pOld = 0.0;
    m_dOld = Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};
    m_dOldOld = Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};

    m_u.getSolution() = Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};
    m_p.getSolution() = 0.0;
    m_d.getSolution() = Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};

    m_one = 1.0;
    m_mu.getSolution() = 0.0;
    m_subOld = Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};
    m_tauOld = 0.0;

    m_xdmf.add("FluidVelocity", m_u.getSolution());
    m_xdmf.add("FluidPressure", m_p.getSolution());
    m_xdmf.add("Displacement", m_d.getSolution());
    m_xdmf.add("BloodViscosity", m_mu.getSolution());
    m_xdmf.add("SubgridVelocity", m_sub.getSolution());
    m_xdmf.add("StabilizationTau", m_tau.getSolution());

    m_wk.clear();
    for (const Attribute tag : m_cfg.outlets)
      m_wk.emplace(tag, m_cfg.defaultRCR);

    m_flowSolver.setTolerances(1e-10, 1e-8, 1e-10, 50, 10000)
        .setStateUpdate(
            [this](const PETSc::Math::Vector &x)
            {
              const size_t uOffset = 0;
              const size_t pOffset = uOffset + m_uh.getSize();
              const size_t dOffset = pOffset + m_ph.getSize();

              m_u.getSolution().setData(x, uOffset);
              m_p.getSolution().setData(x, pOffset);
              m_d.getSolution().setData(x, dOffset);
            });
  }

  void CoupledLV0DCoronary3D::saveReferenceVertices()
  {
    m_referenceVertices.resize(m_mesh.getVertexCount());

    for (auto it = m_mesh.getVertex(); it; ++it)
      m_referenceVertices[it->getIndex()] =
          m_mesh.getVertexCoordinates(it->getIndex());
  }

  void CoupledLV0DCoronary3D::moveMeshToDisplacement(
      const DisplacementGridFunctionType &d)
  {
    assert(m_mesh.getVertexCount() == m_referenceVertices.size());

    const size_t dim = m_mesh.getSpaceDimension();

    // Absolute ALE positioning: x = x_reference + d. Reapplying the same
    // displacement is idempotent, so this is safe to call repeatedly.
    for (auto it = m_mesh.getVertex(); it; ++it)
    {
      const Index vertex = it->getIndex();
      auto x = m_referenceVertices[vertex];

      for (Index c = 0; c < static_cast<Index>(dim); ++c)
        x(c) += d[m_dh.getGlobalIndex({0, vertex}, c)];

      m_mesh.setVertexCoordinates(vertex, x);
    }

    m_mesh.flush();
  }

  void CoupledLV0DCoronary3D::setupDiagnostics()
  {
    if (!isRoot())
      return;

    m_csv.open(m_cfg.csvPath);

    if (!m_csv)
      throw std::runtime_error("Failed to open coupled CSV file: " +
                               m_cfg.csvPath);

    writeCSVHeader();
  }

  void CoupledLV0DCoronary3D::printInitialState() const
  {
    if (!isRoot())
      return;

    const auto &s = m_model.getState();

    Alert::Info() << "Initial 0D state:" << Alert::NewLine << "y     = " << s.y
                  << Alert::NewLine << "v     = " << s.v << Alert::NewLine
                  << "pv    = " << s.pv << Alert::NewLine << "par   = " << s.par
                  << Alert::NewLine << "pd    = " << s.pd << Alert::NewLine
                  << "ec    = " << s.ec << Alert::NewLine << "gamma = " << s.gamma
                  << Alert::NewLine << "beta  = " << s.beta << Alert::NewLine
                  << "w     = " << s.w << Alert::NewLine << "kc    = " << s.kc
                  << Alert::NewLine << "tauc  = " << s.tauc << Alert::NewLine
                  << "3D flow mode: " << flowModeName(m_cfg.flowMode)
                  << Alert::Raise;
  }

  bool CoupledLV0DCoronary3D::isRoot() const
  {
    return m_mesh.getContext().getCommunicator().rank() == RootRank;
  }

  bool CoupledLV0DCoronary3D::advance0D()
  {
    if (isRoot())
      ZeroDInfo() << "Advancing LV model ..." << Alert::Raise;

    const auto rep = m_model.step(m_cfg.dt);

    if (isRoot())
    {
      ZeroDInfo() << "Newton: " << (rep.converged ? "converged" : "NOT converged")
                  << "  iter = " << rep.iterations
                  << "  |F| = " << rep.finalResidual
                  << "  |dx| = " << rep.finalStepNorm << Alert::Raise;
    }

    if (isRoot() && rep.converged)
    {
      const auto &s = m_model.getState();
      ZeroDInfo() << "pv  = " << s.pv << " Pa"
                  << "  par = " << s.par << " Pa"
                  << "  pd  = " << s.pd << " Pa" << Alert::Raise;
    }

    return rep.converged;
  }

  void CoupledLV0DCoronary3D::solveStatic()
  {
    if (isRoot())
      ThreeDInfo() << "Solving static initialization ..." << Alert::Raise;

    const auto normal = BoundaryNormal(m_mesh);
    const Attribute outlet0 = m_cfg.outlets[0];
    const Attribute outlet1 = m_cfg.outlets[1];
    const Attribute outlet2 = m_cfg.outlets[2];
    const Attribute outlet3 = m_cfg.outlets[3];
    const Attribute outlet4 = m_cfg.outlets[4];
    const Attribute outlet5 = m_cfg.outlets[5];

    const auto &s = m_model.getState();
    const Real pin = s.par;

    const Real mu = m_cfg.viscosity.mu0;

    const auto symDU = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));

    const auto symV = 0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));

    const auto &uLag = m_uOld;
    const auto oseenConvection = Mult(Jacobian(m_u), uLag); // (uLag·∇)du
    const auto divLag = Div(uLag);
    const auto oseenTemam = divLag * Dot(m_u, m_v); // Temam stabilization

    // Backflow damping (keeps BCs stable even in steady solve)
    const auto outletBeta = Max(-Dot(m_uOld, normal), 0.0);
    const auto outletBackflowDamping =
        0.5 * m_cfg.outletBackflowStabilization * m_cfg.rho * outletBeta;

    m_flow =
        2.0 * Integral(mu * symDU, symV)

        - Integral(m_p, Div(m_v)) + Integral(Div(m_u), m_q) +
        m_cfg.eps * Integral(m_p, m_q)

        + BoundaryIntegral(pin * Dot(m_v, normal)).over(m_cfg.inlet)

        +
        BoundaryIntegral(m_wk.at(outlet0).pout * Dot(m_v, normal)).over(outlet0) +
        BoundaryIntegral(m_wk.at(outlet1).pout * Dot(m_v, normal)).over(outlet1) +
        BoundaryIntegral(m_wk.at(outlet2).pout * Dot(m_v, normal)).over(outlet2) +
        BoundaryIntegral(m_wk.at(outlet3).pout * Dot(m_v, normal)).over(outlet3) +
        BoundaryIntegral(m_wk.at(outlet4).pout * Dot(m_v, normal)).over(outlet4) +
        BoundaryIntegral(m_wk.at(outlet5).pout * Dot(m_v, normal)).over(outlet5)

        + BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet0) +
        BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet1) +
        BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet2) +
        BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet3) +
        BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet4) +
        BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet5)

        + DirichletBC(m_u, Zero(m_mesh.getSpaceDimension())).on(m_cfg.wall);

    m_flow.assemble();

    if (!m_flowFieldSplitsSet)
    {
      m_flow.setFieldSplits();
      m_flowFieldSplitsSet = true;
    }

    m_flow.solve(m_flowKSP);

    ::KSPConvergedReason reason;
    PetscErrorCode ierr = KSPGetConvergedReason(m_flowKSP.getHandle(), &reason);
    assert(ierr == PETSC_SUCCESS);

    PetscInt iterations = 0;
    ierr = KSPGetIterationNumber(m_flowKSP.getHandle(), &iterations);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;

    if (isRoot()) {
      KSPInfo() << "Static solve: "
                << (reason > 0 ? "Converged" : "Did NOT converge")
                << "  iterations = " << iterations << Alert::Raise;
    }
  }

  bool CoupledLV0DCoronary3D::solve3D()
  {
    const auto setup3DStart = CoronaryClock::now();

    const auto &s = m_model.getState();

    // Inlet/outlet pressures: use the 0D par / RCR pout unless a constant
    // diagnostic override is configured. A pressure ramp eases the whole
    // luminal load on smoothly to avoid the impulsive startup shock.
    const Real rampFactor =
        (m_cfg.pressureRampTime > 0.0)
            ? std::min(s.t / m_cfg.pressureRampTime, Real(1.0))
            : Real(1.0);
    const Real pin = rampFactor * m_cfg.inletPressureOverride.value_or(s.par);
    auto outletPressure = [&](const Attribute tag) -> Real
    {
      return rampFactor * m_cfg.outletPressureOverride.value_or(m_wk.at(tag).pout);
    };

    const auto normal = BoundaryNormal(m_mesh);
    const Attribute outlet0 = m_cfg.outlets[0];
    const Attribute outlet1 = m_cfg.outlets[1];
    const Attribute outlet2 = m_cfg.outlets[2];
    const Attribute outlet3 = m_cfg.outlets[3];
    const Attribute outlet4 = m_cfg.outlets[4];
    const Attribute outlet5 = m_cfg.outlets[5];

    // ---- SOLID PARAMETERS (linear elasticity + weak FSI coupling)
    //
    // Hooke / linear isotropic elasticity:
    //   sigma_s(d) = lambdaS (div d) I + 2 muS eps(d).
    const Real solidRho = m_cfg.solidRho;
    const Real solidNu = m_cfg.solidNu;
    const Real solidE = m_cfg.solidE;
    const Real lambdaS =
        solidE * solidNu / ((1.0 + solidNu) * (1.0 - 2.0 * solidNu));
    const Real muS = solidE / (2.0 * (1.0 + solidNu));
    const Real rayleighAlpha = m_cfg.solidRayleighAlpha;
    const Real rayleighBeta = m_cfg.solidRayleighBeta;

    const Real inactivePenalty = m_cfg.inactivePenalty;

    // Weak FSI kinematic penalty:  gammaTilde (u - (d - dOld)/dt) on Gamma_fsi.
    const Real gammaTilde = m_cfg.fsiPenalty / m_cfg.dt;
    // ------------------------------------------

    const auto &uState = m_u.getSolution();
    const auto &pState = m_p.getSolution();
    const auto &uLag = m_uOld;

    const auto gradDU = Jacobian(m_u);
    const auto gradU = Jacobian(uState);
    const auto gradLag = Jacobian(uLag);

    const auto newtonConvection = Mult(gradDU, uState) + Mult(gradU, m_u);

    const auto stateConvection = Mult(gradU, uState);

    /*
     * ALE transport velocity for the Oseen fluid block.
     *
     * The fluid is advected on a moving mesh. The mesh velocity is the lagged
     * rate of the displacement field on the fluid nodes,
     *
     *   w_mesh = (d^n - d^{n-1}) / dt,
     *
     * and the convective transport uses the relative velocity
     *
     *   c = u^n - w_mesh.
     *
     * With c lagged, the convection (c . grad) u^{n+1} stays Oseen-linear.
     */
    const auto meshVelocity = (1.0 / m_cfg.dt) * (m_dOld - m_dOldOld);
    const auto aleTransport = uLag - meshVelocity;
    const auto oseenConvectionJacobian = Mult(gradDU, aleTransport);

    const auto divDU = Div(m_u);
    const auto divU = Div(uState);

    /*
     * Div does not deduce on a Sum expression, so write div(c) algebraically:
     *   div(c) = div(u^n) - (1/dt) (div d^n - div d^{n-1}).
     */
    const auto divTransport = Div(uLag) -
                              (1.0 / m_cfg.dt) * Div(m_dOld) +
                              (1.0 / m_cfg.dt) * Div(m_dOldOld);

    const auto temamJacobian1 = Dot(divDU * uState, m_v);

    const auto temamJacobian2 = divU * Dot(m_u, m_v);

    const auto temamResidual = divU * Dot(uState, m_v);

    const auto oseenTemamJacobian = divTransport * Dot(m_u, m_v);

    // ALE transport for the Newton fluid block: c = uState - w_mesh, with the
    // lagged mesh velocity (mirroring the Oseen branch). Convection is then
    // linearized fully with respect to the velocity correction m_u.
    const auto transportState = uState - meshVelocity;
    const auto newtonConvectionALE =
        Mult(gradDU, transportState) + Mult(gradU, m_u);
    const auto stateConvectionALE = Mult(gradU, transportState);
    const auto divTransportState = divU -
                                   (1.0 / m_cfg.dt) * Div(m_dOld) +
                                   (1.0 / m_cfg.dt) * Div(m_dOldOld);
    const auto temamJacobian2ALE = divTransportState * Dot(m_u, m_v);
    const auto temamResidualALE = divTransportState * Dot(uState, m_v);

    const auto outletBeta = Max(-Dot(m_uOld, normal), 0.0);
    const auto inletBeta = Max(Dot(m_uOld, normal), 0.0);
    const auto outletBackflowDamping =
        0.5 * m_cfg.outletBackflowStabilization * m_cfg.rho * outletBeta;
    const auto inletBackflowDamping =
        0.5 * m_cfg.inletBackflowStabilization * m_cfg.rho * inletBeta;
    const auto symDU = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));

    const auto symV = 0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));

    // Symmetric gradients of the solid displacement trial / test functions,
    // and of the previous displacement (for stiffness-proportional damping).
    const auto symD = 0.5 * (Jacobian(m_d) + Transpose(Jacobian(m_d)));
    const auto symW = 0.5 * (Jacobian(m_w) + Transpose(Jacobian(m_w)));
    const auto symDOld =
        0.5 * (Jacobian(m_dOld) + Transpose(Jacobian(m_dOld)));

    const auto symU = 0.5 * (Jacobian(uState) + Transpose(Jacobian(uState)));

    const auto symLag = 0.5 * (gradLag + Transpose(gradLag));

    const auto duNormal = Dot(m_u, normal) * normal;

    const auto duTangential = m_u - duNormal;

    const auto uStateNormal = Dot(uState, normal) * normal;

    const auto uStateTangential = uState - uStateNormal;

    const auto &cy = m_cfg.viscosity;
    const Real gammaReg = cy.gammaRegularization;
    const Real mu0 = cy.mu0;
    const Real muInf = cy.muInf;
    const Real lambda = cy.lambda;
    const Real nCY = cy.n;
    const Real yasuda = cy.yasuda;
    const Real deltaMu = mu0 - muInf;

    const auto gamma = Sqrt(gammaReg * gammaReg + 2.0 * Dot(symU, symU));

    const auto gammaLag = Sqrt(gammaReg * gammaReg + 2.0 * Dot(symLag, symLag));

    const auto carreauBase = 1.0 + Pow(lambda * gamma, yasuda);

    const auto carreauBaseLag = 1.0 + Pow(lambda * gammaLag, yasuda);

    const auto mu = muInf + deltaMu * Pow(carreauBase, (nCY - 1.0) / yasuda);

    const auto muLag =
        muInf + deltaMu * Pow(carreauBaseLag, (nCY - 1.0) / yasuda);

    const auto dgamma = 2.0 * Dot(symU, symDU) / gamma;

    const auto dmu = deltaMu * (nCY - 1.0) *
                     Pow(carreauBase, (nCY - 1.0 - yasuda) / yasuda) *
                     std::pow(lambda, yasuda) * Pow(gamma, yasuda - 1.0) * dgamma;

    /*
     * Convective residual target:
     *
     *   a^n = (grad u^n) u^n.
     */
    const auto convectionTarget = Mult(Jacobian(m_uOld), m_uOld);

    /*
     * L2 projection of the convective term:
     *
     *   int proj_conv_h · v = int u^n · grad u^n  · v.
     */

    // Fluid-only L2 projection; zero in the solid (solid mass term, no RHS).
    m_l2ConvU = Integral(m_up, m_vp).over(m_cfg.fluidVolume) -
                Integral(convectionTarget, m_vp).over(m_cfg.fluidVolume) +
                Integral(m_up, m_vp).over(m_cfg.solidVolume);
    m_l2ConvU.assemble();
    m_l2ConvUSolver.solve();

    RealFunction tau = [=, this](const Point &p) -> Real {
      const auto uOld = m_uOld.getValue(p);
      const Real mu = m_mu.getSolution().getValue(p);
      const Real hK = std::pow(p.getPolytope().getMeasure(),
                               1.0 / p.getPolytope().getDimension());
      const Real order = 2;
      const Real speed = std::sqrt(dot(uOld, uOld));

      const Real Tau =
          1. / (4.0 * std::pow(order, 4.) * mu / (m_cfg.rho * std::pow(hK, 2.0)) +
                2.0 * order * speed / hK);

      return 1. / (m_cfg.rho / m_cfg.dt + m_cfg.rho / Tau);
    };

    // Fluid-only projection; zero in the solid (solid mass term, no RHS).
    m_tauProjection = Integral(m_tau, m_t).over(m_cfg.fluidVolume) -
                      Integral(tau, m_t).over(m_cfg.fluidVolume) +
                      Integral(m_tau, m_t).over(m_cfg.solidVolume);
    m_tauProjection.assemble();
    m_tauProjectionSolver.solve();

    auto subUpdate = VectorFunction(
        m_mesh.getSpaceDimension(),
        [=, this](const Point &p) -> Math::SpatialVector<Real>
        {
          const auto conv = convectionTarget.getValue(p);
          const auto proj = m_up.getSolution().getValue(p);
          const auto old = m_subOld.getValue(p);
          const auto tau = m_tau.getSolution().getValue(p);

          Math::SpatialVector<Real> out(m_mesh.getSpaceDimension());

          for (size_t c = 0; c < out.size(); ++c) {
            out[c] =
                tau * m_cfg.rho * (1. / m_cfg.dt * old[c] - (conv[c] - proj[c]));
          }

          return out;
        });

    /*
     * L2 projection of the dynamic subscale into m_sub:
     *
     *   int sub_h · v = int subUpdate · v.
     */

    // Fluid-only projection; zero in the solid (solid mass term, no RHS).
    m_subProjection = Integral(m_sub, m_vp).over(m_cfg.fluidVolume) -
                      Integral(subUpdate, m_vp).over(m_cfg.fluidVolume) +
                      Integral(m_sub, m_vp).over(m_cfg.solidVolume);
    m_subProjection.assemble();
    m_subProjectionSolver.solve();

    if (m_cfg.flowMode == FlowMode::Newton)
    {
      // Solid state and nonlinear NeoHookean operators, evaluated at the
      // current displacement iterate dState.
      const auto &dState = m_d.getSolution();
      const auto symDState =
          0.5 * (Jacobian(dState) + Transpose(Jacobian(dState)));

      Rodin::Solid::NeoHookean solidLaw(lambdaS, muS);
      Rodin::Solid::MaterialTangent solidTangent(solidLaw, m_d, m_w, dState);
      Rodin::Solid::InternalForce solidInternal(solidLaw, m_w, dState);

      m_flow =
          /*
           * =========================
           * Fluid Newton tangent (ALE Navier-Stokes over the fluid)
           * =========================
           * Unknowns are Newton corrections m_u, m_p, m_d.
           */
          (m_cfg.rho / m_cfg.dt) * Integral(m_u, m_v).over(m_cfg.fluidVolume)

          // ALE convection: rho ((c . grad) du + (du . grad) uState, v).
          + m_cfg.rho *
                Integral(Dot(newtonConvectionALE, m_v)).over(m_cfg.fluidVolume)

          // Temam/skew correction tangent.
          + 0.5 * m_cfg.rho * Integral(temamJacobian1).over(m_cfg.fluidVolume) +
          0.5 * m_cfg.rho * Integral(temamJacobian2ALE).over(m_cfg.fluidVolume)

          // Nonlinear Carreau-Yasuda viscous tangent.
          + 2.0 * Integral(mu * symDU, symV).over(m_cfg.fluidVolume) +
          2.0 * Integral(dmu * symU, symV).over(m_cfg.fluidVolume)

          // Stokes pressure/divergence block.
          - Integral(m_p, Div(m_v)).over(m_cfg.fluidVolume) +
          Integral(Div(m_u), m_q).over(m_cfg.fluidVolume) +
          m_cfg.eps * Integral(m_p, m_q).over(m_cfg.fluidVolume)

          // Inlet impedance / tangential / backflow tangents (boundary).
          + m_cfg.inletImpedance *
                BoundaryIntegral(Dot(Dot(m_u, normal) * normal, m_v))
                    .over(m_cfg.inlet)
          + m_cfg.inletTangentialDamping *
                BoundaryIntegral(Dot(duTangential, m_v)).over(m_cfg.inlet)
          +
          BoundaryIntegral(inletBackflowDamping * Dot(m_u, m_v)).over(m_cfg.inlet)
          +
          BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet0) +
          BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet1) +
          BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet2) +
          BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet3) +
          BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet4) +
          BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet5)

          /*
           * =========================
           * Fluid Newton residual (uses the state uState, pState)
           * =========================
           */
          + (m_cfg.rho / m_cfg.dt) * Integral(uState, m_v).over(m_cfg.fluidVolume) -
          (m_cfg.rho / m_cfg.dt) * Integral(m_uOld, m_v).over(m_cfg.fluidVolume)

          + m_cfg.rho *
                Integral(Dot(stateConvectionALE, m_v)).over(m_cfg.fluidVolume)

          + 0.5 * m_cfg.rho * Integral(temamResidualALE).over(m_cfg.fluidVolume)

          + 2.0 * Integral(mu * symU, symV).over(m_cfg.fluidVolume)

          - Integral(pState, Div(m_v)).over(m_cfg.fluidVolume) +
          Integral(Div(uState), m_q).over(m_cfg.fluidVolume) +
          m_cfg.eps * Integral(pState, m_q).over(m_cfg.fluidVolume)

          + BoundaryIntegral(pin * Dot(m_v, normal)).over(m_cfg.inlet)
          + m_cfg.inletImpedance *
                BoundaryIntegral(Dot(Dot(uState, normal) * normal, m_v))
                    .over(m_cfg.inlet)
          + m_cfg.inletTangentialDamping *
                BoundaryIntegral(Dot(uStateTangential, m_v)).over(m_cfg.inlet)

          + BoundaryIntegral(outletPressure(outlet0) * Dot(m_v, normal)).over(outlet0) +
          BoundaryIntegral(outletPressure(outlet1) * Dot(m_v, normal)).over(outlet1) +
          BoundaryIntegral(outletPressure(outlet2) * Dot(m_v, normal)).over(outlet2) +
          BoundaryIntegral(outletPressure(outlet3) * Dot(m_v, normal)).over(outlet3) +
          BoundaryIntegral(outletPressure(outlet4) * Dot(m_v, normal)).over(outlet4) +
          BoundaryIntegral(outletPressure(outlet5) * Dot(m_v, normal)).over(outlet5)

          + BoundaryIntegral(inletBackflowDamping * Dot(uState, m_v)).over(m_cfg.inlet)
          + BoundaryIntegral(outletBackflowDamping * Dot(uState, m_v)).over(outlet0) +
          BoundaryIntegral(outletBackflowDamping * Dot(uState, m_v)).over(outlet1) +
          BoundaryIntegral(outletBackflowDamping * Dot(uState, m_v)).over(outlet2) +
          BoundaryIntegral(outletBackflowDamping * Dot(uState, m_v)).over(outlet3) +
          BoundaryIntegral(outletBackflowDamping * Dot(uState, m_v)).over(outlet4) +
          BoundaryIntegral(outletBackflowDamping * Dot(uState, m_v)).over(outlet5)

          /*
           * =========================
           * Solid NeoHookean elastodynamics over the solid (Newton form)
           * =========================
           */
          // Inertia + Rayleigh damping tangent.
          + (solidRho / (m_cfg.dt * m_cfg.dt)) *
                Integral(m_d, m_w).over(m_cfg.solidVolume)
          + (rayleighAlpha * solidRho / m_cfg.dt) *
                Integral(m_d, m_w).over(m_cfg.solidVolume)
          + (rayleighBeta / m_cfg.dt) * lambdaS *
                Integral(Div(m_d), Div(m_w)).over(m_cfg.solidVolume)
          + (rayleighBeta / m_cfg.dt) * 2.0 * muS *
                Integral(symD, symW).over(m_cfg.solidVolume)
          // NeoHookean material tangent.
          + solidTangent.over(m_cfg.solidVolume)

          // Inertia + Rayleigh damping residual.
          + (solidRho / (m_cfg.dt * m_cfg.dt)) *
                Integral(dState, m_w).over(m_cfg.solidVolume)
          - 2.0 * (solidRho / (m_cfg.dt * m_cfg.dt)) *
                Integral(m_dOld, m_w).over(m_cfg.solidVolume)
          + (solidRho / (m_cfg.dt * m_cfg.dt)) *
                Integral(m_dOldOld, m_w).over(m_cfg.solidVolume)
          + (rayleighAlpha * solidRho / m_cfg.dt) *
                Integral(dState, m_w).over(m_cfg.solidVolume)
          - (rayleighAlpha * solidRho / m_cfg.dt) *
                Integral(m_dOld, m_w).over(m_cfg.solidVolume)
          + (rayleighBeta / m_cfg.dt) * lambdaS *
                Integral(Div(dState), Div(m_w)).over(m_cfg.solidVolume)
          + (rayleighBeta / m_cfg.dt) * 2.0 * muS *
                Integral(symDState, symW).over(m_cfg.solidVolume)
          - (rayleighBeta / m_cfg.dt) * lambdaS *
                Integral(Div(m_dOld), Div(m_w)).over(m_cfg.solidVolume)
          - (rayleighBeta / m_cfg.dt) * 2.0 * muS *
                Integral(symDOld, symW).over(m_cfg.solidVolume)
          // NeoHookean internal force residual.
          + solidInternal.over(m_cfg.solidVolume)

          /*
           * ALE harmonic extension of the displacement into the fluid.
           */
          + m_cfg.aleStiffness *
                Integral(Jacobian(m_d), Jacobian(m_w)).over(m_cfg.fluidVolume)
          + m_cfg.aleStiffness *
                Integral(Jacobian(dState), Jacobian(m_w)).over(m_cfg.fluidVolume)

          /*
           * Weak FSI interface coupling on Gamma_fsi (no DirichletBC).
           * Penalized kinematic constraint t = gammaTilde (u - (d - dOld)/dt),
           * action on fluid (+t.v), reaction on solid (-t.w), Newton form.
           */
          // Fluid tangent.
          + gammaTilde * InterfaceIntegral(Dot(m_u, m_v)).over(m_cfg.fsi)
          - (gammaTilde / m_cfg.dt) * InterfaceIntegral(Dot(m_d, m_v)).over(m_cfg.fsi)
          // Fluid residual.
          + gammaTilde * InterfaceIntegral(Dot(uState, m_v)).over(m_cfg.fsi)
          - (gammaTilde / m_cfg.dt) * InterfaceIntegral(Dot(dState, m_v)).over(m_cfg.fsi)
          + (gammaTilde / m_cfg.dt) * InterfaceIntegral(Dot(m_dOld, m_v)).over(m_cfg.fsi)
          // Solid reaction tangent.
          - gammaTilde * InterfaceIntegral(Dot(m_u, m_w)).over(m_cfg.fsi)
          + (gammaTilde / m_cfg.dt) * InterfaceIntegral(Dot(m_d, m_w)).over(m_cfg.fsi)
          // Solid reaction residual.
          - gammaTilde * InterfaceIntegral(Dot(uState, m_w)).over(m_cfg.fsi)
          + (gammaTilde / m_cfg.dt) * InterfaceIntegral(Dot(dState, m_w)).over(m_cfg.fsi)
          - (gammaTilde / m_cfg.dt) * InterfaceIntegral(Dot(m_dOld, m_w)).over(m_cfg.fsi)

          /*
           * Inactive-domain regularization for the fluid fields in the solid.
           */
          + inactivePenalty * Integral(m_u, m_v).over(m_cfg.solidVolume)
          + inactivePenalty * Integral(uState, m_v).over(m_cfg.solidVolume)
          + inactivePenalty * Integral(m_p, m_q).over(m_cfg.solidVolume)
          + inactivePenalty * Integral(pState, m_q).over(m_cfg.solidVolume)

          /*
           * Solid ring clamps and fluid cap anchors as exact Dirichlet rows
           * in Newton correction form: dd = -dState  =>  d + dd = 0.
           *
           * These replace the former 1e12 penalty terms, which dominated the
           * matrix conditioning and set the achievable residual floor.
           */
          + DirichletBC(m_d, -dState)
                .on(m_cfg.solidRings[0], m_cfg.solidRings[1],
                    m_cfg.solidRings[2], m_cfg.solidRings[3],
                    m_cfg.solidRings[4], m_cfg.solidRings[5])

          + DirichletBC(m_d, -dState)
                .on(m_cfg.inlet, outlet0, outlet1, outlet2, outlet3, outlet4,
                    outlet5);
    }
    else
    {
      m_flow =
          /*
           * =========================
           * Oseen linear operator
           * =========================
           *
           * Unknowns are directly:
           *   m_u : new velocity
           *   m_p : new pressure
           *
           * There is no separate nonlinear residual/tangent split here.
           */

          (m_cfg.rho / m_cfg.dt) * Integral(m_u, m_v).over(m_cfg.fluidVolume)

          /*
           * Oseen convection with lagged convecting velocity uLag = m_uOld.
           */
          + m_cfg.rho * Integral(Dot(oseenConvectionJacobian, m_v)).over(m_cfg.fluidVolume)

          /*
           * Lagged Temam correction.
           */
          + 0.5 * m_cfg.rho * Integral(oseenTemamJacobian).over(m_cfg.fluidVolume)

          /*
           * Lagged viscosity.
           */
          + 2.0 * Integral(muLag * symDU, symV).over(m_cfg.fluidVolume)

          /*
           * Stokes pressure/divergence block.
           */
          - Integral(m_p, Div(m_v)).over(m_cfg.fluidVolume)
          + Integral(Div(m_u), m_q).over(m_cfg.fluidVolume)
          + m_cfg.eps * Integral(m_p, m_q).over(m_cfg.fluidVolume)

          /*
           * VMS bilinear contribution:
           *
           *   int_K tau_K rho^2
           *     ((grad u^{n+1}) u^n)
           *     ·
           *     ((grad v) u^n).
           */

          // + VMSConvectionBilinearIntegrator(m_u, m_v, m_uOld, m_tau.getSolution(),
          //                                   m_cfg.rho)

          /*
           * VMS linear contribution subtracted from the residual:
           *
           *   - int_K rho
           *       (tau_K rho u_proj + u'^{n+1})
           *       ·
           *       ((grad v) u^n).
           */

          // - VMSConvectionLinearIntegrator(m_v, m_sub.getSolution(), m_uOld,
          //                                 m_up.getSolution(), m_tau.getSolution(),
          //                                 m_cfg.rho, m_cfg.dt)

          /*
           * Inlet normal impedance.
           */
          + m_cfg.inletImpedance *
                BoundaryIntegral(Dot(Dot(m_u, normal) * normal, m_v))
                    .over(m_cfg.inlet)

          /*
           * Inlet tangential damping.
           */
          + m_cfg.inletTangentialDamping *
                BoundaryIntegral(Dot(duTangential, m_v)).over(m_cfg.inlet)

          /*
           * Backflow stabilization.
           *
           * Keep only one inlet term.
           */
          + BoundaryIntegral(inletBackflowDamping * Dot(m_u, m_v)).over(m_cfg.inlet)

          + BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet0)
          + BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet1)
          + BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet2)
          + BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet3)
          + BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet4)
          + BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet5)

          /*
           * =========================
           * Oseen right-hand side
           * =========================
           *
           * The old velocity term appears with a minus sign because the problem
           * is assembled in residual form:
           *
           *   A u^{n+1} - rho/dt u^n + boundary loads = 0
           */

          - (m_cfg.rho / m_cfg.dt) * Integral(m_uOld, m_v).over(m_cfg.fluidVolume)

          /*
           * Inlet pressure source.
           */
          + BoundaryIntegral(pin * Dot(m_v, normal)).over(m_cfg.inlet)

          /*
           * Explicit outlet pressures.
           */
          + BoundaryIntegral(outletPressure(outlet0) * Dot(m_v, normal)).over(outlet0)
          + BoundaryIntegral(outletPressure(outlet1) * Dot(m_v, normal)).over(outlet1)
          + BoundaryIntegral(outletPressure(outlet2) * Dot(m_v, normal)).over(outlet2)
          + BoundaryIntegral(outletPressure(outlet3) * Dot(m_v, normal)).over(outlet3)
          + BoundaryIntegral(outletPressure(outlet4) * Dot(m_v, normal)).over(outlet4)
          + BoundaryIntegral(outletPressure(outlet5) * Dot(m_v, normal)).over(outlet5)

          // SOLID BLOCK ---------------------------------------------------
          /*
           * Linear elastodynamics (Hooke) over the solid volume.
           *
           *   rho_s/dt^2 (d^{n+1} - 2 d^n + d^{n-1}, w)
           *   + lambda_s (div d, div w) + 2 mu_s (eps(d), eps(w)) = 0.
           *
           * Second-order inertia uses the two previous displacements
           * d^n = m_dOld and d^{n-1} = m_dOldOld.
           */
          + (solidRho / (m_cfg.dt * m_cfg.dt)) *
                Integral(m_d, m_w).over(m_cfg.solidVolume)

          // Mass-proportional Rayleigh damping:  (alpha rho_s / dt)(d - dOld, w).
          + (rayleighAlpha * solidRho / m_cfg.dt) *
                Integral(m_d, m_w).over(m_cfg.solidVolume)

          + lambdaS * Integral(Div(m_d), Div(m_w)).over(m_cfg.solidVolume)

          + 2.0 * muS * Integral(symD, symW).over(m_cfg.solidVolume)

          // Stiffness-proportional Rayleigh damping:  (beta/dt) K (d - dOld).
          + (rayleighBeta / m_cfg.dt) * lambdaS *
                Integral(Div(m_d), Div(m_w)).over(m_cfg.solidVolume)

          + (rayleighBeta / m_cfg.dt) * 2.0 * muS *
                Integral(symD, symW).over(m_cfg.solidVolume)

          - 2.0 * (solidRho / (m_cfg.dt * m_cfg.dt)) *
                Integral(m_dOld, m_w).over(m_cfg.solidVolume)

          - (rayleighAlpha * solidRho / m_cfg.dt) *
                Integral(m_dOld, m_w).over(m_cfg.solidVolume)

          - (rayleighBeta / m_cfg.dt) * lambdaS *
                Integral(Div(m_dOld), Div(m_w)).over(m_cfg.solidVolume)

          - (rayleighBeta / m_cfg.dt) * 2.0 * muS *
                Integral(symDOld, symW).over(m_cfg.solidVolume)

          + (solidRho / (m_cfg.dt * m_cfg.dt)) *
                Integral(m_dOldOld, m_w).over(m_cfg.solidVolume)

          /*
           * Weak FSI interface coupling on Gamma_fsi (no DirichletBC).
           *
           * Penalized kinematic constraint with action-reaction, which also
           * transfers the interface traction between the two subproblems:
           *
           *   t_Gamma = gammaTilde (u^{n+1} - (d^{n+1} - d^n)/dt),
           *   + int_Gamma t_Gamma . v     (fluid sees +t_Gamma)
           *   - int_Gamma t_Gamma . w     (solid sees the reaction -t_Gamma).
           *
           * gammaTilde = fsiPenalty / dt. As fsiPenalty -> infinity this
           * enforces u = d_t on Gamma. The (d^n, .) terms are the lagged loads.
           */
          // NOTE: Gamma_fsi is an INTERIOR interface (fluid tets / solid
          // wedges share these facets), so the coupling must use
          // InterfaceIntegral. BoundaryIntegral only sees exterior faces and
          // would integrate over nothing here.

          // Fluid momentum: + gammaTilde (u, v) - (gammaTilde/dt) (d, v)
          + gammaTilde * InterfaceIntegral(Dot(m_u, m_v)).over(m_cfg.fsi)

          - (gammaTilde / m_cfg.dt) *
                InterfaceIntegral(Dot(m_d, m_v)).over(m_cfg.fsi)

          + (gammaTilde / m_cfg.dt) *
                InterfaceIntegral(Dot(m_dOld, m_v)).over(m_cfg.fsi)

          // Solid reaction: - gammaTilde (u, w) + (gammaTilde/dt) (d, w)
          - gammaTilde * InterfaceIntegral(Dot(m_u, m_w)).over(m_cfg.fsi)

          + (gammaTilde / m_cfg.dt) *
                InterfaceIntegral(Dot(m_d, m_w)).over(m_cfg.fsi)

          - (gammaTilde / m_cfg.dt) *
                InterfaceIntegral(Dot(m_dOld, m_w)).over(m_cfg.fsi)

          /*
           * ALE harmonic extension of the displacement into the fluid.
           *
           * On the fluid, d is the mesh-motion field: a vector-Laplacian
           * extends the interface motion smoothly through the fluid.
           */
          + m_cfg.aleStiffness *
                Integral(Jacobian(m_d), Jacobian(m_w)).over(m_cfg.fluidVolume)

          /*
           * Inactive-domain algebraic regularization for the globally-defined
           * fluid fields where they carry no physical operator.
           */
          + inactivePenalty * Integral(m_u, m_v).over(m_cfg.solidVolume)

          + inactivePenalty * Integral(m_p, m_q).over(m_cfg.solidVolume)

          /*
           * Solid ring clamps and fluid cap anchors as exact Dirichlet rows.
           *
           * These replace the former 1e12 penalty terms, which dominated the
           * matrix conditioning and set the achievable residual floor.
           */
          + DirichletBC(m_d, Zero(m_mesh.getSpaceDimension()))
                .on(m_cfg.solidRings[0], m_cfg.solidRings[1],
                    m_cfg.solidRings[2], m_cfg.solidRings[3],
                    m_cfg.solidRings[4], m_cfg.solidRings[5])

          + DirichletBC(m_d, Zero(m_mesh.getSpaceDimension()))
                .on(m_cfg.inlet, outlet0, outlet1, outlet2, outlet3, outlet4,
                    outlet5);
    }

    m_stepTiming.setup3DForm = secondsSince(setup3DStart);

    const bool isNewtonFlow = m_cfg.flowMode == FlowMode::Newton;
    const bool needsInitialSystem = !m_flowFieldSplitsSet;

    if (needsInitialSystem)
    {
      if (isRoot())
        ThreeDInfo() << "Initializing flow system ..." << Alert::Raise;
    }

    if (!isNewtonFlow || needsInitialSystem)
    {
      const auto assemble3DStart = CoronaryClock::now();
      m_flow.assemble();
      m_stepTiming.assemble3D = secondsSince(assemble3DStart);
    }

    if (needsInitialSystem)
    {
      const auto fieldSplitsStart = CoronaryClock::now();
      m_flow.setFieldSplits();
      m_stepTiming.fieldSplits = secondsSince(fieldSplitsStart);
      m_flowFieldSplitsSet = true;
    }

    if (isRoot())
    {
      ThreeDInfo() << "Solving with PETSc " << (isNewtonFlow ? "SNES" : "KSP")
                   << " ..." << Alert::Raise;
    }

    const auto solve3DStart = CoronaryClock::now();
    if (isNewtonFlow)
      m_flowSolver.solve();
    else
      m_flow.solve(m_flowKSP);
    m_stepTiming.solve3D = secondsSince(solve3DStart);

    if (isNewtonFlow)
    {
      if (isRoot())
      {
        SNESInfo() << (m_flowSolver.converged() ? "Converged"
                                                : "Did NOT converge")
                   << "  iterations = " << m_flowSolver.getIterationNumber()
                   << Alert::Raise;
      }
      return m_flowSolver.converged();
    }

    ::KSPConvergedReason reason;
    PetscErrorCode ierr = KSPGetConvergedReason(m_flowKSP.getHandle(), &reason);
    assert(ierr == PETSC_SUCCESS);

    PetscInt iterations = 0;
    ierr = KSPGetIterationNumber(m_flowKSP.getHandle(), &iterations);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;

    if (isRoot()) {
      KSPInfo() << (reason > 0 ? "Converged" : "Did NOT converge")
                << "  iterations = " << iterations << Alert::Raise;
    }
    return reason > 0;
  }

  void CoupledLV0DCoronary3D::computeFluxes()
  {
    const auto normal = BoundaryNormal(m_mesh);
    const auto &s = m_model.getState();

    m_flux = BoundaryIntegral(Dot(m_u.getSolution(), normal), m_qFlux)
                 .over(m_cfg.inlet);

    m_flux.assemble();
    m_stepData.qIn = m_flux(m_one);

    m_stepData.qOut.clear();
    m_stepData.qOutSum = 0.0;

    for (const Attribute tag : m_cfg.outlets)
    {
      m_flux =
          BoundaryIntegral(Dot(m_u.getSolution(), normal), m_qFlux).over(tag);

      m_flux.assemble();

      const Real qOut = m_flux(m_one);

      m_stepData.qOut[tag] = qOut;
      m_stepData.qOutSum += qOut;
    }

    m_stepData.flowBalance = m_stepData.qIn + m_stepData.qOutSum;
    m_stepData.t = s.t;
  }

  void CoupledLV0DCoronary3D::printFSIDiagnostics(int step)
  {
    auto maxAbs = [](const auto &gf) -> Real
    {
      PetscReal v = 0.0;
      PetscErrorCode ierr = VecNorm(gf.getData(), NORM_INFINITY, &v);
      assert(ierr == PETSC_SUCCESS);
      (void)ierr;
      return v;
    };

    const Real maxU = maxAbs(m_u.getSolution());
    const Real maxP = maxAbs(m_p.getSolution());
    const Real maxD = maxAbs(m_d.getSolution());

    // Max ALE mesh velocity:  max| (d - dOld)/dt |. A blow-up here (tracking
    // |u|) is the signature of an added-mass / lagged-geometry instability.
    Real maxMeshVel = 0.0;
    {
      ::Vec dtmp = PETSC_NULLPTR;
      PetscErrorCode ierr = VecDuplicate(m_d.getSolution().getData(), &dtmp);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecWAXPY(dtmp, -1.0, m_dOld.getData(), m_d.getSolution().getData());
      assert(ierr == PETSC_SUCCESS);
      PetscReal v = 0.0;
      ierr = VecNorm(dtmp, NORM_INFINITY, &v);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecDestroy(&dtmp);
      assert(ierr == PETSC_SUCCESS);
      (void)ierr;
      maxMeshVel = v / m_cfg.dt;
    }

    // Subdomain volumes on the (moved) mesh. Growth or collapse signals mesh
    // inflation / element tangling from the ALE motion.
    const size_t D = m_mesh.getDimension();
    const Real volFluid = m_mesh.getMeasure(D, m_cfg.fluidVolume);
    const Real volSolid = m_mesh.getMeasure(D, m_cfg.solidVolume);

    // Interface area:  int_fsi 1. Gamma_fsi is an interior interface, so this
    // MUST use InterfaceIntegral (BoundaryIntegral returns 0 here). A non-zero
    // value confirms the coupling actually sees the interface facets.
    m_flux = InterfaceIntegral(m_one, m_qFlux).over(m_cfg.fsi);
    m_flux.assemble();
    const Real areaFsi = m_flux(m_one);

    if (isRoot())
    {
      Alert::Info() << "[FSI diag] step " << step
                    << "  area(fsi) = " << areaFsi << "  vol(f) = " << volFluid
                    << "  vol(s) = " << volSolid << "  max|u| = " << maxU
                    << "  max|meshVel| = " << maxMeshVel
                    << "  max|p| = " << maxP << "  max|d| = " << maxD
                    << Alert::Raise;
    }
  }

  void CoupledLV0DCoronary3D::updateHistory()
  {
    m_uOld.setData(m_u.getSolution().getData());
    m_pOld.setData(m_p.getSolution().getData());
    m_dOldOld.setData(m_dOld.getData());
    m_dOld.setData(m_d.getSolution().getData());

    m_subOld.setData(m_sub.getSolution().getData());
    m_tauOld.setData(m_tau.getSolution().getData());
  }

  void CoupledLV0DCoronary3D::writeOutputs()
  {
    const auto &cy = m_cfg.viscosity;
    const auto symU = 0.5 * (Jacobian(m_u.getSolution()) +
                             Transpose(Jacobian(m_u.getSolution())));
    const auto gamma = Sqrt(cy.gammaRegularization * cy.gammaRegularization +
                            2.0 * Dot(symU, symU));
    const auto carreauBase = 1.0 + Pow(cy.lambda * gamma, cy.yasuda);
    const auto mu = cy.muInf + (cy.mu0 - cy.muInf) *
                                   Pow(carreauBase, (cy.n - 1.0) / cy.yasuda);

    // Fluid-only projection; zero in the solid (solid mass term, no RHS).
    m_viscosityProjection = Integral(m_mu, m_r).over(m_cfg.fluidVolume) -
                            Integral(mu, m_r).over(m_cfg.fluidVolume) +
                            Integral(m_mu, m_r).over(m_cfg.solidVolume);
    m_viscosityProjection.solve(m_viscosityProjectionKSP);

    m_xdmf.write(m_model.getState().t - m_cfg.dt).flush();
  }

  CoupledLV0DCoronary3D::StepData CoupledLV0DCoronary3D::collectStepData() const
  {
    StepData d;

    const auto &s = m_model.getState();

    d.t = s.t;
    d.pat = m_input.pAt(s.t);
    d.psv = m_input.pSv(s.t);

    d.y = s.y;
    d.v = s.v;
    d.radius = m_input.R0 + s.y;

    d.lvVolume =
        (4.0 / 3.0) * std::numbers::pi_v<Real> * d.radius * d.radius * d.radius;

    d.lvFlow = 4.0 * std::numbers::pi_v<Real> * d.radius * d.radius * s.v;

    d.pv = s.pv;
    d.par = s.par;
    d.pd = s.pd;

    d.ec = s.ec;
    d.gamma = s.gamma;
    d.beta = s.beta;
    d.w = s.w;
    d.kc = s.kc;
    d.tauc = s.tauc;

    d.qIn = m_stepData.qIn;
    d.qOutSum = m_stepData.qOutSum;
    d.qDistalSum = 0.0;
    d.qCapChargingSum = 0.0;
    d.flowBalance = m_stepData.flowBalance;
    d.qOut = m_stepData.qOut;

    for (const auto &[tag, bc] : m_wk)
    {
      const auto qOutIt = m_stepData.qOut.find(tag);
      const Real qOut = (qOutIt == m_stepData.qOut.end()) ? 0.0 : qOutIt->second;

      d.qDistal[tag] = bc.qd;
      d.qDistalSum += bc.qd;
      d.qCapChargingSum += qOut - bc.qd;
      d.pc[tag] = bc.pc;
      d.pOut[tag] = bc.pout;
    }

    return d;
  }

  void CoupledLV0DCoronary3D::writeCSVHeader()
  {
    if (!isRoot())
      return;

    m_csv << "t,"
          << "LeftAtriumPressure,"
          << "VenaCavaPressure,"
          << "LeftVentricleDisplacement,"
          << "LeftVentricleVelocity,"
          << "LeftVentricleRadius,"
          << "LeftVentricleVolume,"
          << "LeftVentriclePressure,"
          << "AortaPressure,"
          << "DistalPressure,"
          << "LeftVentricleFlow,"
          << "CoronaryInletFlux,"
          << "CoronaryOutlet4Flux,"
          << "CoronaryOutlet5Flux,"
          << "CoronaryOutlet6Flux,"
          << "CoronaryOutlet7Flux,"
          << "CoronaryOutlet8Flux,"
          << "CoronaryOutlet9Flux,"
          << "CoronaryOutletFluxTotal,"
          << "CoronaryOutlet4DistalFlux,"
          << "CoronaryOutlet5DistalFlux,"
          << "CoronaryOutlet6DistalFlux,"
          << "CoronaryOutlet7DistalFlux,"
          << "CoronaryOutlet8DistalFlux,"
          << "CoronaryOutlet9DistalFlux,"
          << "CoronaryDistalFluxTotal,"
          << "CoronaryCapChargingFluxTotal,"
          << "CoronaryOutlet4CapPressure,"
          << "CoronaryOutlet5CapPressure,"
          << "CoronaryOutlet6CapPressure,"
          << "CoronaryOutlet7CapPressure,"
          << "CoronaryOutlet8CapPressure,"
          << "CoronaryOutlet9CapPressure,"
          << "CoronaryOutlet4Pressure,"
          << "CoronaryOutlet5Pressure,"
          << "CoronaryOutlet6Pressure,"
          << "CoronaryOutlet7Pressure,"
          << "CoronaryOutlet8Pressure,"
          << "CoronaryOutlet9Pressure,"
          << "FlowBalance,"
          << "ec,"
          << "gamma,"
          << "beta,"
          << "w,"
          << "kc,"
          << "tauc\n";

    m_csv.flush();
  }

  void CoupledLV0DCoronary3D::writeCSVRow()
  {
    if (!isRoot())
      return;

    const StepData d = collectStepData();

    auto get = [](const std::map<Attribute, Real> &m, Attribute a) -> Real
    {
      const auto it = m.find(a);
      return (it == m.end()) ? 0.0 : it->second;
    };

    m_csv << d.t << ',' << d.pat << ',' << d.psv << ',' << d.y << ',' << d.v
          << ',' << d.radius << ',' << d.lvVolume << ',' << d.pv << ',' << d.par
          << ',' << d.pd << ',' << d.lvFlow << ',' << d.qIn << ','
          << get(d.qOut, 4) << ',' << get(d.qOut, 5) << ',' << get(d.qOut, 6)
          << ',' << get(d.qOut, 7) << ',' << get(d.qOut, 8) << ','
          << get(d.qOut, 9) << ',' << d.qOutSum << ',' << get(d.qDistal, 4) << ','
          << get(d.qDistal, 5) << ',' << get(d.qDistal, 6) << ','
          << get(d.qDistal, 7) << ',' << get(d.qDistal, 8) << ','
          << get(d.qDistal, 9) << ',' << d.qDistalSum << ',' << d.qCapChargingSum
          << ',' << get(d.pc, 4) << ',' << get(d.pc, 5) << ',' << get(d.pc, 6)
          << ',' << get(d.pc, 7) << ',' << get(d.pc, 8) << ',' << get(d.pc, 9)
          << ',' << get(d.pOut, 4) << ',' << get(d.pOut, 5) << ','
          << get(d.pOut, 6) << ',' << get(d.pOut, 7) << ',' << get(d.pOut, 8)
          << ',' << get(d.pOut, 9) << ',' << d.flowBalance << ',' << d.ec << ','
          << d.gamma << ',' << d.beta << ',' << d.w << ',' << d.kc << ','
          << d.tauc << '\n';

    m_csv.flush();
  }

  void CoupledLV0DCoronary3D::printStepTiming(int step) const
  {
    const auto &comm = m_mesh.getContext().getCommunicator();
    auto maxTime = [&](Real local)
    {
      return boost::mpi::all_reduce(comm, local, MaxReal{});
    };

    const Real total = maxTime(m_stepTiming.total);
    const Real advance0D = maxTime(m_stepTiming.advance0D);
    const Real setup3DForm = maxTime(m_stepTiming.setup3DForm);
    const Real assemble3D = maxTime(m_stepTiming.assemble3D);
    const Real fieldSplits = maxTime(m_stepTiming.fieldSplits);
    const Real solve3D = maxTime(m_stepTiming.solve3D);
    const Real fluxes = maxTime(m_stepTiming.fluxes);
    const Real outletRCR = maxTime(m_stepTiming.outletRCR);
    const Real csv = maxTime(m_stepTiming.csv);
    const Real history = maxTime(m_stepTiming.history);
    const Real output = maxTime(m_stepTiming.output);

    if (isRoot())
    {
      TimingInfo() << "step = " << step << "  total = " << total << " s"
                   << "  0D = " << advance0D << " s"
                   << "  3D-form = " << setup3DForm << " s"
                   << "  3D-assemble = " << assemble3D << " s"
                   << "  field-splits = " << fieldSplits << " s"
                   << "  3D-solve = " << solve3D << " s"
                   << "  fluxes = " << fluxes << " s"
                   << "  RCR = " << outletRCR << " s"
                   << "  csv = " << csv << " s"
                   << "  history = " << history << " s"
                   << "  output = " << output << " s" << Alert::Raise;
    }
  }

  int CoupledLV0DCoronary3D::run()
  {
    if (!m_initialized)
      initialize();

    const Real baseDt = m_cfg.dt;
    const Real factor = m_cfg.timeAdaptivityReductionFactor;
    const int maxLevels = m_cfg.timeAdaptivityMaxLevels;
    const Real startTime = m_model.getState().t;
    const Real finalTime = startTime + static_cast<Real>(m_cfg.nsteps) * baseDt;

    if (!(factor > 0.0 && factor < 1.0))
    {
      Alert::Exception() << "Invalid time adaptivity reduction factor: " << factor
                         << ". Expected a value in (0, 1)." << Alert::Raise;
      return 1;
    }

    if (maxLevels < 0)
    {
      Alert::Exception() << "Invalid time adaptivity maximum levels: "
                         << maxLevels << ". Expected a non-negative value."
                         << Alert::Raise;
      return 1;
    }

    Real nextDt = baseDt;
    int acceptedStep = 0;

    {
      // solveStatic();

      // computeFluxes();
      // for (const Attribute tag : m_cfg.outlets)
      // updateRCRNonNew(m_cfg, m_model, m_wk[tag], m_stepData.qOut.at(tag),
      //  m_cfg.dt);

      // updateHistory();
    }

    while (m_model.getState().t < finalTime - 0.5 * std::numeric_limits<Real>::epsilon())
    {
      const Real t_current = m_model.getState().t;
      const Real remaining = finalTime - t_current;
      const Real physicalDt = std::min(nextDt, remaining);
      Real solverDt = physicalDt;
      int level = 0;

      const auto savedRCR = m_wk;
      const StepData savedStepData = m_stepData;
      const VecSnapshot savedU(m_u.getSolution().getData());
      const VecSnapshot savedP(m_p.getSolution().getData());
      const VecSnapshot savedD(m_d.getSolution().getData());

      m_cfg.dt = physicalDt;
      m_stepTiming = StepTiming{};
      const auto stepStart = CoronaryClock::now();

      if (isRoot())
      {
        Alert::Info() << "━━━ Step " << (acceptedStep + 1)
                      << "  t = " << t_current << " s"
                      << " / " << finalTime << " s"
                      << "  (dt = " << physicalDt << " s)"
                      << " ━━━" << Alert::Raise;
      }

      const auto advance0DStart = CoronaryClock::now();
      const bool advanced0D = advance0D();
      m_stepTiming.advance0D = secondsSince(advance0DStart);

      if (!advanced0D)
      {
        Alert::Exception() << "0D solver failed to converge at step "
                           << (acceptedStep + 1)
                           << ", t = " << m_model.getState().t << Alert::Raise;
        return 1;
      }

      const auto advancedState = m_model.getState();
      const auto advancedHistory = m_model.getHistory();
      const auto advancedUnknowns = m_model.getUnknowns();
      const auto advancedReport = m_model.getReport();

      bool accepted = false;

      // ALE: freeze the geometry for this step at the previous converged
      // displacement. The linear Oseen system is assembled on this mesh.
      moveMeshToDisplacement(m_dOld);

      while (!accepted)
      {
        m_cfg.dt = solverDt;

        if (!solve3D())
        {
          m_model.restore(advancedState, advancedHistory, advancedUnknowns,
                          advancedReport);
          m_wk = savedRCR;
          m_stepData = savedStepData;
          savedU.restore(m_u.getSolution().getData());
          savedP.restore(m_p.getSolution().getData());
          savedD.restore(m_d.getSolution().getData());

          if (level >= maxLevels)
          {
            Alert::Exception() << "3D flow solver failed to converge at step "
                               << (acceptedStep + 1) << " after " << (level + 1)
                               << " attempt(s). Minimum solver dt = " << solverDt
                               << " s." << Alert::Raise;
            return 1;
          }

          solverDt *= factor;
          ++level;

          if (isRoot())
          {
            Alert::Info() << "[3D] Retrying step " << (acceptedStep + 1)
                          << " with reduced solver dt = " << solverDt
                          << " s  (adapt level = " << level << " / " << maxLevels
                          << ")" << Alert::Raise;
          }

          continue;
        }

        // ALE: advance the mesh to the newly computed displacement so that
        // fluxes and output use the deformed configuration.
        moveMeshToDisplacement(m_d.getSolution());

        const auto fluxesStart = CoronaryClock::now();
        computeFluxes();
        m_stepTiming.fluxes = secondsSince(fluxesStart);

        printFSIDiagnostics(acceptedStep + 1);

        const auto outletRCRStart = CoronaryClock::now();
        m_cfg.dt = physicalDt;
        for (const Attribute tag : m_cfg.outlets)
          updateRCRNonNew(m_cfg, tag, m_model, m_wk[tag], m_stepData.qOut.at(tag),
                          m_cfg.dt);
        m_stepTiming.outletRCR = secondsSince(outletRCRStart);

        const auto csvStart = CoronaryClock::now();
        writeCSVRow();
        m_stepTiming.csv = secondsSince(csvStart);

        const auto historyStart = CoronaryClock::now();
        updateHistory();
        m_stepTiming.history = secondsSince(historyStart);

        const auto outputStart = CoronaryClock::now();
        writeOutputs();
        m_stepTiming.output = secondsSince(outputStart);

        m_stepTiming.total = secondsSince(stepStart);
        printStepTiming(acceptedStep + 1);

        accepted = true;
        ++acceptedStep;
        nextDt = std::min(baseDt, solverDt / factor);
      }
    }

    m_cfg.dt = baseDt;
    m_xdmf.close();
    return 0;
  }
} // namespace Rodin::Examples::Heart
