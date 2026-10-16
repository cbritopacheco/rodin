/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file AtriumPTAG.cpp
 * @brief Driver for the 3D atrium with the P-T-A-G coagulation model.
 *
 * Run (from the build directory, after `git lfs pull` for the mesh):
 *   mpirun -n 8 ./examples/Heart/AtriumPTAG -ptag_dt 1e-3
 *
 * Options (same names as LeftAtrium2DPTAG where the meaning is the same):
 *   -ptag_mesh, -ptag_mesh_scale, -ptag_dt, -ptag_pin, -ptag_flow_cycles,
 *   -ptag_species_cycles, -ptag_output_every, -ptag_vms, -ptag_vms_scale,
 *   -ptag_graddiv_scale, -ptag_pspg_scale, -ptag_quasistatic,
 *   -ptag_pspg_residual, -ptag_species_temam, -ptag_crosswind,
 *   -ptag_transport_iterative, -ptag_replay, -ptag_replay_flow_cycles,
 *   -ptag_replay_refresh, -ptag_species_substeps, -ptag_kinetics,
 *   -ptag_brinkman, -ptag_drag, -ptag_fibrinogen, -ptag_wall_flux,
 *   -ptag_layer, -ptag_k_a, -ptag_th_in.
 */
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

#include <petscsys.h>
#include <petscvec.h>

#include <Rodin/Alert.h>
#include <Rodin/Configure.h>

#ifdef RODIN_USE_SCOTCH
#include <Rodin/Scotch/MeshPartitioner.h>
#endif

#include "AtriumPTAG.h"
#include "CoronaryArtery/CoronaryArteryAlerts.h"
#include "CoronaryArtery/CoronaryArteryTiming.h"

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

    void setPrefixedDefault(
      const std::string& prefix, const char* suffix, const char* value)
    {
      const std::string name = "-" + prefix + suffix;
      PetscBool set = PETSC_FALSE;
      PetscErrorCode ierr =
        PetscOptionsHasName(PETSC_NULLPTR, PETSC_NULLPTR, name.c_str(), &set);
      assert(ierr == PETSC_SUCCESS);
      if (!set)
      {
        ierr = PetscOptionsSetValue(PETSC_NULLPTR, name.c_str(), value);
        assert(ierr == PETSC_SUCCESS);
      }
      (void)ierr;
    }

    /// @brief The projections invert mass matrices: CG + Jacobi, never the
    ///        global direct solver used by the coupled flow system.
    void configureMassSolver(Rodin::Solver::KSP& ksp, const std::string& prefix)
    {
      setPrefixedDefault(prefix, "ksp_type", "cg");
      setPrefixedDefault(prefix, "pc_type", "jacobi");
      ksp.setPrefix(prefix);
    }
  }

  // ==========================================================================
  // Construction
  // ==========================================================================
  AtriumPTAG::AtriumPTAG(const Context::MPI& context, const Config& cfg)
    : m_cfg(cfg),
      m_input(makeInput(m_cfg)),
      m_model(m_input),
      m_mesh(makeMesh(context, m_cfg)),
      m_xdmf(context.getCommunicator(), m_cfg.xdmfBasename),
      m_inletSet(m_cfg.inlets.begin(), m_cfg.inlets.end()),
      m_vh(std::integral_constant<size_t, 1>{}, m_mesh, m_mesh.getSpaceDimension()),
      m_sh(std::integral_constant<size_t, 1>{}, m_mesh),
      m_u(m_vh),
      m_p(m_sh),
      m_v(m_vh),
      m_q(m_sh),
      m_uOld(m_vh),
      m_sTrial(m_sh),
      m_sTest(m_sh),
      m_wTrial(m_vh),
      m_wTest(m_vh),
      m_tauFn([this](const Point& p) { return vmsTauAt(p); }),
      m_tauCFn([this](const Point& p) { return tauCAt(p); }),
      m_sqrtTauCFn([this](const Point& p) { return sqrtTauCAt(p); }),
      m_tauPFn([this](const Point& p) { return tauPAt(p); }),
      m_dragFn([this](const Point& p) { return dragAt(p); }),
      m_pspgDragFn([this](const Point& p) { return tauPAt(p) * dragAt(p); }),
      m_piTilde(m_sh),
      m_convProjection(m_vh),
      m_sub(m_vh),
      m_subOld(m_vh),
      m_pT(m_sh), m_tT(m_sh), m_aT(m_sh), m_gT(m_sh),
      m_vP(m_sh), m_vT(m_sh), m_vA(m_sh), m_vG(m_sh),
      m_pCur(m_sh), m_tCur(m_sh), m_aCur(m_sh), m_gCur(m_sh),
      m_pPrev(m_sh), m_tPrev(m_sh), m_aPrev(m_sh), m_gPrev(m_sh),
      m_pNext(m_sh), m_tNext(m_sh), m_aNext(m_sh), m_gNext(m_sh),
      m_fibrin(m_sh), m_tat(m_sh), m_drag(m_sh), m_drift(m_sh),
      m_eT(m_sh), m_vE(m_sh), m_layer(m_sh),
      m_wallArea(m_sh),
      m_wssRec(m_vh),
      m_wss(m_vh),
      m_gradRec0(m_vh),
      m_gradRec1(m_vh),
      m_gradRec2(m_vh),
      m_netShear(m_vh),
      m_absShear(m_sh),
      m_shearMagnitude(m_sh),
      m_tawss(m_sh),
      m_osi(m_sh),
      m_activation(m_sh),
      m_qFlux(m_sh),
      m_one(m_sh),
      m_flux(m_qFlux),
      m_flow(m_u, m_p, m_v, m_q),
      m_flowKSP(m_flow),
      m_transportP(m_pT, m_vP), m_transportT(m_tT, m_vT),
      m_transportA(m_aT, m_vA), m_transportG(m_gT, m_vG),
      m_kspP(m_transportP), m_kspT(m_transportT), m_kspA(m_transportA),
      m_kspG(m_transportG),
      m_extension(m_eT, m_vE),
      m_extensionKSP(m_extension),
      m_scalarProjection(m_sTrial, m_sTest),
      m_scalarProjectionKSP(m_scalarProjection),
      m_vectorProjection(m_wTrial, m_wTest),
      m_vectorProjectionKSP(m_vectorProjection),
      m_wssTrial(m_vh),
      m_wssTest(m_vh),
      m_wssProjection(m_wssTrial, m_wssTest),
      m_wssKSP(m_wssProjection)
  {
    // Queried on every rank, printed on one: nothing that may reduce over the
    // communicator belongs inside a root-only branch.
    const auto cellCount = m_mesh.getCellCount();
    const auto vertexCount = m_mesh.getVertexCount();
    const auto velocityDOFs = m_vh.getSize();
    const auto pressureDOFs = m_sh.getSize();

    if (isRoot())
    {
      Alert::Info() << "[mesh] cells=" << cellCount << " vertices=" << vertexCount
                    << " velocity DOFs=" << velocityDOFs
                    << " pressure DOFs=" << pressureDOFs << Alert::Raise;
    }
  }

  AtriumPTAG::~AtriumPTAG()
  {
    for (auto& v : m_cycleU)
      if (v)
      {
        PetscErrorCode ierr = VecDestroy(&v);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
      }
  }

  bool AtriumPTAG::isRoot() const
  {
    return m_mesh.getContext().getCommunicator().rank() == RootRank;
  }

  AtriumPTAG::MeshType AtriumPTAG::makeMesh(const Context::MPI& context, const Config& cfg)
  {
    const auto& comm = context.getCommunicator();

    Rodin::MPI::Sharder sharder(context);
    if (comm.rank() == RootRank)
    {
      Geometry::Mesh<Context::Local> mesh;
      mesh.load(cfg.meshPath, IO::FileFormat::MEDIT);

      if (mesh.getSpaceDimension() != 3)
        throw std::runtime_error("AtriumPTAG expects a 3D tetrahedral mesh.");

      const size_t D = mesh.getDimension();
      mesh.getConnectivity().compute(D, D);
      mesh.getConnectivity().compute(D, 0);
      mesh.getConnectivity().compute(D, D - 1);
      mesh.getConnectivity().compute(D - 1, D);
      mesh.getConnectivity().compute(D - 1, 0);
      mesh.getConnectivity().compute(D - 1, 1);
      mesh.getConnectivity().compute(1, 0);

#ifdef RODIN_USE_SCOTCH
      Scotch::Partitioner partitioner(mesh);
#else
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

    return mesh;
  }

  AtriumPTAG::Model::Input AtriumPTAG::makeInput(const Config& cfg)
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
    input.localMaxIterations = static_cast<size_t>(cfg.lv.localMaxIterations);
    input.localDamping = cfg.lv.localDamping;
    input.absRegularization = cfg.lv.absRegularization;

    input.pSv = [p = cfg.lv.systemicVenousPressure](Real) { return p; };
    input.pAt = [p = cfg.atrialPressure](Real t) { return atrialWave(p, t); };
    input.u = [a = cfg.activation](Real t) { return activationWave(a, t); };
    input.m0 = [low = cfg.lv.relaxationM0Low, high = cfg.lv.relaxationM0High,
                 lowEc = cfg.lv.relaxationM0LowEc,
                 highEc = cfg.lv.relaxationM0HighEc](Real ec) {
      if (highEc <= lowEc || ec >= highEc)
        return high;
      if (ec <= lowEc)
        return low;
      const Real s = (ec - lowEc) / (highEc - lowEc);
      return (1.0 - s) * low + s * high;
    };
    input.dm0 = [low = cfg.lv.relaxationM0Low, high = cfg.lv.relaxationM0High,
                  lowEc = cfg.lv.relaxationM0LowEc,
                  highEc = cfg.lv.relaxationM0HighEc](Real ec) {
      if (highEc <= lowEc || ec <= lowEc || ec >= highEc)
        return 0.0;
      return (high - low) / (highEc - lowEc);
    };

    {
      using PassiveEnergy = std::decay_t<decltype(input.passiveEnergy)>;
      typename PassiveEnergy::Parameters hp;
      hp.C0 = cfg.lv.passiveC0;
      hp.C1 = cfg.lv.passiveC1;
      hp.C2 = cfg.lv.passiveC2;
      hp.C3 = cfg.lv.passiveC3;
      input.passiveEnergy = PassiveEnergy(hp);
    }

    return input;
  }

  AtriumPTAG::Real AtriumPTAG::activationWave(const Activation& cfg, Real t)
  {
    const Real T = cfg.period;
    const Real tau = t - T * std::floor(t / T);

    auto ss = [](Real s) {
      s = s < 0.0 ? 0.0 : (s > 1.0 ? 1.0 : s);
      return s * s * (3.0 - 2.0 * s);
    };

    if (tau < cfg.tRampStart)
      return 0.0;
    if (tau < cfg.tRampEnd)
      return cfg.positiveValue *
        ss((tau - cfg.tRampStart) / (cfg.tRampEnd - cfg.tRampStart));
    if (tau < cfg.tPlateauEnd)
      return cfg.positiveValue;
    if (tau < cfg.tRelaxEnd)
      return cfg.positiveValue +
        (cfg.negativeValue - cfg.positiveValue) *
        ss((tau - cfg.tPlateauEnd) / (cfg.tRelaxEnd - cfg.tPlateauEnd));
    if (tau < cfg.tNegativeEnd)
      return cfg.negativeValue;
    return cfg.negativeValue *
      (1.0 - ss((tau - cfg.tNegativeEnd) / (T - cfg.tNegativeEnd)));
  }

  AtriumPTAG::Real AtriumPTAG::atrialWave(const AtrialPressure& cfg, Real t)
  {
    const Real T = cfg.period;
    const Real tau = t - T * std::floor(t / T);

    auto ss = [](Real s) {
      s = s < 0.0 ? 0.0 : (s > 1.0 ? 1.0 : s);
      return s * s * (3.0 - 2.0 * s);
    };
    auto ramp = [&](Real a, Real b, Real s) { return a + (b - a) * ss(s); };

    if (tau < cfg.t1)
      return ramp(cfg.minValue, cfg.maxValue, tau / cfg.t1);
    if (tau < cfg.t2)
      return cfg.maxValue;
    if (tau < cfg.t3)
      return ramp(cfg.maxValue, cfg.minValue, (tau - cfg.t2) / (cfg.t3 - cfg.t2));
    if (tau < cfg.t4)
      return ramp(cfg.minValue, cfg.secondThreshold, (tau - cfg.t3) / (cfg.t4 - cfg.t3));
    if (tau < cfg.t5)
      return cfg.secondThreshold;
    if (tau < cfg.t6)
      return ramp(cfg.secondThreshold, cfg.minValue, (tau - cfg.t5) / (cfg.t6 - cfg.t5));
    return cfg.minValue;
  }

  AtriumPTAG::Real AtriumPTAG::cellSize(const Point& p)
  {
    return std::pow(p.getPolytope().getMeasure(), 1.0 / p.getPolytope().getDimension());
  }

  void AtriumPTAG::axpy(Real a, const ::Vec& x, ::Vec& y)
  {
    PetscErrorCode ierr = VecAXPY(y, a, x);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGhostUpdateBegin(y, INSERT_VALUES, SCATTER_FORWARD);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGhostUpdateEnd(y, INSERT_VALUES, SCATTER_FORWARD);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
  }

  // ==========================================================================
  // Pointwise coefficients
  // ==========================================================================
  AtriumPTAG::Real AtriumPTAG::viscosityAt(const Point& p) const
  {
    const auto& cy = m_cfg.viscosity;

    // Built from the same expression the momentum equation uses, so the two
    // cannot drift apart.
    const auto sym = 0.5 * (Jacobian(m_uOld) + Transpose(Jacobian(m_uOld)));
    const Real shear = std::sqrt(
      cy.gammaRegularization * cy.gammaRegularization + 2.0 * Dot(sym, sym).getValue(p));

    return cy.muInf +
      (cy.mu0 - cy.muInf) *
      std::pow(1.0 + std::pow(cy.lambda * shear, cy.yasuda), (cy.n - 1.0) / cy.yasuda);
  }

  AtriumPTAG::Real AtriumPTAG::fibrinAt(const Point& p) const
  {
    const Real G0 = m_cfg.kinetics.fibrinogen0;
    return std::clamp<Real>(G0 - m_gCur.getValue(p), 0.0, G0);
  }

  AtriumPTAG::Real AtriumPTAG::dragAt(const Point& p) const
  {
    return m_cfg.brinkman(fibrinAt(p));
  }

  // tau_1 = [ 4 nu/h^2 + 2|u|/h + sigma_B/rho ]^-1: the Brinkman drag is a
  // reaction and enters the stabilisation scale like one.
  AtriumPTAG::Real AtriumPTAG::tau1At(const Point& p) const
  {
    const auto uc = m_uOld.getValue(p);
    const Real h = cellSize(p);
    const Real nu = viscosityAt(p) / m_cfg.rho;
    return 1.0 / (4.0 * nu / (h * h) + 2.0 * std::sqrt(Math::dot(uc, uc)) / h +
                  dragAt(p) / m_cfg.rho);
  }

  AtriumPTAG::Real AtriumPTAG::vmsTauAt(const Point& p) const
  {
    if (!m_cfg.useVMS)
      return 0.0;
    if (m_cfg.quasiStaticSubscales)
      return m_cfg.vmsScale * tau1At(p) / m_cfg.rho;
    return m_cfg.vmsScale / (m_cfg.rho / m_cfg.dt + m_cfg.rho / tau1At(p));
  }

  AtriumPTAG::Real AtriumPTAG::sqrtTauCAt(const Point& p) const
  {
    if (!m_cfg.useVMS)
      return 0.0;
    const Real h = cellSize(p);
    return std::sqrt(m_cfg.gradDivScale * m_cfg.rho * h * h / (4.0 * tau1At(p)));
  }

  AtriumPTAG::Real AtriumPTAG::tauCAt(const Point& p) const
  {
    const Real s = sqrtTauCAt(p);
    return s * s;
  }

  AtriumPTAG::Real AtriumPTAG::tauPAt(const Point& p) const
  {
    return m_cfg.pspgScale * tau1At(p) / m_cfg.rho;
  }

  AtriumPTAG::Real AtriumPTAG::crosswind(const Point& p, Real cur, Real prev,
    const Math::SpatialVector<Real>& gradient, Real diffusivity) const
  {
    const Real gn = std::sqrt(Math::dot(gradient, gradient));
    if (gn < 1.0e-14)
      return 0.0;
    const auto uc = m_uOld.getValue(p);
    const Real residual = (cur - prev) / m_dtChem + Math::dot(uc, gradient);
    return std::max<Real>(0.0,
      m_cfg.crosswindC * std::abs(residual) * cellSize(p) / (2.0 * gn) - diffusivity);
  }

  // ==========================================================================
  // Setup
  // ==========================================================================
  AtriumPTAG& AtriumPTAG::initialize()
  {
    if (m_cfg.speciesSubsteps < 1)
      throw std::runtime_error("speciesSubsteps must be >= 1");
    if (m_cfg.flowCycles < 1)
      throw std::runtime_error("flowCycles must be >= 1: the activation needs one closed cycle");
    m_dtChem = m_cfg.speciesSubsteps * m_cfg.dt;

    setupSpaces();
    setupFlow();
    setupSpecies();
    setupWallShear();
    setupLayer();

    if (isRoot())
    {
      m_csv.open(m_cfg.csvPath);
      if (!m_csv)
        throw std::runtime_error("Failed to open " + m_cfg.csvPath);
      writeCSVHeader();
    }

    m_initialized = true;
    return *this;
  }

  void AtriumPTAG::setupSpaces()
  {
    Model::State s0;
    s0.t = 0.0;
    s0.y = m_cfg.lv.initialY;
    s0.v = m_cfg.lv.initialV;
    s0.pv = m_input.pAt(0.0) + m_cfg.lv.initialPvOffset;
    s0.par = m_cfg.lv.initialPar;
    s0.pd = m_cfg.lv.initialPd;
    s0.w = m_input.m0(s0.ec);
    m_model.setMaxIterations(m_cfg.lv.maxIterations)
      .setAbsoluteTolerance(m_cfg.lv.absoluteTolerance)
      .setRelativeTolerance(m_cfg.lv.relativeTolerance)
      .setStepTolerance(m_cfg.lv.stepTolerance)
      .setDampingFactor(m_cfg.lv.dampingFactor);
    m_model.initialize(s0);

    const auto zeroVector = Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};
    for (auto* g : { &m_uOld, &m_subOld, &m_sub, &m_convProjection, &m_wss, &m_wssRec,
                     &m_netShear, &m_gradRec0, &m_gradRec1, &m_gradRec2 })
      *g = zeroVector;
    for (auto* g : { &m_piTilde, &m_absShear, &m_shearMagnitude, &m_tawss, &m_osi,
                     &m_activation, &m_fibrin, &m_tat, &m_drag, &m_drift, &m_layer })
      *g = Real(0);
    m_one = Real(1);

    // Inflow composition everywhere at t = 0: no thrombin, full pools.
    const auto& k = m_cfg.kinetics;
    for (auto* g : { &m_pCur, &m_pPrev, &m_pNext }) *g = Real(k.prothrombin0);
    for (auto* g : { &m_tCur, &m_tPrev, &m_tNext }) *g = Real(k.inletThrombin);
    for (auto* g : { &m_aCur, &m_aPrev, &m_aNext }) *g = Real(k.antithrombin0);
    for (auto* g : { &m_gCur, &m_gPrev, &m_gNext }) *g = Real(k.fibrinogen0);

    configureMassSolver(m_scalarProjectionKSP, "ptag_sproj_");
    configureMassSolver(m_vectorProjectionKSP, "ptag_vproj_");
    configureMassSolver(m_wssKSP, "ptag_wss_");
    configureMassSolver(m_extensionKSP, "ptag_layer_");

    // The four transport systems inherit the global direct solver unless
    // asked for GMRES + block-Jacobi/ILU, which is enough for these
    // well-conditioned SUPG systems and much cheaper than MUMPS.
    const auto configureTransport = [this](Rodin::Solver::KSP& ksp, const char* prefix) {
      if (m_cfg.iterativeTransport)
      {
        setPrefixedDefault(prefix, "ksp_type", "gmres");
        setPrefixedDefault(prefix, "pc_type", "bjacobi");
        setPrefixedDefault(prefix, "ksp_rtol", "1e-10");
      }
      ksp.setPrefix(prefix);
    };
    configureTransport(m_kspP, "ptag_P_");
    configureTransport(m_kspT, "ptag_T_");
    configureTransport(m_kspA, "ptag_A_");
    configureTransport(m_kspG, "ptag_G_");

    m_u.setName("velocity");
    m_p.setName("pressure");
    m_xdmf.setMesh(m_mesh);
    m_xdmf.add("velocity", m_u.getSolution());
    m_xdmf.add("pressure", m_p.getSolution());
    m_xdmf.add("prothrombin", m_pCur);
    m_xdmf.add("thrombin", m_tCur);
    m_xdmf.add("antithrombin", m_aCur);
    m_xdmf.add("fibrinogen", m_gCur);
    m_xdmf.add("fibrin", m_fibrin);
    m_xdmf.add("TAT", m_tat);
    m_xdmf.add("brinkmanDrag", m_drag);
    m_xdmf.add("TAWSS", m_tawss);
    m_xdmf.add("OSI", m_osi);
    m_xdmf.add("activation", m_activation);
    m_xdmf.add("activationLayer", m_layer);
    m_xdmf.add("shearStress", m_wss);
    m_xdmf.add("wallArea", m_wallArea);

    m_outletArea = boundaryArea(m_cfg.outlet);

    // Nodal wall indicator from the assembled wall measure (ADD-assembled
    // across ranks, then ghost-refreshed by setData): consistent on every
    // rank, which a Region::Boundary projection is not.
    m_flux = BoundaryIntegral(m_one, m_qFlux).over(m_cfg.wall);
    m_flux.assemble();
    m_wallArea.setData(m_flux.getVector());
    // Interpolating a P1 field at its own node leaves O(1e-21) round-off at
    // the nodes where it is zero, so the indicator is tested against a
    // threshold well below any wall nodal area and well above that noise.
    m_wallAreaThreshold = 1.0e-6 * m_wallArea.max();

    m_pIn = m_cfg.inletPressureMean;
    m_pOut = m_model.getState().pv;
    m_outletPressure = m_pIn;
    m_mitralZ = m_cfg.mitralResistanceOpen * m_outletArea;

    const Real Tth = k.inhibitionRate * k.antithrombin0 * k.amplificationHalf /
      (k.amplificationRate * k.prothrombin0 / k.prothrombinRef -
       k.inhibitionRate * k.antithrombin0);
    if (isRoot())
    {
      Alert::Info() << "[outlet] area = " << m_outletArea << " m^2  pv0 = " << m_pOut
                    << " Pa" << Alert::Raise;
      Alert::Info() << "[kinetics] k_on A_0=" << (k.inhibitionRate * k.antithrombin0)
                    << " 1/s, activation threshold T_th=" << (Tth * 1e6)
                    << " nM, wall layer sqrt(D/(k_on A_0))="
                    << std::sqrt(k.diffusivityProtein / (k.inhibitionRate * k.antithrombin0))
                    << " m, reaction layer delta=" << k.layerThickness << " m"
                    << Alert::Raise;
    }
  }

  bool AtriumPTAG::onWallNode(const Point& p) const
  {
    return m_wallArea.getValue(p) > m_wallAreaThreshold;
  }

  AtriumPTAG::Real AtriumPTAG::boundaryArea(Attribute tag)
  {
    m_flux = BoundaryIntegral(m_one, m_qFlux).over(tag);
    m_flux.assemble();
    return std::max<Real>(m_flux(m_one), 1e-12);
  }

  void AtriumPTAG::setupFlow()
  {
    const size_t dim = m_mesh.getSpaceDimension();
    const auto normal = BoundaryNormal(m_mesh);
    const auto& cy = m_cfg.viscosity;
    const Real deltaMu = cy.mu0 - cy.muInf;
    const Real rho = m_cfg.rho;
    const Real dt = m_cfg.dt;
    const Real pspgR = m_cfg.pspgResidualScale * rho;

    const auto symU = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));
    const auto symV = 0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));
    const auto symLag = 0.5 * (Jacobian(m_uOld) + Transpose(Jacobian(m_uOld)));

    const auto shearLag =
      Sqrt(cy.gammaRegularization * cy.gammaRegularization + 2.0 * Dot(symLag, symLag));
    const auto muLag = cy.muInf +
      deltaMu * Pow(1.0 + Pow(cy.lambda * shearLag, cy.yasuda), (cy.n - 1.0) / cy.yasuda);

    const auto convU = Mult(Jacobian(m_u), m_uOld);
    const auto temam = Div(m_uOld) * Dot(m_u, m_v);

    const auto uNormal = Dot(m_u, normal) * normal;
    const auto uTangential = m_u - uNormal;

    // Directional do-nothing on every pressure boundary (see Atrium.cpp):
    // + (rho/2) max(-u^n.n, 0) u.v cancels the incoming kinetic-energy flux,
    // the only term of the discrete energy balance that is not signed.
    const auto inletBackflow =
      0.5 * rho * m_cfg.inletBackflowStabilization * Max(-Dot(m_uOld, normal), 0.0);
    const auto outletBackflow =
      0.5 * rho * m_cfg.outletBackflowStabilization * Max(-Dot(m_uOld, normal), 0.0);

    // Time-dependent scalars enter through these, so the form is assigned once
    // and only reassembled.
    RealFunction pInFn = [this](const Point&) { return m_pIn; };
    RealFunction pOutFn = [this](const Point&) { return m_pOut; };
    RealFunction mitralFn = [this](const Point&) { return m_mitralZ; };

    // Backward Euler, skew-symmetric convection (Temam), VMS orthogonal
    // subscales, lagged Carreau-Yasuda viscosity, PSPG, and the Brinkman drag
    //   + int sigma_B(F^n) u^{n+1}.v,
    // implicit in u and lagged in F. sigma_B >= 0 makes it a positive
    // semidefinite mass term: it can only remove kinetic energy.
    m_flow = (rho / dt) * Integral(m_u, m_v) - (rho / dt) * Integral(m_uOld, m_v)

      + rho * Integral(Dot(convU, m_v)) + 0.5 * rho * Integral(temam)

      + VMSConvectionBilinearIntegrator(m_u, m_v, m_uOld, m_tauFn, rho) -
      VMSConvectionLinearIntegrator(
        m_v, m_sub, m_uOld, m_convProjection, m_tauFn, rho, dt)

      + VMSGradDivBilinearIntegrator(m_u, m_v, m_tauCFn) -
      VMSGradDivLinearIntegrator(m_v, m_piTilde, m_sqrtTauCFn)

      + 2.0 * Integral(muLag * symU, symV)

      + Integral(m_dragFn * m_u, m_v)

      - Integral(m_p, Div(m_v)) + Integral(Div(m_u), m_q) +
      m_cfg.pressurePenalty * Integral(m_p, m_q)

      // PSPG: tau_p grad q . [ grad p + s (rho (u - u^n)/dt + rho (grad u) u^n
      //                        + sigma_B u) ], s = pspgResidualScale.
      + Integral(m_tauPFn * Grad(m_p), Grad(m_q))
      + (pspgR / dt) * Integral(m_tauPFn * m_u, Grad(m_q))
      - (pspgR / dt) * Integral(m_tauPFn * m_uOld, Grad(m_q))
      + pspgR * Integral(m_tauPFn * convU, Grad(m_q))
      + m_cfg.pspgResidualScale * Integral(m_pspgDragFn * m_u, Grad(m_q))

      + BoundaryIntegral(pInFn * Dot(m_v, normal)).over(m_inletSet) +
      BoundaryIntegral(pOutFn * Dot(m_v, normal)).over(m_cfg.outlet)

      // Source impedance of the pressure inlets and mitral resistance,
      // both positive semidefinite and assembled implicitly.
      + m_cfg.inletImpedance * BoundaryIntegral(Dot(uNormal, m_v)).over(m_inletSet) +
      BoundaryIntegral(mitralFn * Dot(uNormal, m_v)).over(m_cfg.outlet)

      + m_cfg.inletTangentialDamping *
        BoundaryIntegral(Dot(uTangential, m_v)).over(m_inletSet)

      + BoundaryIntegral(inletBackflow * Dot(m_u, m_v)).over(m_inletSet) +
      BoundaryIntegral(outletBackflow * Dot(m_u, m_v)).over(m_cfg.outlet)

      + DirichletBC(m_u, Zero(dim)).on(m_cfg.wall);
  }

  void AtriumPTAG::assignTransport(ScalarProblemType& problem,
    ScalarTrialFunctionType& c, ScalarTestFunctionType& v,
    ScalarGridFunctionType& cur, ScalarGridFunctionType& prev, Real D, Real inlet,
    Real robinScale, const ScalarGridFunctionType* load)
  {
    const Real dt = m_dtChem;
    const Real scale = m_cfg.kinetics.stabilizationScale;
    const Real temam = 0.5 * m_cfg.speciesTemamScale;

    // Cell Peclet ~ 1e6: tau is set by the transient and advective scales.
    RealFunction tauFn = [this, D, scale, dt](const Point& p) {
      const auto uc = m_uOld.getValue(p);
      return scale * supgTau(cellSize(p), std::sqrt(Math::dot(uc, uc)), D, 0.0, dt);
    };
    RealFunction kdcFn = [this, D, &cur, &prev](const Point& p) {
      return crosswind(p, cur.getValue(p), prev.getValue(p), Grad(cur).getValue(p), D);
    };
    RealFunction invSpeedSqFn = [this](const Point& p) {
      const auto uc = m_uOld.getValue(p);
      const Real u2 = Math::dot(uc, uc);
      return (u2 > 1.0e-20) ? 1.0 / u2 : 0.0;
    };
    // Wall exchange j = j_0 E P/P_ref: a Robin sink on P (implicit), the same
    // amount as a load on T (with the P just solved). E is zero until the
    // first cycle has closed.
    RealFunction robinFn = [this, robinScale](const Point& p) {
      return robinScale * m_activation.getValue(p);
    };
    RealFunction loadFn = [this, robinScale, load](const Point& p) {
      return load ? robinScale * m_activation.getValue(p) * load->getValue(p) : 0.0;
    };

    const auto pv = Dot(m_uOld, Grad(v));
    // Temam correction of the transport: a discrete correction, not part of
    // the strong residual, so it stays out of the SUPG and crosswind terms.
    const auto divU = Div(m_uOld);

    problem =
        (1.0 / dt) * Integral(c, v) - (1.0 / dt) * Integral(cur, v)
      + D * Integral(Grad(c), Grad(v))
      + Integral(Dot(m_uOld, Grad(c)), v)
      + temam * Integral(divU * c, v)
      + BoundaryIntegral(robinFn * c, v).over(m_cfg.wall)
      - BoundaryIntegral(loadFn * v).over(m_cfg.wall)

      // SUPG on the transient and advective residual.
      + (1.0 / dt) * Integral(tauFn * c, pv) - (1.0 / dt) * Integral(tauFn * cur, pv)
      + Integral(tauFn * Dot(m_uOld, Grad(c)), pv)

      // Codina crosswind: SUPG leaves the crosswind direction free.
      + Integral(kdcFn * Grad(c), Grad(v))
      - Integral(kdcFn * invSpeedSqFn * Dot(m_uOld, Grad(c)), Dot(m_uOld, Grad(v)))

      // Inflow composition on the pulmonary veins; the mitral outlet is free.
      + DirichletBC(c, RealFunction(Real(inlet))).on(m_inletSet);
  }

  void AtriumPTAG::setupSpecies()
  {
    const auto& k = m_cfg.kinetics;
    // Wall exchange as a Robin flux only when there is no reaction layer;
    // otherwise initiation is the volumetric channel of reactionStep().
    const Real robin = (k.layerThickness > 0.0) ? 0.0 : k.wallFlux / k.prothrombinRef;
    assignTransport(m_transportP, m_pT, m_vP, m_pCur, m_pPrev, k.diffusivityProtein,
      k.prothrombin0, robin, nullptr);
    assignTransport(m_transportT, m_tT, m_vT, m_tCur, m_tPrev, k.diffusivityProtein,
      k.inletThrombin, robin, (k.layerThickness > 0.0) ? nullptr : &m_pCur);
    assignTransport(m_transportA, m_aT, m_vA, m_aCur, m_aPrev, k.diffusivityProtein,
      k.antithrombin0, 0.0, nullptr);
    assignTransport(m_transportG, m_gT, m_vG, m_gCur, m_gPrev, k.diffusivityFibrinogen,
      k.fibrinogen0, 0.0, nullptr);
  }

  void AtriumPTAG::setupWallShear()
  {
    const auto normal = BoundaryNormal(m_mesh);
    const auto& cy = m_cfg.viscosity;
    const auto& uSol = m_u.getSolution();

    const auto symU = 0.5 * (Jacobian(uSol) + Transpose(Jacobian(uSol)));
    const auto shear =
      Sqrt(cy.gammaRegularization * cy.gammaRegularization + 2.0 * Dot(symU, symU));
    const auto mu = cy.muInf +
      (cy.mu0 - cy.muInf) *
        Pow(1.0 + Pow(cy.lambda * shear, cy.yasuda), (cy.n - 1.0) / cy.yasuda);

    // Traction mu (grad u) n from the recovered rows of grad u, projected on
    // the wall only; on a no-slip wall its tangential part equals that of
    // 2 mu eps(u) n.
    const auto traction = VectorFunction(mu * Dot(m_gradRec0, normal),
      mu * Dot(m_gradRec1, normal), mu * Dot(m_gradRec2, normal));
    const auto wallStress = traction - Dot(traction, normal) * normal;

    const Real reg = 1.0e-3;
    m_wssProjection = BoundaryIntegral(Dot(m_wssTrial, m_wssTest)).over(m_cfg.wall) +
      reg * Integral(Dot(m_wssTrial, m_wssTest)) -
      BoundaryIntegral(Dot(wallStress, m_wssTest)).over(m_cfg.wall);
  }

  void AtriumPTAG::setupLayer()
  {
    const Real delta = m_cfg.kinetics.layerThickness;
    if (delta <= 0.0)
      return;
    // Screened Poisson extension of the wall activation: E~ decays into the
    // lumen over the layer thickness; symmetric positive definite, CG.
    m_extension = (delta * delta) * Integral(Grad(m_eT), Grad(m_vE))
                + Integral(m_eT, m_vE)
                + DirichletBC(m_eT, RealFunction([this](const Point& p) -> Real {
                    return m_activation.getValue(p); })).on(m_cfg.wall);
  }

  void AtriumPTAG::extendActivation()
  {
    const Real delta = m_cfg.kinetics.layerThickness;
    if (delta <= 0.0)
      return;
    m_extension.assemble();
    m_extension.solve(m_extensionKSP);
    m_layer.setData(m_eT.getSolution().getData());

    // Dose: int (j_0/delta) scale E~ dV must equal j_0 int_Gw E dA. On a flat
    // wall scale = 1; in pockets narrower than delta (the appendage) the
    // tails overlap and scale < 1.
    m_flux = BoundaryIntegral(m_activation, m_qFlux).over(m_cfg.wall);
    m_flux.assemble();
    const Real wallDose = m_flux(m_one);
    m_flux = Integral(m_layer, m_qFlux);
    m_flux.assemble();
    const Real layerDose = m_flux(m_one);
    m_layerScale = (layerDose > 0.0) ? delta * wallDose / layerDose : 0.0;

    if (isRoot())
      Alert::Info() << "[layer] int E dA = " << wallDose << " m^2, int E~ dV = "
                    << layerDose << " m^3, dose scale = " << m_layerScale
                    << Alert::Raise;
  }

  // ==========================================================================
  // Time stepping
  // ==========================================================================
  bool AtriumPTAG::advance0D()
  {
    const auto rep = m_model.step(m_cfg.dt);
    m_pOut = m_model.getState().pv;

    if (isRoot())
    {
      const auto& s = m_model.getState();
      ZeroDInfo() << (rep.converged ? "converged" : "NOT converged")
                  << "  iter=" << rep.iterations << "  |F|=" << rep.finalResidual
                  << "  |  pv=" << s.pv << " par=" << s.par << " pd=" << s.pd << " Pa"
                  << Alert::Raise;
    }

    return rep.converged;
  }

  bool AtriumPTAG::solveFlow()
  {
    const Real rho = m_cfg.rho;
    const Real dt = m_cfg.dt;

    // The first steps are traced projection by projection: each of these
    // assembles a mass matrix over the whole mesh, so on a fine mesh they are
    // where the time goes, and a stall must be visible where it happens.
    const bool trace = isRoot() && m_step < 3;
    const auto phase = [trace](const char* what) {
      if (trace)
      {
        ThreeDInfo() << what << " ..." << Alert::Raise;
        std::cout.flush();
      }
    };

    const auto vmsStart = CoronaryClock::now();

    if (m_cfg.useVMS)
    {
      phase("VMS: projecting the convective term");
      projectVector(Mult(Jacobian(m_uOld), m_uOld), m_convProjection);

      // Quasi-static: no subscale history, m_sub stays at zero.
      if (!m_cfg.quasiStaticSubscales)
      {
        phase("VMS: projecting the dynamic subscale");
        const size_t dim = m_mesh.getSpaceDimension();
        projectVector(VectorFunction(dim,
                        [this, rho, dt, dim](const Point& p) -> Math::SpatialVector<Real> {
                          const auto conv = Mult(Jacobian(m_uOld), m_uOld).getValue(p);
                          const auto proj = m_convProjection.getValue(p);
                          const auto old = m_subOld.getValue(p);
                          const Real tau = vmsTauAt(p);
                          Math::SpatialVector<Real> out(dim);
                          for (Index c = 0; c < static_cast<Index>(dim); ++c)
                            out(c) = tau * rho * (old(c) / dt - (conv(c) - proj(c)));
                          return out;
                        }),
          m_sub);
      }

      phase("VMS: projecting the grad-div residual");
      project(m_sqrtTauCFn * Div(m_uOld), m_piTilde);
    }

    m_timing.vms = secondsSince(vmsStart);

    if (trace)
    {
      ThreeDInfo() << "Assembling the flow system (" << (m_vh.getSize() + m_sh.getSize())
                   << " unknowns) ..." << Alert::Raise;
      std::cout.flush();
    }

    const auto assemblyStart = CoronaryClock::now();
    m_flow.assemble();
    m_timing.assembly = secondsSince(assemblyStart);

    if (!m_flowFieldSplitsSet)
    {
      if (isRoot())
        ThreeDInfo() << "Creating field splits ..." << Alert::Raise;
      m_flow.setFieldSplits();
      m_flowFieldSplitsSet = true;
    }

    if (trace)
    {
      ThreeDInfo() << "Solving with PETSc KSP"
                   << (m_step == 0 ? " (first solve: the direct factorization "
                                     "is built here and is the slowest one)"
                                   : "")
                   << " ..." << Alert::Raise;
      std::cout.flush();
    }

    const auto solveStart = CoronaryClock::now();
    m_flow.solve(m_flowKSP);
    m_timing.solve = secondsSince(solveStart);

    ::KSPConvergedReason reason;
    PetscErrorCode ierr = KSPGetConvergedReason(m_flowKSP.getHandle(), &reason);
    assert(ierr == PETSC_SUCCESS);
    PetscInt iterations = 0;
    ierr = KSPGetIterationNumber(m_flowKSP.getHandle(), &iterations);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;

    if (isRoot())
      KSPInfo() << (reason > 0 ? "Converged" : "Did NOT converge")
                << "  iterations = " << iterations << "  (" << m_timing.solve << " s)"
                << Alert::Raise;

    m_uOld.setData(m_u.getSolution().getData());
    if (m_cfg.useVMS)
      m_subOld.setData(m_sub.getData());

    // A direct solve always "succeeds": check the physics, not the solver.
    m_speed = std::max(std::abs(m_uOld.max()), std::abs(m_uOld.min()));

    return reason > 0 && std::isfinite(m_speed) && m_speed <= m_cfg.maxVelocity;
  }

  void AtriumPTAG::storeSnapshot(int k)
  {
    PetscErrorCode ierr;
    if (!m_cycleU[k])
    {
      ierr = VecDuplicate(m_uOld.getData(), &m_cycleU[k]);
      assert(ierr == PETSC_SUCCESS);
    }
    ierr = VecCopy(m_uOld.getData(), m_cycleU[k]);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
  }

  void AtriumPTAG::loadSnapshot(int k)
  {
    m_uOld.setData(m_cycleU[k]);
    m_u.getSolution().setData(m_cycleU[k]);
  }

  AtriumPTAG::Real AtriumPTAG::fibrinMass()
  {
    m_flux = Integral(m_fibrin, m_qFlux);
    m_flux.assemble();
    return m_flux(m_one);
  }

  void AtriumPTAG::solveSpecies()
  {
    const auto speciesStart = CoronaryClock::now();

    if (isRoot() && m_step < 3)
    {
      ThreeDInfo() << "Transporting P, T, A, G and reacting ..." << Alert::Raise;
      std::cout.flush();
    }

    // Transport: P first, so that the T load reads P^{n+1}.
    const auto advance = [](ScalarProblemType& problem, Rodin::Solver::KSP& ksp,
      ScalarTrialFunctionType& trial, ScalarGridFunctionType& cur,
      ScalarGridFunctionType& prev) {
      prev.setData(cur.getData());
      problem.assemble();
      problem.solve(ksp);
      cur.setData(trial.getSolution().getData());
    };
    advance(m_transportP, m_kspP, m_pT, m_pCur, m_pPrev);
    advance(m_transportT, m_kspT, m_tT, m_tCur, m_tPrev);
    advance(m_transportA, m_kspA, m_aT, m_aCur, m_aPrev);
    advance(m_transportG, m_kspG, m_gT, m_gCur, m_gPrev);

    reactionStep();

    m_timing.species = secondsSince(speciesStart);
  }

  void AtriumPTAG::reactionStep()
  {
    const auto& k = m_cfg.kinetics;
    const Real dt = m_dtChem;

    // Nodal interpolation of the updated state; every species reads the same
    // pre-reaction values, so the four evaluations are consistent.
    // Layer initiation rate s/P = (j_0/delta) scale E~ / P_ref, a P -> T
    // channel inside the same conservative nodal step.
    const Real layerRate = (k.layerThickness > 0.0)
      ? k.wallFlux / k.layerThickness * m_layerScale / k.prothrombinRef : 0.0;
    const auto state = [this, &k, dt, layerRate](const Point& p) {
      return k.advance({ m_pCur.getValue(p), m_tCur.getValue(p),
                         m_aCur.getValue(p), m_gCur.getValue(p) }, dt,
                       layerRate * std::max<Real>(m_layer.getValue(p), 0.0));
    };
    m_pNext.project(RealFunction([state](const Point& p) -> Real { return state(p).P; }));
    m_tNext.project(RealFunction([state](const Point& p) -> Real { return state(p).T; }));
    m_aNext.project(RealFunction([state](const Point& p) -> Real { return state(p).A; }));
    m_gNext.project(RealFunction([state](const Point& p) -> Real { return state(p).G; }));

    m_pCur.setData(m_pNext.getData());
    m_tCur.setData(m_tNext.getData());
    m_aCur.setData(m_aNext.getData());
    m_gCur.setData(m_gNext.getData());
  }

  void AtriumPTAG::updateDiagnostics()
  {
    const auto& k = m_cfg.kinetics;
    const Real invariant = k.prothrombin0 - k.antithrombin0;
    m_fibrin.project(RealFunction([this](const Point& p) -> Real { return fibrinAt(p); }));
    m_tat.project(RealFunction([this, &k](const Point& p) -> Real {
      return k.antithrombin0 - m_aCur.getValue(p); }));
    m_drag.project(RealFunction([this](const Point& p) -> Real { return dragAt(p); }));
    m_drift.project(RealFunction([this, invariant](const Point& p) -> Real {
      return m_pCur.getValue(p) + m_tCur.getValue(p) - m_aCur.getValue(p) - invariant; }));
  }

  void AtriumPTAG::computeWallShear()
  {
    const auto shearStart = CoronaryClock::now();
    const auto& uSol = m_u.getSolution();

    if (isRoot() && m_step < 3)
    {
      ThreeDInfo() << "Recovering the wall shear stress ..." << Alert::Raise;
      std::cout.flush();
    }

    // Rows of grad u recovered onto the nodes (grad u_h is cellwise constant
    // on P1).
    projectVector(VectorFunction(Component(Jacobian(uSol), 0, 0),
                    Component(Jacobian(uSol), 0, 1), Component(Jacobian(uSol), 0, 2)),
      m_gradRec0);
    projectVector(VectorFunction(Component(Jacobian(uSol), 1, 0),
                    Component(Jacobian(uSol), 1, 1), Component(Jacobian(uSol), 1, 2)),
      m_gradRec1);
    projectVector(VectorFunction(Component(Jacobian(uSol), 2, 0),
                    Component(Jacobian(uSol), 2, 1), Component(Jacobian(uSol), 2, 2)),
      m_gradRec2);

    m_wssProjection.assemble();
    m_wssProjection.solve(m_wssKSP);

    // Keep the wall only: off it the regularised projection carries a
    // sign-alternating tail that would poison the cycle accumulators (see
    // Atrium.cpp for the argument). Projected over the cells with the nodal
    // wall indicator, so every rank writes every local node consistently.
    const size_t sdim = m_mesh.getSpaceDimension();
    m_wssRec.setData(m_wssTrial.getSolution().getData());
    m_wss.project(VectorFunction(sdim,
      [this, sdim](const Point& p) -> Math::SpatialVector<Real> {
        Math::SpatialVector<Real> out(sdim);
        if (!onWallNode(p))
        {
          out.setZero();
          return out;
        }
        const auto w = m_wssRec.getValue(p);
        for (Index c = 0; c < static_cast<Index>(sdim); ++c)
          out(c) = w(c);
        return out;
      }));

    // Cycle accumulators: int tau_w dt (vector) and int |tau_w| dt (nodal
    // magnitude, never an L2 projection, which can go negative).
    axpy(m_cfg.dt, m_wss.getData(), m_netShear.getData());
    m_shearMagnitude.project(Sqrt(Dot(m_wss, m_wss)));
    axpy(m_cfg.dt, m_shearMagnitude.getData(), m_absShear.getData());

    m_timing.shear = secondsSince(shearStart);
  }

  void AtriumPTAG::closeCycle(Real elapsed)
  {
    if (elapsed <= 0.0)
      return;

    // TAWSS, OSI and E live on the wall and nowhere else; all three are nodal
    // and are projected over the cells with the wall indicator (never over
    // Region::Boundary, see m_wallArea).
    m_tawss.project(RealFunction([this, elapsed](const Point& p) -> Real {
      return onWallNode(p) ? m_absShear.getValue(p) / elapsed : 0.0;
    }));

    m_osi.project(RealFunction([this](const Point& p) -> Real {
      if (!onWallNode(p))
        return 0.0;
      const Real abs = m_absShear.getValue(p);
      if (abs <= 0.0)
        return 0.0;
      const auto net = m_netShear.getValue(p);
      const Real mag = std::sqrt(Math::dot(net, net));
      return std::clamp<Real>(0.5 * (1.0 - mag / abs), 0.0, 0.5);
    }));

    // E = sigma_low(TAWSS) * min(2 OSI, 1): a logistic in the measured
    // endothelial switch (~0.4 Pa) times the oscillatory index. Bounded in
    // [0,1], so the wall exchange can never become a sink of thrombin.
    const Real tauA = m_cfg.kinetics.activationShearStress;
    const Real width = std::max<Real>(m_cfg.kinetics.activationShearWidth, 1e-12);
    m_activation.project(RealFunction([this, tauA, width](const Point& p) -> Real {
      if (!onWallNode(p))
        return 0.0;
      const Real low = 1.0 / (1.0 + std::exp((m_tawss.getValue(p) - tauA) / width));
      const Real osi = std::clamp<Real>(2.0 * m_osi.getValue(p), 0.0, 1.0);
      return std::clamp<Real>(low * osi, 0.0, 1.0);
    }));

    // The ghost entries must be zeroed too, or the next accumulation reads a
    // stale halo.
    const auto zeroWithGhosts = [](::Vec& v) {
      PetscErrorCode e = VecZeroEntries(v);
      assert(e == PETSC_SUCCESS);
      e = VecGhostUpdateBegin(v, INSERT_VALUES, SCATTER_FORWARD);
      assert(e == PETSC_SUCCESS);
      e = VecGhostUpdateEnd(v, INSERT_VALUES, SCATTER_FORWARD);
      assert(e == PETSC_SUCCESS);
      (void)e;
    };
    zeroWithGhosts(m_netShear.getData());
    zeroWithGhosts(m_absShear.getData());

    extendActivation();
    m_indicesReady = true;
  }

  void AtriumPTAG::computeFluxes()
  {
    const auto normal = BoundaryNormal(m_mesh);
    const auto& uSol = m_u.getSolution();

    m_flux = BoundaryIntegral(Dot(uSol, normal), m_qFlux).over(m_inletSet);
    m_flux.assemble();
    m_qIn = m_flux(m_one);

    m_flux = BoundaryIntegral(Dot(uSol, normal), m_qFlux).over(m_cfg.outlet);
    m_flux.assemble();
    m_qOut = m_flux(m_one);

    m_flux = BoundaryIntegral(m_p.getSolution(), m_qFlux).over(m_cfg.outlet);
    m_flux.assemble();
    m_outletPressure = m_flux(m_one) / m_outletArea;
  }

  void AtriumPTAG::writeCSVHeader()
  {
    m_csv << "t,cycle,pIn,pv,par,pd,valveOpen,qIn,qOut,maxU,maxTAWSS,maxOSI,maxE,"
             "maxT,minP,minA,maxTAT,maxF,maxDrag,invariantDrift\n";
  }

  void AtriumPTAG::writeCSVRow(int cycle)
  {
    // Collective reductions on every rank, then a root-only write.
    const auto& s = m_model.getState();
    const Real tawss = m_tawss.max(), osi = m_osi.max(), act = m_activation.max();
    const Real T = m_tCur.max(), P = m_pCur.min(), A = m_aCur.min();
    const Real tat = m_tat.max(), F = m_fibrin.max(), drag = m_drag.max();
    const Real drift = std::max(std::abs(m_drift.max()), std::abs(m_drift.min()));
    if (!isRoot())
      return;
    m_csv << m_t << ',' << cycle << ',' << m_pIn << ',' << s.pv << ',' << s.par << ','
          << s.pd << ',' << (m_valveOpen ? 1 : 0) << ',' << m_qIn << ',' << m_qOut << ','
          << m_speed << ',' << tawss << ',' << osi << ',' << act << ',' << T << ',' << P
          << ',' << A << ',' << tat << ',' << F << ',' << drag << ',' << drift << '\n';
    m_csv.flush();
  }

  int AtriumPTAG::run()
  {
    if (!m_initialized)
      initialize();

    const Real PI = std::numbers::pi_v<Real>;
    const int stepsPerCycle = static_cast<int>(m_cfg.period / m_cfg.dt + 0.5);
    const int totalCycles = m_cfg.flowCycles + m_cfg.speciesCycles;
    const int totalSteps = totalCycles * stepsPerCycle;
    const int sub = m_cfg.speciesSubsteps;
    if (stepsPerCycle % sub != 0)
      throw std::runtime_error("speciesSubsteps must divide the steps per cycle");
    if (m_cfg.replay)
      m_cycleU.assign(stepsPerCycle, PETSC_NULLPTR);

    if (isRoot())
      Alert::Info() << "[run] " << totalCycles << " cycles of " << stepsPerCycle
                    << " steps (" << totalSteps << " total); species from cycle "
                    << (m_cfg.flowCycles + 1) << " on (never before the first OSI); "
                    << "Brinkman " << (m_cfg.brinkman.enabled ? "on" : "off")
                    << "; replay " << (m_cfg.replay ? "on" : "off") << "; XDMF every "
                    << m_cfg.outputEvery << " steps and at every cycle boundary"
                    << Alert::Raise;

    Real cycleElapsed = 0.0;
    Real storedFibrinMass = -1.0;
    int refreshCyclesLeft = 0;
    int flowSteps = 0;
    int chemSteps = 0;
    const auto runStart = CoronaryClock::now();

    for (int step = 0; step < totalSteps; ++step)
    {
      m_step = step;
      m_timing = Timing{};
      const auto stepStart = CoronaryClock::now();

      m_t += m_cfg.dt;
      const int cycle = step / stepsPerCycle;
      const int k = step % stepsPerCycle;
      const bool endOfCycle = (k == stepsPerCycle - 1);
      const bool warmup = cycle < m_cfg.flowCycles;
      const bool solveFlowNow = !m_cfg.replay || warmup || refreshCyclesLeft > 0;
      // Species advance on every sub-th step with dt_q = sub dt; a replayed
      // step that carries neither flow nor species costs nothing.
      const bool chemStep = ((step + 1) % sub == 0);

      if (isRoot() && (solveFlowNow || chemStep))
      {
        Alert::Info() << "---- Step " << (step + 1) << "/" << totalSteps << "  cycle "
                      << (cycle + 1) << "/" << totalCycles << "  t = " << m_t
                      << " s  (dt = " << m_cfg.dt << " s"
                      << (warmup ? ", warm-up" : "") << (solveFlowNow ? "" : ", replay")
                      << ") ----" << Alert::Raise;
        std::cout.flush();
      }

      if (solveFlowNow)
      {
        // 1. 0D ventricle, mitral diode and the flow with the drag of the
        //    current fibrin field; wall shear accumulated and the indices
        //    closed once per cycle.
        const auto zeroDStart = CoronaryClock::now();
        const bool advanced = advance0D();
        m_timing.zeroD = secondsSince(zeroDStart);
        if (!advanced)
        {
          Alert::Exception() << "[0D] Newton did not converge at t = " << m_t
                             << Alert::Raise;
          return 1;
        }

        m_pIn = m_cfg.inletPressureMean +
          m_cfg.inletPressureAmplitude * std::sin(2.0 * PI * m_t / m_cfg.period);

        // Mitral diode, switched on the previous step's outlet pressure.
        m_valveOpen = (m_outletPressure > m_pOut);
        m_mitralZ =
          (m_valveOpen ? m_cfg.mitralResistanceOpen : m_cfg.mitralResistanceClosed) *
          m_outletArea;

        if (!solveFlow())
        {
          Alert::Exception() << "[flow] diverged at step " << (step + 1)
                             << ": max|u| = " << m_speed << " m/s" << Alert::Raise;
          return 1;
        }
        ++flowSteps;

        computeWallShear();
        cycleElapsed += m_cfg.dt;

        const auto fluxStart = CoronaryClock::now();
        computeFluxes();
        m_timing.fluxes = secondsSince(fluxStart);

        if (m_cfg.replay)
          storeSnapshot(k);

        if (endOfCycle)
        {
          closeCycle(cycleElapsed);
          cycleElapsed = 0.0;
          if (refreshCyclesLeft > 0)
            --refreshCyclesLeft;

          // VecMax is collective: every rank must reach it.
          const Real maxTawss = m_tawss.max();
          const Real maxOsi = m_osi.max();
          const Real maxActivation = m_activation.max();
          if (isRoot())
            Alert::Info() << "[cycle " << (cycle + 1) << "/" << totalCycles
                          << "] maxTAWSS=" << maxTawss << " Pa  maxOSI=" << maxOsi
                          << "  maxE=" << maxActivation << Alert::Raise;
        }
      }
      else if (chemStep)
      {
        // 1'. Replay: the stored cycle stands in for the flow and the 0D
        //     ventricle; E is frozen.
        loadSnapshot(k);
        computeFluxes();
      }

      // 2. Transport (four linear solves) and nodal reaction.
      if (m_cfg.solveKinetics && m_indicesReady && !warmup && chemStep)
      {
        solveSpecies();
        ++chemSteps;
      }

      if (solveFlowNow || chemStep)
      {
        updateDiagnostics();
        writeCSVRow(cycle);
      }

      // 3. Refresh test, at the end of a replayed cycle: the flow is solved
      //    again only once the gel is dense enough to act on it and has grown
      //    since the stored cycle.
      if (m_cfg.replay && !warmup && endOfCycle && refreshCyclesLeft == 0 &&
          m_cfg.brinkman.enabled)
      {
        const Real fmax = m_fibrin.max();
        const Real mass = fibrinMass();
        if (fmax > 0.5 * m_cfg.brinkman.gelThreshold &&
            (storedFibrinMass < 0.0 ||
             mass > (1.0 + m_cfg.replayRefresh) * storedFibrinMass))
        {
          refreshCyclesLeft = m_cfg.replayFlowCycles;
          storedFibrinMass = mass;
          m_subOld = Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};
          if (isRoot())
            Alert::Info() << "[replay] cycle " << (cycle + 1) << ": max F = " << fmax
                          << ", fibrin mass = " << mass << " -> refreshing the flow for "
                          << refreshCyclesLeft << " cycles" << Alert::Raise;
        }
      }

      const bool output =
        endOfCycle || (m_cfg.outputEvery > 0 && step % m_cfg.outputEvery == 0);
      if (output)
      {
        const auto outputStart = CoronaryClock::now();
        m_xdmf.write(m_t).flush();
        m_timing.output = secondsSince(outputStart);
      }

      m_timing.total = secondsSince(stepStart);

      if (isRoot() && (solveFlowNow || chemStep))
      {
        const Real elapsed = secondsSince(runStart);
        const Real perStep = elapsed / static_cast<Real>(step + 1);

        Alert::Info() << "[3D] max|u|=" << m_speed << " m/s  pIn=" << m_pIn
                      << " Pa  pv=" << m_pOut
                      << " Pa  valve=" << (m_valveOpen ? "open" : "closed")
                      << "  qIn=" << m_qIn << "  qOut=" << m_qOut << " m^3/s"
                      << Alert::Raise;

        Alert::Info() << "[timing] 0D=" << m_timing.zeroD << "  vms=" << m_timing.vms
                      << "  asm=" << m_timing.assembly << "  ksp=" << m_timing.solve
                      << "  wss=" << m_timing.shear << "  flux=" << m_timing.fluxes
                      << "  species=" << m_timing.species << "  out=" << m_timing.output
                      << "  total=" << m_timing.total << " s  |  ETA "
                      << (perStep * (totalSteps - step - 1) / 60.0) << " min"
                      << Alert::Raise;
        std::cout.flush();
      }
    }

    if (isRoot())
      Alert::Info() << "[run] of " << totalSteps << " clock steps: flow solved in "
                    << flowSteps << ", species solved in " << chemSteps
                    << " (dt_q = " << m_dtChem << " s)" << Alert::Raise;

    m_xdmf.close();
    if (isRoot())
      m_csv.close();

    return 0;
  }
}

int main(int argc, char** argv)
{
  PetscInitialize(&argc, &argv, PETSC_NULLPTR, PETSC_NULLPTR);

  const auto setPETScDefault = [](const char* name, const char* value) {
    PetscBool set = PETSC_FALSE;
    PetscErrorCode ierr = PetscOptionsHasName(PETSC_NULLPTR, PETSC_NULLPTR, name, &set);
    if (ierr == PETSC_SUCCESS && !set)
      ierr = PetscOptionsSetValue(PETSC_NULLPTR, name, value);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
  };

  setPETScDefault("-ksp_type", "preonly");
  setPETScDefault("-pc_type", "lu");
  setPETScDefault("-pc_factor_mat_solver_type", "mumps");
  setPETScDefault("-mat_mumps_icntl_20", "0");
  setPETScDefault("-mat_mumps_icntl_21", "0");

  boost::mpi::environment env(argc, argv);
  boost::mpi::communicator world(PETSC_COMM_WORLD, boost::mpi::comm_attach);
  Rodin::Context::MPI context(env, world);

  try
  {
    int status = 0;

    {
      using Rodin::Examples::Heart::AtriumPTAG;
      AtriumPTAG::Config cfg;

      const auto getString = [](const char* name, std::string& out) {
        char buffer[512];
        PetscBool got = PETSC_FALSE;
        PetscOptionsGetString(PETSC_NULLPTR, PETSC_NULLPTR, name, buffer, sizeof(buffer), &got);
        if (got) out = buffer;
      };
      const auto getReal = [](const char* name, Rodin::Real& out) {
        PetscReal value = 0.0;
        PetscBool got = PETSC_FALSE;
        PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, name, &value, &got);
        if (got) out = value;
      };
      const auto getInt = [](const char* name, int& out) {
        PetscInt value = 0;
        PetscBool got = PETSC_FALSE;
        PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, name, &value, &got);
        if (got) out = static_cast<int>(value);
      };
      const auto getBool = [](const char* name, bool& out) {
        PetscBool value = PETSC_FALSE;
        PetscBool got = PETSC_FALSE;
        PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, name, &value, &got);
        if (got) out = (value == PETSC_TRUE);
      };

      // Geometry, forcing and time stepping (Atrium).
      getString("-ptag_mesh", cfg.meshPath);
      getReal("-ptag_mesh_scale", cfg.meshScale);
      getReal("-ptag_dt", cfg.dt);
      getReal("-ptag_pin", cfg.inletPressureMean);
      getInt("-ptag_flow_cycles", cfg.flowCycles);
      getInt("-ptag_species_cycles", cfg.speciesCycles);
      getInt("-ptag_output_every", cfg.outputEvery);

      // Stabilisation.
      getBool("-ptag_vms", cfg.useVMS);
      getReal("-ptag_vms_scale", cfg.vmsScale);
      getReal("-ptag_graddiv_scale", cfg.gradDivScale);
      getReal("-ptag_pspg_scale", cfg.pspgScale);
      getBool("-ptag_quasistatic", cfg.quasiStaticSubscales);
      getReal("-ptag_pspg_residual", cfg.pspgResidualScale);
      getReal("-ptag_species_temam", cfg.speciesTemamScale);
      getReal("-ptag_crosswind", cfg.crosswindC);
      getBool("-ptag_transport_iterative", cfg.iterativeTransport);

      // Replay and species clock.
      getBool("-ptag_replay", cfg.replay);
      getInt("-ptag_replay_flow_cycles", cfg.replayFlowCycles);
      getReal("-ptag_replay_refresh", cfg.replayRefresh);
      getInt("-ptag_species_substeps", cfg.speciesSubsteps);

      // Kinetics and drag.
      getBool("-ptag_kinetics", cfg.solveKinetics);
      getBool("-ptag_brinkman", cfg.brinkman.enabled);
      getReal("-ptag_drag", cfg.brinkman.coefficient);
      getReal("-ptag_fibrinogen", cfg.kinetics.fibrinogen0);
      getReal("-ptag_wall_flux", cfg.kinetics.wallFlux);
      getReal("-ptag_layer", cfg.kinetics.layerThickness);
      getReal("-ptag_k_a", cfg.kinetics.amplificationRate);
      getReal("-ptag_th_in", cfg.kinetics.inletThrombin);

      AtriumPTAG simulation(context, cfg);
      status = simulation.initialize().run();
    }

    PetscFinalize();
    return status;
  }
  catch (const std::exception& e)
  {
    std::cerr << "Fatal error: " << e.what() << "\n";
    PetscFinalize();
    return 1;
  }
}
