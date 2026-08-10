// Atrium.cpp
//
// Run:
//   mpirun -n 8 ./examples/Heart/Atrium -atrium_mesh <path> -atrium_dt 1e-3
//
// Options: -atrium_mesh, -atrium_mesh_scale, -atrium_dt, -atrium_flow_cycles,
//          -atrium_species_cycles, -atrium_pin.
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

#include "Atrium.h"
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

    void setPrefixedDefault(const std::string& prefix, const char* suffix,
      const char* value)
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

  Atrium::Atrium(const Context::MPI& context, const Config& cfg)
    : m_cfg(cfg),
      m_input(makeInput(m_cfg)),
      m_model(m_input),
      m_mesh(makeMesh(context, m_cfg)),
      m_xdmf(context.getCommunicator(), m_cfg.xdmfBasename),
      m_inletSet(m_cfg.inlets.begin(), m_cfg.inlets.end()),
      m_vh(std::integral_constant<size_t, 1>{}, m_mesh, m_mesh.getSpaceDimension()),
      m_sh(std::integral_constant<size_t, 1>{}, m_mesh),
      m_u(m_vh), m_p(m_sh), m_v(m_vh), m_q(m_sh), m_uOld(m_vh),
      m_sTrial(m_sh), m_sTest(m_sh), m_wTrial(m_vh), m_wTest(m_vh),
      m_tauFn([this](const Point& p) { return vmsTauAt(p); }),
      m_tauCFn([this](const Point& p) { return tauCAt(p); }),
      m_sqrtTauCFn([this](const Point& p) { return sqrtTauCAt(p); }),
      m_tauPFn([this](const Point& p) { return tauPAt(p); }),
      m_piTilde(m_sh),
      m_convProjection(m_vh), m_sub(m_vh), m_subOld(m_vh),
      m_th(m_sh), m_fg(m_sh), m_fn(m_sh),
      m_vth(m_sh), m_vfg(m_sh), m_vfn(m_sh),
      m_thCur(m_sh), m_fgCur(m_sh), m_fnCur(m_sh),
      m_thPrev(m_sh), m_fgPrev(m_sh), m_fnPrev(m_sh),
      m_wss(m_vh), m_gradRec0(m_vh), m_gradRec1(m_vh), m_gradRec2(m_vh),
      m_netShear(m_vh), m_absShear(m_sh), m_shearMagnitude(m_sh),
      m_tawss(m_sh), m_osi(m_sh), m_activation(m_sh),
      m_qFlux(m_sh), m_one(m_sh), m_flux(m_qFlux),
      m_flow(m_u, m_p, m_v, m_q),
      m_flowKSP(m_flow),
      m_species(m_th, m_fg, m_fn, m_vth, m_vfg, m_vfn),
      m_speciesKSP(m_species),
      m_scalarProjection(m_sTrial, m_sTest),
      m_scalarProjectionKSP(m_scalarProjection),
      m_vectorProjection(m_wTrial, m_wTest),
      m_vectorProjectionKSP(m_vectorProjection),
      m_wssTrial(m_vh), m_wssTest(m_vh),
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

  Atrium::~Atrium() = default;

  bool Atrium::isRoot() const
  {
    return m_mesh.getContext().getCommunicator().rank() == RootRank;
  }

  Atrium::MeshType Atrium::makeMesh(const Context::MPI& context, const Config& cfg)
  {
    const auto& comm = context.getCommunicator();

    Rodin::MPI::Sharder sharder(context);
    if (comm.rank() == RootRank)
    {
      Geometry::Mesh<Context::Local> mesh;
      mesh.load(cfg.meshPath, IO::FileFormat::MEDIT);

      if (mesh.getSpaceDimension() != 3)
        throw std::runtime_error("Atrium expects a 3D tetrahedral mesh.");

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

  Atrium::Model::Input Atrium::makeInput(const Config& cfg)
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

  Atrium::Real Atrium::activationWave(const Activation& cfg, Real t)
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
      return cfg.positiveValue + (cfg.negativeValue - cfg.positiveValue) *
        ss((tau - cfg.tPlateauEnd) / (cfg.tRelaxEnd - cfg.tPlateauEnd));
    if (tau < cfg.tNegativeEnd)
      return cfg.negativeValue;
    return cfg.negativeValue *
      (1.0 - ss((tau - cfg.tNegativeEnd) / (T - cfg.tNegativeEnd)));
  }

  Atrium::Real Atrium::atrialWave(const AtrialPressure& cfg, Real t)
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

  Atrium::Real Atrium::cellSize(const Point& p)
  {
    return std::pow(p.getPolytope().getMeasure(), 1.0 / p.getPolytope().getDimension());
  }

  Atrium::Real Atrium::viscosityAt(const Point& p) const
  {
    const auto& cy = m_cfg.viscosity;

    // Built from the same expression the momentum equation uses, so the two
    // cannot drift apart.
    const auto sym = 0.5 * (Jacobian(m_uOld) + Transpose(Jacobian(m_uOld)));
    const Real shear = std::sqrt(cy.gammaRegularization * cy.gammaRegularization +
                                 2.0 * Dot(sym, sym).getValue(p));

    return cy.muInf + (cy.mu0 - cy.muInf) *
      std::pow(1.0 + std::pow(cy.lambda * shear, cy.yasuda),
               (cy.n - 1.0) / cy.yasuda);
  }

  Atrium::Real Atrium::tau1At(const Point& p) const
  {
    const auto uc = m_uOld.getValue(p);
    const Real h = cellSize(p);
    const Real nu = viscosityAt(p) / m_cfg.rho;
    return 1.0 / (4.0 * nu / (h * h) + 2.0 * std::sqrt(Math::dot(uc, uc)) / h);
  }

  Atrium::Real Atrium::vmsTauAt(const Point& p) const
  {
    return m_cfg.vmsScale / (m_cfg.rho / m_cfg.dt + m_cfg.rho / tau1At(p));
  }

  Atrium::Real Atrium::sqrtTauCAt(const Point& p) const
  {
    const Real h = cellSize(p);
    return std::sqrt(m_cfg.gradDivScale * m_cfg.rho * h * h / (4.0 * tau1At(p)));
  }

  Atrium::Real Atrium::tauCAt(const Point& p) const
  {
    const Real s = sqrtTauCAt(p);
    return s * s;
  }

  Atrium::Real Atrium::tauPAt(const Point& p) const
  {
    return m_cfg.pspgScale * tau1At(p) / m_cfg.rho;
  }

  void Atrium::axpy(Real a, const ::Vec& x, ::Vec& y)
  {
    PetscErrorCode ierr = VecAXPY(y, a, x);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGhostUpdateBegin(y, INSERT_VALUES, SCATTER_FORWARD);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGhostUpdateEnd(y, INSERT_VALUES, SCATTER_FORWARD);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
  }

  Atrium& Atrium::initialize()
  {
    setupSpaces();
    setupFlow();
    setupSpecies();
    setupWallShear();

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

  void Atrium::setupSpaces()
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

    m_uOld = zeroVector;
    m_subOld = zeroVector;
    m_sub = zeroVector;
    m_convProjection = zeroVector;
    m_wss = zeroVector;
    m_netShear = zeroVector;
    m_gradRec0 = zeroVector;
    m_gradRec1 = zeroVector;
    m_gradRec2 = zeroVector;

    m_piTilde = Real(0);
    m_absShear = Real(0);
    m_shearMagnitude = Real(0);
    m_tawss = Real(0);
    m_osi = Real(0);
    m_activation = Real(0);
    m_one = Real(1);

    const Real fg0 =
      m_cfg.thrombosis.fibrinogenSinusRhythm / m_cfg.thrombosis.fibrinogenMolarMass;
    m_thCur = Real(0);
    m_thPrev = Real(0);
    m_fnCur = Real(0);
    m_fnPrev = Real(0);
    m_fgCur = fg0;
    m_fgPrev = fg0;

    configureMassSolver(m_scalarProjectionKSP, "atrium_sproj_");
    configureMassSolver(m_vectorProjectionKSP, "atrium_vproj_");
    configureMassSolver(m_wssKSP, "atrium_wss_");

    m_u.setName("velocity");
    m_p.setName("pressure");
    m_thCur.setName("thrombin");
    m_fgCur.setName("fibrinogen");
    m_fnCur.setName("fibrin");
    m_tawss.setName("TAWSS");
    m_osi.setName("OSI");
    m_activation.setName("activation");
    m_wss.setName("shearStress");

    m_xdmf.setMesh(m_mesh);
    m_xdmf.add("velocity", m_u.getSolution());
    m_xdmf.add("pressure", m_p.getSolution());
    m_xdmf.add("thrombin", m_thCur);
    m_xdmf.add("fibrinogen", m_fgCur);
    m_xdmf.add("fibrin", m_fnCur);
    m_xdmf.add("TAWSS", m_tawss);
    m_xdmf.add("OSI", m_osi);
    m_xdmf.add("activation", m_activation);
    m_xdmf.add("shearStress", m_wss);

    m_outletArea = boundaryArea(m_cfg.outlet);
    m_pIn = m_cfg.inletPressureMean;
    m_pOut = m_model.getState().pv;
    m_outletPressure = m_pIn;
    m_mitralZ = m_cfg.mitralResistanceOpen * m_outletArea;

    if (isRoot())
      Alert::Info() << "[outlet] area = " << m_outletArea << " m^2  pv0 = "
                    << m_pOut << " Pa" << Alert::Raise;
  }

  Atrium::Real Atrium::boundaryArea(Attribute tag)
  {
    m_flux = BoundaryIntegral(m_one, m_qFlux).over(tag);
    m_flux.assemble();
    return std::max<Real>(m_flux(m_one), 1e-12);
  }

  void Atrium::setupFlow()
  {
    const size_t dim = m_mesh.getSpaceDimension();
    const auto normal = BoundaryNormal(m_mesh);
    const auto& cy = m_cfg.viscosity;
    const Real deltaMu = cy.mu0 - cy.muInf;
    const Real rho = m_cfg.rho;
    const Real dt = m_cfg.dt;

    const auto symU = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));
    const auto symV = 0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));
    const auto symLag = 0.5 * (Jacobian(m_uOld) + Transpose(Jacobian(m_uOld)));

    const auto shearLag = Sqrt(cy.gammaRegularization * cy.gammaRegularization +
                               2.0 * Dot(symLag, symLag));
    const auto muLag = cy.muInf +
      deltaMu * Pow(1.0 + Pow(cy.lambda * shearLag, cy.yasuda),
                    (cy.n - 1.0) / cy.yasuda);

    const auto convU = Mult(Jacobian(m_u), m_uOld);
    const auto temam = Div(m_uOld) * Dot(m_u, m_v);

    const auto uNormal = Dot(m_u, normal) * normal;
    const auto uTangential = m_u - uNormal;

    const auto inletBackflow = 0.5 * rho * m_cfg.inletBackflowStabilization *
      Max(Dot(m_uOld, normal), 0.0);
    const auto outletBackflow = 0.5 * rho * m_cfg.outletBackflowStabilization *
      Max(-Dot(m_uOld, normal), 0.0);

    // Time-dependent scalars enter through these, so the form is assigned once
    // and only reassembled.
    RealFunction pInFn = [this](const Point&) { return m_pIn; };
    RealFunction pOutFn = [this](const Point&) { return m_pOut; };
    RealFunction mitralFn = [this](const Point&) { return m_mitralZ; };

    m_flow = (rho / dt) * Integral(m_u, m_v)
           - (rho / dt) * Integral(m_uOld, m_v)

           + rho * Integral(Dot(convU, m_v))
           + 0.5 * rho * Integral(temam)

           + VMSConvectionBilinearIntegrator(m_u, m_v, m_uOld, m_tauFn, rho)
           - VMSConvectionLinearIntegrator(
               m_v, m_sub, m_uOld, m_convProjection, m_tauFn, rho, dt)

           + VMSGradDivBilinearIntegrator(m_u, m_v, m_tauCFn)
           - VMSGradDivLinearIntegrator(m_v, m_piTilde, m_sqrtTauCFn)

           + 2.0 * Integral(muLag * symU, symV)
           - Integral(m_p, Div(m_v)) + Integral(Div(m_u), m_q)
           + m_cfg.pressurePenalty * Integral(m_p, m_q)

           // PSPG. Required by the equal-order pair.
           + Integral(m_tauPFn * Grad(m_p), Grad(m_q))

           + BoundaryIntegral(pInFn * Dot(m_v, normal)).over(m_inletSet)
           + BoundaryIntegral(pOutFn * Dot(m_v, normal)).over(m_cfg.outlet)

           // Source impedance of the pressure inlets and mitral resistance,
           // both positive semidefinite and assembled implicitly.
           + m_cfg.inletImpedance *
               BoundaryIntegral(Dot(uNormal, m_v)).over(m_inletSet)
           + BoundaryIntegral(mitralFn * Dot(uNormal, m_v)).over(m_cfg.outlet)

           + m_cfg.inletTangentialDamping *
               BoundaryIntegral(Dot(uTangential, m_v)).over(m_inletSet)

           + BoundaryIntegral(inletBackflow * Dot(m_u, m_v)).over(m_inletSet)
           + BoundaryIntegral(outletBackflow * Dot(m_u, m_v)).over(m_cfg.outlet)

           + DirichletBC(m_u, Zero(dim)).on(m_cfg.wall);
  }

  Atrium::Real Atrium::crosswind(const Point& p, Real cur, Real prev,
    const Math::SpatialVector<Real>& gradient, Real diffusivity, Real reaction) const
  {
    const Real gn = std::sqrt(Math::dot(gradient, gradient));
    if (gn < 1.0e-14)
      return 0.0;
    const auto uc = m_uOld.getValue(p);
    const Real residual =
      (cur - prev) / m_cfg.dt + Math::dot(uc, gradient) + reaction;
    return std::max<Real>(
      0.0, m_cfg.crosswindC * std::abs(residual) * cellSize(p) / (2.0 * gn) -
        diffusivity);
  }

  void Atrium::setupSpecies()
  {
    const Real Dth = m_cfg.thrombosis.diffusivityThrombin;
    const Real Dfg = m_cfg.thrombosis.diffusivityFibrinogen;
    const Real Dfn = m_cfg.thrombosis.diffusivityFibrin;
    const Real keff = m_cfg.thrombosis.reactionRate;
    const Real scale = m_cfg.thrombosis.stabilizationScale;
    const Real dt = m_cfg.dt;
    const Real fg0 =
      m_cfg.thrombosis.fibrinogenSinusRhythm / m_cfg.thrombosis.fibrinogenMolarMass;

    // Cell Peclet is of order 1e6 here: these are essentially pure advection
    // problems and their tau is set by the advective and transient scales, not
    // by the flow's tau.
    RealFunction tauThFn = [this, Dth, scale, dt](const Point& p) {
      const auto uc = m_uOld.getValue(p);
      return scale * supgTau(cellSize(p), std::sqrt(Math::dot(uc, uc)), Dth, 0.0, dt);
    };
    RealFunction tauFgFn = [this, Dfg, keff, scale, dt](const Point& p) {
      const auto uc = m_uOld.getValue(p);
      return scale * supgTau(cellSize(p), std::sqrt(Math::dot(uc, uc)), Dfg,
        keff * std::abs(m_thCur.getValue(p)), dt);
    };
    RealFunction tauFnFn = [this, Dfn, scale, dt](const Point& p) {
      const auto uc = m_uOld.getValue(p);
      return scale * supgTau(cellSize(p), std::sqrt(Math::dot(uc, uc)), Dfn, 0.0, dt);
    };

    // Codina crosswind: SUPG leaves the crosswind direction free, which is what
    // produces negative concentration bands across shear layers.
    RealFunction invSpeedSqFn = [this](const Point& p) {
      const auto uc = m_uOld.getValue(p);
      const Real u2 = Math::dot(uc, uc);
      return (u2 > 1.0e-20) ? 1.0 / u2 : 0.0;
    };
    RealFunction kdcThFn = [this, Dth](const Point& p) {
      return crosswind(p, m_thCur.getValue(p), m_thPrev.getValue(p),
        Grad(m_thCur).getValue(p), Dth, 0.0);
    };
    RealFunction kdcFgFn = [this, Dfg, keff](const Point& p) {
      return crosswind(p, m_fgCur.getValue(p), m_fgPrev.getValue(p),
        Grad(m_fgCur).getValue(p), Dfg,
        keff * m_thCur.getValue(p) * m_fgCur.getValue(p));
    };
    RealFunction kdcFnFn = [this, Dfn, keff](const Point& p) {
      return crosswind(p, m_fnCur.getValue(p), m_fnPrev.getValue(p),
        Grad(m_fnCur).getValue(p), Dfn,
        -keff * m_thCur.getValue(p) * m_fgCur.getValue(p));
    };

    const auto pTh = Dot(m_uOld, Grad(m_vth));
    const auto pFg = Dot(m_uOld, Grad(m_vfg));
    const auto pFn = Dot(m_uOld, Grad(m_vfn));

    m_species =
        (1.0 / dt) * Integral(m_th, m_vth) - (1.0 / dt) * Integral(m_thCur, m_vth)
      + Dth * Integral(Grad(m_th), Grad(m_vth))
      + Integral(Dot(m_uOld, Grad(m_th)), m_vth)

      + (1.0 / dt) * Integral(m_fg, m_vfg) - (1.0 / dt) * Integral(m_fgCur, m_vfg)
      + Dfg * Integral(Grad(m_fg), Grad(m_vfg))
      + Integral(Dot(m_uOld, Grad(m_fg)), m_vfg)
      + keff * Integral(m_thCur * m_fg, m_vfg)

      + (1.0 / dt) * Integral(m_fn, m_vfn) - (1.0 / dt) * Integral(m_fnCur, m_vfn)
      + Dfn * Integral(Grad(m_fn), Grad(m_vfn))
      + Integral(Dot(m_uOld, Grad(m_fn)), m_vfn)
      - keff * Integral(m_thCur * m_fg, m_vfn)

      // Endothelial thrombin flux: a surface flux, armed by the cycle indices.
      - m_cfg.thrombosis.thrombinWallFlux *
          BoundaryIntegral(m_activation * m_vth).over(m_cfg.wall)

      // SUPG. Only fibrinogen has a sink proportional to its own unknown.
      + (1.0 / dt) * Integral(tauThFn * m_th, pTh)
      - (1.0 / dt) * Integral(tauThFn * m_thCur, pTh)
      + Integral(tauThFn * Dot(m_uOld, Grad(m_th)), pTh)

      + (1.0 / dt) * Integral(tauFgFn * m_fg, pFg)
      - (1.0 / dt) * Integral(tauFgFn * m_fgCur, pFg)
      + Integral(tauFgFn * Dot(m_uOld, Grad(m_fg)), pFg)
      + keff * Integral(tauFgFn * m_thCur * m_fg, pFg)

      + (1.0 / dt) * Integral(tauFnFn * m_fn, pFn)
      - (1.0 / dt) * Integral(tauFnFn * m_fnCur, pFn)
      + Integral(tauFnFn * Dot(m_uOld, Grad(m_fn)), pFn)
      - keff * Integral(tauFnFn * m_thCur * m_fg, pFn)

      // Codina crosswind, one pair per species.
      + Integral(kdcThFn * Grad(m_th), Grad(m_vth))
      - Integral(kdcThFn * invSpeedSqFn * Dot(m_uOld, Grad(m_th)),
                 Dot(m_uOld, Grad(m_vth)))

      + Integral(kdcFgFn * Grad(m_fg), Grad(m_vfg))
      - Integral(kdcFgFn * invSpeedSqFn * Dot(m_uOld, Grad(m_fg)),
                 Dot(m_uOld, Grad(m_vfg)))

      + Integral(kdcFnFn * Grad(m_fn), Grad(m_vfn))
      - Integral(kdcFnFn * invSpeedSqFn * Dot(m_uOld, Grad(m_fn)),
                 Dot(m_uOld, Grad(m_vfn)))

      // Pure advection needs a Dirichlet condition on the inflow boundary;
      // fresh blood carries plasma fibrinogen and no thrombin or fibrin.
      + DirichletBC(m_th, RealFunction(Real(0))).on(m_inletSet)
      + DirichletBC(m_fn, RealFunction(Real(0))).on(m_inletSet)
      + DirichletBC(m_fg, RealFunction(fg0)).on(m_inletSet);
  }

  void Atrium::setupWallShear()
  {
    const auto normal = BoundaryNormal(m_mesh);
    const auto& cy = m_cfg.viscosity;
    const auto& uSol = m_u.getSolution();

    const auto symU = 0.5 * (Jacobian(uSol) + Transpose(Jacobian(uSol)));
    const auto shear = Sqrt(cy.gammaRegularization * cy.gammaRegularization +
                            2.0 * Dot(symU, symU));
    const auto mu = cy.muInf + (cy.mu0 - cy.muInf) *
      Pow(1.0 + Pow(cy.lambda * shear, cy.yasuda), (cy.n - 1.0) / cy.yasuda);

    const auto traction = VectorFunction(
      mu * Dot(m_gradRec0, normal),
      mu * Dot(m_gradRec1, normal),
      mu * Dot(m_gradRec2, normal));
    const auto wallStress = traction - Dot(traction, normal) * normal;

    const Real reg = 1.0e-3;
    m_wssProjection = BoundaryIntegral(Dot(m_wssTrial, m_wssTest)).over(m_cfg.wall)
                    + reg * Integral(Dot(m_wssTrial, m_wssTest))
                    - BoundaryIntegral(Dot(wallStress, m_wssTest)).over(m_cfg.wall);
  }

  bool Atrium::advance0D()
  {
    const auto rep = m_model.step(m_cfg.dt);
    m_pOut = m_model.getState().pv;

    if (isRoot())
    {
      const auto& s = m_model.getState();
      ZeroDInfo() << (rep.converged ? "converged" : "NOT converged")
                  << "  iter=" << rep.iterations << "  |F|=" << rep.finalResidual
                  << "  |  pv=" << s.pv << " par=" << s.par << " pd=" << s.pd
                  << " Pa" << Alert::Raise;
    }

    return rep.converged;
  }

  bool Atrium::solveFlow()
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

    // Only the two orthogonal-subscale projections remain: Pi[(grad u)u] and
    // Pi[sqrt(tau_C) div u]. Those ARE L2 projections -- subtracting the
    // projection of the residual is what the method does. The taus are not:
    // they are evaluated pointwise by m_tauFn and friends.
    if (m_cfg.useVMS)
    {
      phase("VMS: projecting the convective term");
      projectVector(Mult(Jacobian(m_uOld), m_uOld), m_convProjection);

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
        }), m_sub);

      // Same sqrt(tau_C) that multiplies div(v) in the linear term and whose
      // square is the implicit coefficient. Lagging one of the two by a step,
      // as the reference does, leaves the explicit half larger than the
      // implicit one on the step where the viscosity drops.
      phase("VMS: projecting the grad-div residual");
      project(m_sqrtTauCFn * Div(m_uOld), m_piTilde);
    }

    m_timing.vms = secondsSince(vmsStart);

    if (isRoot())
    {
      ThreeDInfo() << "Assembling the flow system ("
                   << (m_vh.getSize() + m_sh.getSize()) << " unknowns) ..."
                   << Alert::Raise;
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

    if (isRoot())
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
                << "  iterations = " << iterations << "  (" << m_timing.solve
                << " s)" << Alert::Raise;

    m_uOld.setData(m_u.getSolution().getData());
    if (m_cfg.useVMS)
      m_subOld.setData(m_sub.getData());

    // A direct solve always "succeeds": check the physics, not the solver.
    m_speed = std::max(std::abs(m_uOld.max()), std::abs(m_uOld.min()));

    return reason > 0 && std::isfinite(m_speed) && m_speed <= m_cfg.maxVelocity;
  }

  void Atrium::computeWallShear()
  {
    const auto shearStart = CoronaryClock::now();
    const auto& uSol = m_u.getSolution();

    if (isRoot() && m_step < 3)
    {
      ThreeDInfo() << "Recovering the wall shear stress ..." << Alert::Raise;
      std::cout.flush();
    }

    projectVector(VectorFunction(Component(Jacobian(uSol), 0, 0),
      Component(Jacobian(uSol), 0, 1), Component(Jacobian(uSol), 0, 2)), m_gradRec0);
    projectVector(VectorFunction(Component(Jacobian(uSol), 1, 0),
      Component(Jacobian(uSol), 1, 1), Component(Jacobian(uSol), 1, 2)), m_gradRec1);
    projectVector(VectorFunction(Component(Jacobian(uSol), 2, 0),
      Component(Jacobian(uSol), 2, 1), Component(Jacobian(uSol), 2, 2)), m_gradRec2);

    m_wssProjection.assemble();
    m_wssProjection.solve(m_wssKSP);
    m_wss.setData(m_wssTrial.getSolution().getData());

    // Cycle accumulators. Two are needed and they are not interchangeable: the
    // vector integral measures how much net direction survives, the scalar one
    // how much shear was applied regardless of direction.
    axpy(m_cfg.dt, m_wss.getData(), m_netShear.getData());
    project(Sqrt(Dot(m_wss, m_wss)), m_shearMagnitude);
    axpy(m_cfg.dt, m_shearMagnitude.getData(), m_absShear.getData());

    m_timing.shear = secondsSince(shearStart);
  }

  void Atrium::closeCycle(Real elapsed)
  {
    if (elapsed <= 0.0)
      return;

    project(RealFunction([this, elapsed](const Point& p) {
      return m_absShear.getValue(p) / elapsed;
    }), m_tawss);

    project(RealFunction([this](const Point& p) {
      const Real abs = m_absShear.getValue(p);
      const auto net = m_netShear.getValue(p);
      const Real mag = std::sqrt(Math::dot(net, net));
      return (abs > 0.0) ? 0.5 * (1.0 - mag / abs) : 0.0;
    }), m_osi);

    // Smooth, bounded activation: a logistic in the measured shear threshold
    // times the oscillatory index mapped onto [0,1]. Unlike ECAP it does not
    // diverge where TAWSS vanishes.
    const Real tauA = m_cfg.thrombosis.activationShearStress;
    const Real width = std::max<Real>(m_cfg.thrombosis.activationShearWidth, 1e-12);
    project(RealFunction([this, tauA, width](const Point& p) {
      const Real low = 1.0 / (1.0 + std::exp((m_tawss.getValue(p) - tauA) / width));
      return low * std::min<Real>(2.0 * m_osi.getValue(p), 1.0);
    }), m_activation);

    // The ghost entries must be zeroed too, or the next accumulation reads a
    // stale halo.
    const auto zeroWithGhosts = [](::Vec& v) {
      PetscErrorCode ierr = VecZeroEntries(v);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecGhostUpdateBegin(v, INSERT_VALUES, SCATTER_FORWARD);
      assert(ierr == PETSC_SUCCESS);
      ierr = VecGhostUpdateEnd(v, INSERT_VALUES, SCATTER_FORWARD);
      assert(ierr == PETSC_SUCCESS);
      (void)ierr;
    };

    zeroWithGhosts(m_netShear.getData());
    zeroWithGhosts(m_absShear.getData());
  }

  void Atrium::solveSpecies()
  {
    const auto speciesStart = CoronaryClock::now();

    if (isRoot() && m_step < 3)
    {
      ThreeDInfo() << "Assembling and solving the coagulation kinetics ..."
                   << Alert::Raise;
      std::cout.flush();
    }

    m_species.assemble();
    m_species.solve(m_speciesKSP);

    // History rotates before the current level is overwritten: the crosswind
    // residual needs both c^n and c^{n-1}.
    m_thPrev.setData(m_thCur.getData());
    m_fgPrev.setData(m_fgCur.getData());
    m_fnPrev.setData(m_fnCur.getData());
    m_thCur.setData(m_th.getSolution().getData());
    m_fgCur.setData(m_fg.getSolution().getData());
    m_fnCur.setData(m_fn.getSolution().getData());

    m_timing.species = secondsSince(speciesStart);
  }

  void Atrium::computeFluxes()
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

  void Atrium::writeCSVHeader()
  {
    m_csv << "t,cycle,pIn,pv,par,pd,valveOpen,qIn,qOut,maxU,"
          << "maxTAWSS,maxOSI,maxActivation,maxThrombin,minFibrinogen,maxFibrin\n";
  }

  void Atrium::writeCSVRow(int cycle)
  {
    const auto& s = m_model.getState();
    const Real tawss = m_tawss.max();
    const Real osi = m_osi.max();
    const Real activation = m_activation.max();
    const Real th = m_thCur.max();
    const Real fg = m_fgCur.min();
    const Real fn = m_fnCur.max();

    if (!isRoot())
      return;

    m_csv << m_t << ',' << cycle << ',' << m_pIn << ',' << s.pv << ',' << s.par
          << ',' << s.pd << ',' << (m_valveOpen ? 1 : 0) << ',' << m_qIn << ','
          << m_qOut << ',' << m_speed << ',' << tawss << ',' << osi << ','
          << activation << ',' << th << ',' << fg << ',' << fn << '\n';
    m_csv.flush();
  }

  int Atrium::run()
  {
    if (!m_initialized)
      initialize();

    const Real PI = std::numbers::pi_v<Real>;
    const int stepsPerCycle = static_cast<int>(m_cfg.period / m_cfg.dt + 0.5);
    const int totalCycles = m_cfg.flowCycles + m_cfg.speciesCycles;
    const int totalSteps = totalCycles * stepsPerCycle;

    if (isRoot())
      Alert::Info() << "[run] " << totalCycles << " cycles of " << stepsPerCycle
                    << " steps (" << totalSteps << " total); species from cycle "
                    << (m_cfg.flowCycles + 1) << " on; XDMF every "
                    << m_cfg.outputEvery << " steps and at every cycle boundary"
                    << Alert::Raise;

    Real cycleElapsed = 0.0;
    const auto runStart = CoronaryClock::now();

    for (int step = 0; step < totalSteps; ++step)
    {
      m_step = step;
      m_timing = Timing{};
      const auto stepStart = CoronaryClock::now();

      m_t += m_cfg.dt;
      const int cycle = step / stepsPerCycle;
      const bool endOfCycle = (step % stepsPerCycle == stepsPerCycle - 1);

      if (isRoot())
      {
        Alert::Info() << "---- Step " << (step + 1) << "/" << totalSteps
                      << "  cycle " << (cycle + 1) << "/" << totalCycles
                      << "  t = " << m_t << " s  (dt = " << m_cfg.dt << " s"
                      << (cycle < m_cfg.flowCycles ? ", warm-up" : "") << ") ----"
                      << Alert::Raise;
        std::cout.flush();
      }

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
      m_mitralZ = (m_valveOpen ? m_cfg.mitralResistanceOpen
                               : m_cfg.mitralResistanceClosed) * m_outletArea;

      if (!solveFlow())
      {
        Alert::Exception() << "[flow] diverged at step " << (step + 1)
                           << ": max|u| = " << m_speed << " m/s" << Alert::Raise;
        return 1;
      }

      computeWallShear();
      cycleElapsed += m_cfg.dt;

      const auto fluxStart = CoronaryClock::now();
      computeFluxes();
      m_timing.fluxes = secondsSince(fluxStart);

      if (endOfCycle)
      {
        closeCycle(cycleElapsed);
        cycleElapsed = 0.0;

        // VecMax is collective: every rank must reach it, so the reductions
        // happen outside the root-only print.
        const Real maxTawss = m_tawss.max();
        const Real maxOsi = m_osi.max();
        const Real maxActivation = m_activation.max();

        if (isRoot())
          Alert::Info() << "[cycle " << (cycle + 1) << "/" << totalCycles
                        << "] maxTAWSS=" << maxTawss << " Pa  maxOSI=" << maxOsi
                        << "  maxActivation=" << maxActivation << Alert::Raise;
      }

      if (m_cfg.solveKinetics && cycle >= m_cfg.flowCycles)
        solveSpecies();

      writeCSVRow(cycle);

      // The flow is written from the very first step, so the warm-up cycles are
      // available before the species are switched on.
      const bool output = endOfCycle ||
        (m_cfg.outputEvery > 0 && step % m_cfg.outputEvery == 0);
      if (output)
      {
        const auto outputStart = CoronaryClock::now();
        m_xdmf.write(m_t).flush();
        m_timing.output = secondsSince(outputStart);
      }

      m_timing.total = secondsSince(stepStart);

      if (isRoot())
      {
        const Real elapsed = secondsSince(runStart);
        const Real perStep = elapsed / static_cast<Real>(step + 1);

        Alert::Info() << "[3D] max|u|=" << m_speed << " m/s  pIn=" << m_pIn
                      << " Pa  pv=" << m_pOut << " Pa  valve="
                      << (m_valveOpen ? "open" : "closed") << "  qIn=" << m_qIn
                      << "  qOut=" << m_qOut << " m^3/s" << Alert::Raise;

        Alert::Info() << "[timing] 0D=" << m_timing.zeroD << "  vms=" << m_timing.vms
                     << "  asm=" << m_timing.assembly << "  ksp=" << m_timing.solve
                     << "  wss=" << m_timing.shear << "  flux=" << m_timing.fluxes
                     << "  species=" << m_timing.species << "  out="
                     << m_timing.output << "  total=" << m_timing.total
                     << " s  |  ETA " << (perStep * (totalSteps - step - 1) / 60.0)
                     << " min" << Alert::Raise;
        std::cout.flush();
      }
    }

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
      Rodin::Examples::Heart::Atrium::Config cfg;

      char buffer[512];
      PetscBool got = PETSC_FALSE;
      PetscOptionsGetString(PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_mesh",
        buffer, sizeof(buffer), &got);
      if (got)
        cfg.meshPath = buffer;

      PetscReal real = 0.0;
      got = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_mesh_scale",
        &real, &got);
      if (got)
        cfg.meshScale = real;

      got = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_dt", &real, &got);
      if (got)
        cfg.dt = real;

      got = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_pin", &real, &got);
      if (got)
        cfg.inletPressureMean = real;

      PetscInt integer = 0;
      got = PETSC_FALSE;
      PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_flow_cycles",
        &integer, &got);
      if (got)
        cfg.flowCycles = static_cast<int>(integer);

      got = PETSC_FALSE;
      PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_species_cycles",
        &integer, &got);
      if (got)
        cfg.speciesCycles = static_cast<int>(integer);

      got = PETSC_FALSE;
      PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_output_every",
        &integer, &got);
      if (got)
        cfg.outputEvery = static_cast<int>(integer);

      got = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_vms_scale",
        &real, &got);
      if (got)
        cfg.vmsScale = real;

      got = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_graddiv_scale",
        &real, &got);
      if (got)
        cfg.gradDivScale = real;

      got = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_pspg_scale",
        &real, &got);
      if (got)
        cfg.pspgScale = real;

      PetscBool flag = PETSC_FALSE;
      got = PETSC_FALSE;
      PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_vms", &flag, &got);
      if (got)
        cfg.useVMS = (flag == PETSC_TRUE);

      got = PETSC_FALSE;
      PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, "-atrium_kinetics",
        &flag, &got);
      if (got)
        cfg.solveKinetics = (flag == PETSC_TRUE);

      Rodin::Examples::Heart::Atrium simulation(context, cfg);
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
