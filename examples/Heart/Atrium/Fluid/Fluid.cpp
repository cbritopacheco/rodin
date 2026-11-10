/*
 *          Copyright Oscar RUZ NUNES 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <iostream>
#include <numbers>
#include <limits>
#include <stdexcept>

#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>

#include <Rodin/Alert.h>
#include <Rodin/Configure.h>

#include "Fluid.h"
#include "FlowIntegral.h"
#include "Output.h"
#include "../Log.h"
#include "../../CoronaryArtery/CoronaryArteryTiming.h"

namespace Rodin::Examples::Heart::Atrium
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

    /**
     * @brief The projections invert mass matrices: CG + Jacobi, never the
     * global direct solver used by the coupled flow system.
     */
    void configureMassSolver(Rodin::Solver::KSP& ksp, const std::string& prefix)
    {
      setPrefixedDefault(prefix, "ksp_type", "cg");
      setPrefixedDefault(prefix, "pc_type", "jacobi");
      ksp.setPrefix(prefix);
    }
  }

  Fluid::Fluid(const ParentMeshType& parent, const Parameters& cfg,
    const RunParameters& run, Attribute material)
    : m_cfg(cfg),
      m_run(run),
      m_input(makeInput(m_cfg, &m_outletPressure)),
      m_model(m_input),
      m_mesh(makeMesh(parent, m_cfg, material)),
      m_output(std::make_unique<Output>(*this, parent.getContext())),
      m_inletSet(m_cfg.inlets.begin(), m_cfg.inlets.end()),
      m_vh(std::integral_constant<size_t, 1>{}, m_mesh, m_mesh.getSpaceDimension()),
      m_sh(std::integral_constant<size_t, 1>{}, m_mesh),
      m_u(m_vh),
      m_p(m_sh),
      m_v(m_vh),
      m_q(m_sh),
      m_uOld(m_vh),
      m_wallVelocity(m_vh),
      m_wallDisplacement(m_vh),
      m_wallTraction(m_vh),
      m_meshVelocity(m_vh),
      m_meshDisplacement(m_vh),
      m_speedMagnitude(m_sh),
      m_sTrial(m_sh),
      m_sTest(m_sh),
      m_wTrial(m_vh),
      m_wTest(m_vh),
      m_tauFn([this](const Point& p) { return vmsTauAt(p); }),
      m_tauCFn([this](const Point& p) { return tauCAt(p); }),
      m_sqrtTauCFn([this](const Point& p) { return sqrtTauCAt(p); }),
      m_tauPFn([this](const Point& p) { return tauPAt(p); }),
      m_viscosityFn([this](const Point& p) { return viscosityAt(p); }),
      m_dragFn([this](const Point& p) { return dragAt(p); }),
      m_piTilde(m_sh),
      m_convProjection(m_vh),
      m_pressureProjection(m_vh),
      m_pT(m_sh),
      m_tT(m_sh),
      m_aT(m_sh),
      m_gT(m_sh),
      m_xT(m_sh),
      m_vP(m_sh),
      m_vT(m_sh),
      m_vA(m_sh),
      m_vG(m_sh),
      m_vX(m_sh),
      m_pCur(m_sh),
      m_tCur(m_sh),
      m_aCur(m_sh),
      m_gCur(m_sh),
      m_xCur(m_sh),
      m_pPrev(m_sh),
      m_tPrev(m_sh),
      m_aPrev(m_sh),
      m_gPrev(m_sh),
      m_xPrev(m_sh),
      m_fibrin(m_sh),
      m_tat(m_sh),
      m_drag(m_sh),
      m_drift(m_sh),
      m_eT(m_sh),
      m_vE(m_sh),
      m_layer(m_sh),
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
      m_transportP(m_pT, m_vP),
      m_transportT(m_tT, m_vT),
      m_transportA(m_aT, m_vA),
      m_transportG(m_gT, m_vG),
      m_transportX(m_xT, m_vX),
      m_kspP(m_transportP),
      m_kspT(m_transportT),
      m_kspA(m_transportA),
      m_kspG(m_transportG),
      m_kspX(m_transportX),
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
      Log("Fluid/Mesh") << "cells=" << cellCount << " vertices=" << vertexCount
                        << " velocity DOFs=" << velocityDOFs
                        << " pressure DOFs=" << pressureDOFs << Alert::Raise;
    }
  }

  Fluid::~Fluid() = default;

  bool Fluid::isRoot() const
  {
    return m_mesh.getContext().getCommunicator().rank() == RootRank;
  }

  Fluid::MeshType Fluid::makeMesh(
    const ParentMeshType& parent, const Parameters& cfg, Attribute material)
  {
    if (!std::isfinite(cfg.meshScale) || cfg.meshScale <= 0)
    {
      Alert::Exception() << "Fluid length conversion must be finite and positive"
                         << Alert::Raise;
    }
    MeshType::Builder builder;
    builder.initialize(parent);
    for (auto cell = parent.getCell(); cell; ++cell)
    {
      if (parent.getShard().isOwned(3, cell->getIndex()) &&
        cell->getAttribute() == material)
        builder.include(3, cell->getIndex());
    }
    auto mesh = builder.finalize();
    if (mesh.getCellCount() == 0)
      Alert::Exception() << "No fluid cells for material " << material << Alert::Raise;
    size_t invalid = 0;
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      if (cell->getGeometry() != Geometry::Polytope::Type::Tetrahedron &&
        mesh.getShard().isOwned(3, cell->getIndex()))
        ++invalid;
    }
    invalid = boost::mpi::all_reduce(
      parent.getContext().getCommunicator(), invalid, std::plus<size_t>());
    if (invalid)
    {
      Alert::Exception() << "The affine-P1 Oseen submesh requires tetrahedra"
                         << Alert::Raise;
    }
    mesh.scale(cfg.meshScale);
    for (size_t d = 1; d <= 3; ++d)
    {
      mesh.getConnectivity().compute(d, 0);
      mesh.getConnectivity().compute(d, d - 1);
      mesh.getConnectivity().compute(d - 1, d);
    }
    return mesh;
  }

  Fluid::Model::Input Fluid::makeInput(const Parameters& cfg, const Real* atrialPressure)
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
    // With feedback the ventricle reads the mean pressure of the mitral plane
    // at the last accepted step; otherwise its own atrial waveform.
    if (cfg.ventricularFeedback)
      input.pAt = [atrialPressure](Real) { return *atrialPressure; };
    else
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

  Fluid::Real Fluid::activationWave(const Activation& cfg, Real t)
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
    {
      return cfg.positiveValue *
        ss((tau - cfg.tRampStart) / (cfg.tRampEnd - cfg.tRampStart));
    }
    if (tau < cfg.tPlateauEnd)
      return cfg.positiveValue;
    if (tau < cfg.tRelaxEnd)
    {
      return cfg.positiveValue +
        (cfg.negativeValue - cfg.positiveValue) *
        ss((tau - cfg.tPlateauEnd) / (cfg.tRelaxEnd - cfg.tPlateauEnd));
    }
    if (tau < cfg.tNegativeEnd)
      return cfg.negativeValue;
    return cfg.negativeValue *
      (1.0 - ss((tau - cfg.tNegativeEnd) / (T - cfg.tNegativeEnd)));
  }

  Fluid::Real Fluid::atrialWave(const AtrialPressure& cfg, Real t)
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

  Fluid::Real Fluid::cellSize(const Point& p)
  {
    return std::pow(p.getPolytope().getMeasure(), 1.0 / p.getPolytope().getDimension());
  }

  const Fluid::CellCoefficients& Fluid::cellCoefficients(const Point& p) const
  {
    const auto& cell = p.getPolytope();
    assert(cell.getDimension() == m_mesh.getDimension());
    if (m_cellCoefficients.empty())
      m_cellCoefficients.resize(
        m_mesh.getShard().getPolytopeCount(m_mesh.getDimension()));
    auto& entry = m_cellCoefficients[cell.getIndex()];
    if (entry.epoch == m_coefficientEpoch)
      return entry;
    const auto& cy = m_cfg.viscosity;

    // Lagged Carreau-Yasuda viscosity of the momentum equation; the strain
    // rate of an affine P1 field is constant on the cell.
    const auto sym = 0.5 * (Jacobian(m_uOld) + Transpose(Jacobian(m_uOld)));
    const Real shear = std::sqrt(
      cy.gammaRegularization * cy.gammaRegularization + 2.0 * Dot(sym, sym).getValue(p));
    entry.viscosity = cy.muInf +
      (cy.mu0 - cy.muInf) *
        std::pow(1.0 + std::pow(cy.lambda * shear, cy.yasuda), (cy.n - 1.0) / cy.yasuda);
    entry.size = cellSize(p);

    // tau_1 = [4 nu/h^2 + 2|u^n - w|/h + sigma_B/rho]^-1 with the relative
    // speed and the drag of the centroid; the Brinkman drag is a reaction and
    // enters the stabilization scale as one.
    const Point center = centroid(p);
    const auto uc = (m_uOld - m_meshVelocity).getValue(center);
    entry.speed = std::sqrt(Math::dot(uc, uc));
    entry.drag = dragAt(center);
    const Real h = entry.size;
    entry.tau1 = 1.0 /
      (4.0 * entry.viscosity / (m_cfg.rho * h * h) + 2.0 * entry.speed / h +
        entry.drag / m_cfg.rho);
    entry.epoch = m_coefficientEpoch;
    return entry;
  }

  Point Fluid::centroid(const Point& p)
  {
    return Point(p.getPolytope(), Math::SpatialPoint{.25, .25, .25});
  }

  Math::SpatialVector<Fluid::Real> Fluid::subscaleAt(const Point& p) const
  {
    if (m_subscale.empty())
      return Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};
    const auto& value = m_subscale[p.getPolytope().getIndex()];
    return Math::SpatialVector<Real>{{value[0], value[1], value[2]}};
  }

  Fluid::Real Fluid::viscosityAt(const Point& p) const
  {
    return cellCoefficients(p).viscosity;
  }

  Fluid::Real Fluid::tau1At(const Point& p) const
  {
    return cellCoefficients(p).tau1;
  }

  Fluid::Real Fluid::vmsTauAt(const Point& p) const
  {
    if (m_cfg.useVMS && m_cfg.quasiStaticSubscales)
      return m_cfg.vmsScale * tau1At(p) / m_cfg.rho;
    return m_cfg.vmsScale / (m_cfg.rho / m_run.timeStep + m_cfg.rho / tau1At(p));
  }

  Fluid::Real Fluid::sqrtTauCAt(const Point& p) const
  {
    const auto& cell = cellCoefficients(p);
    return std::sqrt(
      m_cfg.gradDivScale * m_cfg.rho * cell.size * cell.size / (4.0 * cell.tau1));
  }

  Fluid::Real Fluid::tauCAt(const Point& p) const
  {
    const Real s = sqrtTauCAt(p);
    return s * s;
  }

  Fluid::Real Fluid::tauPAt(const Point& p) const
  {
    if (m_cfg.useVMS)
      return m_cfg.pspgScale * tau1At(p) / m_cfg.rho;
    return m_cfg.pspgScale / (m_cfg.rho / m_run.timeStep + m_cfg.rho / tau1At(p));
  }

  void Fluid::axpy(Real a, const ::Vec& x, ::Vec& y)
  {
    PetscErrorCode ierr = VecAXPY(y, a, x);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGhostUpdateBegin(y, INSERT_VALUES, SCATTER_FORWARD);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGhostUpdateEnd(y, INSERT_VALUES, SCATTER_FORWARD);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
  }

  Fluid& Fluid::initialize()
  {
    if (m_finished || m_failed)
    {
      Alert::Exception() << "Cannot initialize a finished or failed atrium"
                         << Alert::Raise;
    }
    if (m_initialized)
      return *this;
    if (m_cfg.preconditionerLag < 1 || !std::isfinite(m_run.timeStep) ||
      m_run.timeStep <= 0 || !std::isfinite(m_cfg.period) || m_cfg.period <= 0 ||
      m_cfg.period / m_run.timeStep < 0.5 || m_cfg.flowCycles < 0 ||
      !std::isfinite(m_cfg.pressureRamp) || m_cfg.pressureRamp < 0 ||
      !std::isfinite(m_cfg.prescribedOutletPressure) ||
      !std::isfinite(m_cfg.prescribedOutletResistance) ||
      m_cfg.prescribedOutletResistance < 0 ||
      m_cfg.period / m_run.timeStep > std::numeric_limits<int>::max() - 0.5)
    {
      Alert::Exception() << "Invalid atrium time-step or cycle configuration"
                         << Alert::Raise;
    }
    if (m_cfg.useVentricularModel && m_cfg.ventricularFeedback &&
      !(m_cfg.lv.Kat > 0 && m_cfg.lv.Kp > 0 && std::isfinite(1.0 / m_cfg.lv.Kat) &&
        std::isfinite(1.0 / m_cfg.lv.Kp)))
    {
      Alert::Exception() << "The coupled mitral valve requires positive ventricular "
                            "conductances Kat and Kp"
                         << Alert::Raise;
    }
    m_stepsPerCycle = static_cast<int>(m_cfg.period / m_run.timeStep + 0.5);
    m_initialDt = m_run.timeStep;
    m_initialPeriod = m_cfg.period;
    m_initialUseVentricularModel = m_cfg.useVentricularModel;
    m_failed = true;
    setupSpaces();
    setupFlow();
    setupSpecies();
    setupWallShear();
    if (!m_coupled)
      m_output->write(0);

    m_initialized = true;
    m_failed = false;
    return *this;
  }

  void Fluid::setupSpaces()
  {
    // The fluid starts at rest: with a pressure ramp the mitral plane is at
    // zero pressure, which is also the first atrial pressure of the ventricle.
    m_pIn = m_cfg.inletPressureMean;
    m_outletPressure = m_cfg.pressureRamp > 0 ? Real(0) : m_pIn;

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
    m_wallVelocity = zeroVector;
    m_wallDisplacement = zeroVector;
    m_wallTraction = zeroVector;
    m_meshVelocity = zeroVector;
    m_meshDisplacement = zeroVector;
    m_convProjection = zeroVector;
    m_pressureProjection = zeroVector;
    m_subscale.assign(
      m_mesh.getShard().getPolytopeCount(m_mesh.getDimension()), {Real(0), Real(0), Real(0)});
    m_subscaleOld = m_subscale;
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

    const auto& k = m_cfg.kinetics;
    m_pCur = k.prothrombin0;
    m_pPrev = k.prothrombin0;
    m_tCur = k.inletThrombin;
    m_tPrev = k.inletThrombin;
    m_aCur = k.antithrombin0;
    m_aPrev = k.antithrombin0;
    m_gCur = k.fibrinogen0;
    m_gPrev = k.fibrinogen0;
    m_xCur = Real(0);
    m_xPrev = Real(0);
    m_fibrin = Real(0);
    m_tat = Real(0);
    m_drag = Real(0);
    m_drift = Real(0);
    m_layer = Real(0);
    m_wssRec = zeroVector;

    if (m_cfg.useVMS)
      configureMassSolver(m_scalarProjectionKSP, "atrium_sproj_");
    if (m_cfg.solveKinetics)
    {
      // Mass-dominated nonsymmetric transport and a symmetric positive
      // definite layer extension: iterative solvers, independent of the
      // global direct solver used for the flow.
      for (auto* ksp : {&m_kspP, &m_kspT, &m_kspA, &m_kspG, &m_kspX})
      {
        setPrefixedDefault("atrium_ptag_", "ksp_type", "gmres");
        setPrefixedDefault("atrium_ptag_", "pc_type", "bjacobi");
        setPrefixedDefault("atrium_ptag_", "ksp_rtol", "1e-10");
        ksp->setPrefix(std::string("atrium_ptag_"));
      }
      configureMassSolver(m_extensionKSP, "atrium_ptag_layer_");
    }
    configureMassSolver(m_vectorProjectionKSP, "atrium_vproj_");
    configureMassSolver(m_wssKSP, "atrium_wss_");

    m_output->bind();

    m_outletArea = boundaryArea(m_cfg.outlet);

    // Nodal wall indicator from the assembled wall measure: consistent on
    // every rank, unlike a projection restricted to owned boundary faces.
    m_flux = BoundaryIntegral(m_one, m_qFlux).over(m_cfg.wall);
    m_flux.assemble();
    m_wallArea.setData(m_flux.getVector());
    m_wallAreaThreshold = 1.0e-6 * m_wallArea.max();

    m_pOut =
      m_cfg.useVentricularModel ? m_model.getState().pv : m_cfg.prescribedOutletPressure;
    m_mitralZ = (m_cfg.useVentricularModel ? mitralResistance(true)
                                           : m_cfg.prescribedOutletResistance) *
      m_outletArea;

    if (isRoot())
    {
      Log("Fluid/Boundary") << "outlet_area=" << m_outletArea
                            << " m^2 initial_outlet_pressure=" << m_pOut << " Pa"
                            << Alert::Raise;
    }
  }

  Fluid::Real Fluid::boundaryArea(Attribute tag)
  {
    m_flux = BoundaryIntegral(m_one, m_qFlux).over(tag);
    m_flux.assemble();
    return std::max<Real>(m_flux(m_one), 1e-12);
  }

  void Fluid::setupFlow()
  {
    const auto normal = BoundaryNormal(m_mesh);
    const Real rho = m_cfg.rho;
    const Real dt = m_run.timeStep;

    const auto symU = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));
    const auto symV = 0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));

    RealFunction massRatio([this](const Point& p) {
      return m_coupled ? m_massRatios.at(p.getPolytope().getIndex()) : Real(1);
    });
    const auto transport = m_uOld - m_meshVelocity;
    const auto convU = Mult(Jacobian(m_u), transport);
    const auto geometricDivergence = Div(m_uOld) - 2.0 * Div(m_meshVelocity);
    const auto temam = geometricDivergence * Dot(m_u, m_v);
    const auto streamlineTest = rho * Mult(Jacobian(m_v), transport);

    const auto uNormal = Dot(m_u, normal) * normal;
    const auto uTangential = m_u - uNormal;

    // Both pressure openings may reverse flow. With the outward normal,
    // inward relative advection supplies the potentially negative energy flux.
    const auto inletBackflow =
      0.5 * rho * m_cfg.inletBackflowStabilization * Max(-Dot(transport, normal), 0.0);
    const auto outletBackflow =
      0.5 * rho * m_cfg.outletBackflowStabilization * Max(-Dot(transport, normal), 0.0);

    // Time-dependent scalars enter through these, so the form is assigned once
    // and only reassembled.
    RealFunction pInFn = [this](const Point&) { return m_pIn; };
    RealFunction pOutFn = [this](const Point&) { return m_pOut; };
    RealFunction mitralFn = [this](const Point&) { return m_mitralZ; };
    // The affine-tetrahedron volume ratio is constant; P1 history times P1 test is quadratic.
    const auto historyMass = FlowIntegral(massRatio * m_uOld, m_v);

    // Tracked velocity subscale of the previous time level, constant on each cell.
    const VectorFunction subscale(
      m_mesh.getSpaceDimension(), [this](const Point& p) { return subscaleAt(p); });

    // Galerkin terms shared by both stabilizations. The lagged viscosity is
    // read from the cellwise coefficients used by the stabilization scales.
    m_flow = (rho / dt) * FlowIntegral(m_u, m_v) - (rho / dt) * historyMass

      + rho * FlowIntegral(Dot(convU, m_v)) +
      0.5 * rho * FlowIntegral(temam)

      // Streamline and grad-div dissipation along u^n - w.
      + rho * FlowIntegral(m_tauFn * convU, streamlineTest) +
      FlowIntegral(m_tauCFn * Div(m_u), Div(m_v))

      + 2.0 * FlowIntegral(m_viscosityFn * symU, symV) - FlowIntegral(m_p, Div(m_v)) +
      FlowIntegral(Div(m_u), m_q) + m_cfg.pressurePenalty * FlowIntegral(m_p, m_q) +
      FlowIntegral(m_tauPFn * Grad(m_p), Grad(m_q))

      + BoundaryIntegral(pInFn * Dot(m_v, normal)).over(m_inletSet) +
      BoundaryIntegral(pOutFn * Dot(m_v, normal)).over(m_cfg.outlet)

      // Source impedance of the pressure inlets and mitral resistance,
      // both positive semidefinite and assembled implicitly.
      + m_cfg.inletImpedance * BoundaryIntegral(Dot(uNormal, m_v)).over(m_inletSet) +
      BoundaryIntegral(mitralFn * Dot(uNormal, m_v)).over(m_cfg.outlet)

      + m_cfg.inletTangentialDamping *
        BoundaryIntegral(Dot(uTangential, m_v)).over(m_inletSet)

      + BoundaryIntegral(inletBackflow * Dot(m_u, m_v)).over(m_inletSet) +
      BoundaryIntegral(outletBackflow * Dot(m_u, m_v)).over(m_cfg.outlet);
    if (m_cfg.useVMS)
    {
      // Term-by-term orthogonal subscales. The projections of the convection,
      // the divergence, and the pressure gradient belong to the previous
      // time level and are nodal fields carried by the mesh; the dynamic
      // velocity subscale is carried by the cells.
      m_flow -= rho * FlowIntegral(m_tauFn * m_convProjection, streamlineTest);
      if (!m_cfg.quasiStaticSubscales)
        m_flow -= (rho / dt) * FlowIntegral(m_tauFn * subscale, streamlineTest);
      m_flow -= FlowIntegral(m_sqrtTauCFn * m_piTilde, Div(m_v));
      if (m_cfg.projectPressureGradient)
        m_flow -= FlowIntegral(m_tauPFn * m_pressureProjection, Grad(m_q));
    }
    else
    {
      // SUPG tests the remaining affine-P1 momentum residual along the lagged flow.
      m_flow += (rho / dt) * FlowIntegral(m_tauFn * m_u, streamlineTest);
      m_flow -= (rho / dt) * FlowIntegral(m_tauFn * massRatio * m_uOld, streamlineTest);
      m_flow +=
        (0.5 * rho) * FlowIntegral(m_tauFn * geometricDivergence * m_u, streamlineTest);
      m_flow += FlowIntegral(m_tauFn * Grad(m_p), streamlineTest);

      // PSPG tests that same momentum residual with the pressure gradient.
      m_flow += (rho / dt) * FlowIntegral(m_tauPFn * m_u, Grad(m_q));
      m_flow -= (rho / dt) * FlowIntegral(m_tauPFn * massRatio * m_uOld, Grad(m_q));
      m_flow += rho * FlowIntegral(m_tauPFn * convU, Grad(m_q));
      m_flow +=
        (0.5 * rho) * FlowIntegral(m_tauPFn * geometricDivergence * m_u, Grad(m_q));
    }
    if (m_cfg.solveKinetics && m_cfg.brinkman.enabled)
    {
      // Fibrin drag, implicit in the velocity and lagged in the fibrin: a
      // positive semidefinite mass term.
      m_flow += FlowIntegral(m_dragFn * m_u, m_v);
      if (!m_cfg.useVMS)
      {
        // The drag belongs to the momentum residual tested by SUPG and PSPG.
        m_flow += FlowIntegral(m_tauFn * m_dragFn * m_u, streamlineTest);
        m_flow += FlowIntegral(m_tauPFn * m_dragFn * m_u, Grad(m_q));
      }
    }
    if (m_coupled)
    {
      // Robin data live per reference area; integrate on the moving fluid wall.
      RealFunction inverseStretch([this](const Point& point) {
        return m_inverseInterfaceStretches.at(point.getPolytope().getIndex());
      });
      constexpr size_t RobinOrder = 2;
      constexpr size_t InterfaceConvectionOrder = 3;
      auto impedance = BoundaryIntegral(inverseStretch * m_u, m_v);
      auto velocity = BoundaryIntegral(inverseStretch * m_wallVelocity, m_v);
      auto traction = BoundaryIntegral(inverseStretch * m_wallTraction, m_v);
      auto convection = BoundaryIntegral(Dot(transport, normal) * Dot(m_u, m_v));
      impedance.setOrder(RobinOrder);
      velocity.setOrder(RobinOrder);
      traction.setOrder(RobinOrder);
      convection.setOrder(InterfaceConvectionOrder);
      m_flow += m_robinCoefficient * impedance.over(m_cfg.wall);
      m_flow -= m_robinCoefficient * velocity.over(m_cfg.wall);
      m_flow += traction.over(m_cfg.wall);
      m_flow -= (0.5 * rho) * convection.over(m_cfg.wall);
    }
    else
      m_flow += DirichletBC(m_u, m_wallVelocity).on(m_cfg.wall);
  }

  void Fluid::prepareSubscales()
  {
    if (!m_cfg.useVMS)
      return;
    Stage stage("VMS projections", isRoot());
    const auto vmsStart = CoronaryClock::now();
    const Real rho = m_cfg.rho;
    const Real dt = m_run.timeStep;

    // Projections of the previous time level on its own configuration:
    // Pi[(grad u^n)(u^n - w^n)], Pi[sqrt(tau_C) div u^n], and Pi[grad p^n].
    // They remain fixed during the coupling iterations of the step, so the
    // interface tangent carries no derivative of a projection.
    const auto conv = Mult(Jacobian(m_uOld), m_uOld - m_meshVelocity);
    projectVector(conv, m_convProjection);
    if (!m_cfg.quasiStaticSubscales)
    {
      // Dynamic subscale of the convection,
      //   u~^n = tau_s rho [u~^{n-1}/dt - (c^n - Pi c^n)],
      // as its cell mean. A subscale is orthogonal to the finite element
      // space, so it is carried by the cells and never projected onto P1.
      // The convective residual is affine on a cell and its mean is the
      // centroid value.
      for (auto cell = m_mesh.getCell(); cell; ++cell)
      {
        const Index i = cell->getIndex();
        const Point center(*cell, Math::SpatialPoint{.25, .25, .25});
        const auto c = conv.getValue(center);
        const auto proj = m_convProjection.getValue(center);
        const Real tau = vmsTauAt(center);
        for (size_t d = 0; d < 3; ++d)
        {
          m_subscale[i][d] =
            tau * rho * (m_subscaleOld[i][d] / dt - (c(d) - proj(d)));
        }
      }
    }
    project(m_sqrtTauCFn * Div(m_uOld), m_piTilde);
    if (m_cfg.projectPressureGradient)
      projectVector(Grad(m_p.getSolution()), m_pressureProjection);
    m_timing.vms = secondsSince(vmsStart);
  }

  void Fluid::setupWallShear()
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

    const auto traction = VectorFunction(mu * Dot(m_gradRec0, normal),
      mu * Dot(m_gradRec1, normal), mu * Dot(m_gradRec2, normal));
    const auto wallStress = traction - Dot(traction, normal) * normal;

    const Real reg = 1.0e-3;
    m_wssProjection = BoundaryIntegral(Dot(m_wssTrial, m_wssTest)).over(m_cfg.wall) +
      reg * Integral(Dot(m_wssTrial, m_wssTest)) -
      BoundaryIntegral(Dot(wallStress, m_wssTest)).over(m_cfg.wall);
  }

  Fluid::Real Fluid::mitralResistance(bool open) const
  {
    // The coupled valve applies the conductances of the ventricular model,
    // so both sides exchange the same flow.
    if (m_cfg.ventricularFeedback)
      return 1.0 / (open ? m_cfg.lv.Kat : m_cfg.lv.Kp);
    return open ? m_cfg.mitralResistanceOpen : m_cfg.mitralResistanceClosed;
  }

  bool Fluid::advance0D()
  {
    const auto rep = m_model.step(m_run.timeStep);
    m_pOut = m_model.getState().pv;

    if (isRoot())
    {
      const auto& s = m_model.getState();
      Log("Fluid/0D", false, rep.converged ? Log::Status::Success : Log::Status::Failure)
        << (rep.converged ? "converged" : "NOT converged") << "  iter=" << rep.iterations
        << "  |F|=" << rep.finalResidual << "  |  pv=" << s.pv << " par=" << s.par
        << " pd=" << s.pd << " Pa" << Alert::Raise;
    }

    return rep.converged;
  }

  bool Fluid::solveFlow()
  {
    Stage stage("Fluid Oseen solve", isRoot());
    // The lagged velocity or the geometry may have changed since the last solve.
    ++m_coefficientEpoch;

    const auto assemblyStart = CoronaryClock::now();
    {
      Stage assembly("Fluid assembly", isRoot());
      m_flow.assemble();
    }
    m_timing.assembly += secondsSince(assemblyStart);

    if (!m_flowFieldSplitsSet)
    {
      Stage stage("Fluid field splits", isRoot());
      m_flow.setFieldSplits();
      m_flowFieldSplitsSet = true;
    }

    if (isRoot())
    {
      Log("Fluid/KSP") << "solve: Oseen linear system" << Alert::Raise;
      std::cout.flush();
    }

    const auto solveStart = CoronaryClock::now();
    bool refreshed = m_cfg.preconditionerLag == 1 || !m_lastPreconditionerRefreshStep ||
      m_step - *m_lastPreconditionerRefreshStep >= m_cfg.preconditionerLag;
    bool retried = false;
    m_flowKSP.setReusePreconditioner(!refreshed);
    {
      Stage solve("Fluid KSP and preconditioner", isRoot());
      m_flow.solve(m_flowKSP);
    }

    ::KSPConvergedReason reason;
    PetscErrorCode ierr = KSPGetConvergedReason(m_flowKSP.getHandle(), &reason);
    assert(ierr == PETSC_SUCCESS);
    PetscInt iterations = 0;
    ierr = KSPGetIterationNumber(m_flowKSP.getHandle(), &iterations);
    assert(ierr == PETSC_SUCCESS);
    m_linearIterations += static_cast<int>(iterations);

    const char* type = nullptr;
    ierr = KSPGetType(m_flowKSP.getHandle(), &type);
    assert(ierr == PETSC_SUCCESS);
    if (m_cfg.preconditionerLag > 1 && type && std::strcmp(type, KSPPREONLY) == 0)
    {
      Alert::Exception() << "Fluid preconditioner lag greater than one requires "
                            "an iterative KSP, such as GMRES"
                         << Alert::Raise;
    }

    if (reason <= 0 && !refreshed)
    {
      Stage retry("Fluid KSP retry with fresh preconditioner", isRoot());
      // Retry the current linear system before any timestep state is committed.
      ierr = VecSet(m_flow.getLinearSystem().getSolution(), 0.0);
      assert(ierr == PETSC_SUCCESS);
      m_flowKSP.setReusePreconditioner(false);
      m_flow.solve(m_flowKSP);
      ierr = KSPGetConvergedReason(m_flowKSP.getHandle(), &reason);
      assert(ierr == PETSC_SUCCESS);
      ierr = KSPGetIterationNumber(m_flowKSP.getHandle(), &iterations);
      assert(ierr == PETSC_SUCCESS);
      m_linearIterations += static_cast<int>(iterations);
      refreshed = true;
      retried = true;
    }
    m_timing.solve += secondsSince(solveStart);
    m_preconditionerRefreshed |= refreshed;
    m_preconditionerRetried |= retried;
    (void)ierr;

    PetscReal residual = 0;
    ierr = KSPGetResidualNorm(m_flowKSP.getHandle(), &residual);
    assert(ierr == PETSC_SUCCESS);
    if (isRoot())
    {
      Log("Fluid/KSP", false,
        reason <= 0 ? Log::Status::Failure
          : retried ? Log::Status::Retry
                    : Log::Status::Success)
        << "iterations=" << iterations << " residual=" << residual
                       << " reason=" << KSPConvergedReasons[reason]
                       << " pc=" << (refreshed ? "refreshed" : "retained")
                       << " retry=" << (retried ? "yes" : "no")
                       << " elapsed=" << secondsSince(solveStart) << " s (root)"
                       << Alert::Raise;
    }

    if (reason <= 0)
    {
      Alert::Exception() << "Oseen linear solve did not converge, reason "
                         << static_cast<int>(reason) << Alert::Raise;
    }

    if (refreshed)
    {
      m_lastPreconditionerRefreshStep = m_step;
      ++m_preconditionerRefreshes;
    }

    // Observe the trial without committing the Oseen history.
    m_speedMagnitude = Sqrt(Dot(m_u.getSolution(), m_u.getSolution()));
    m_speed = m_speedMagnitude.max();

    const bool accepted = std::isfinite(m_speed) && m_speed <= m_cfg.maxVelocity;
    if (!accepted && isRoot())
    {
      Log("Fluid/Velocity", false, Log::Status::Failure)
        << "max_velocity=" << m_speed << " m/s numerical_guard=" << m_cfg.maxVelocity
        << " m/s" << Alert::Raise;
    }
    return accepted;
  }

  void Fluid::recoverGradient()
  {
    Stage stage("Velocity gradient recovery", isRoot());
    const auto start = CoronaryClock::now();
    const auto& uSol = m_u.getSolution();
    projectVector(VectorFunction(Component(Jacobian(uSol), 0, 0),
                    Component(Jacobian(uSol), 0, 1), Component(Jacobian(uSol), 0, 2)),
      m_gradRec0);
    projectVector(VectorFunction(Component(Jacobian(uSol), 1, 0),
                    Component(Jacobian(uSol), 1, 1), Component(Jacobian(uSol), 1, 2)),
      m_gradRec1);
    projectVector(VectorFunction(Component(Jacobian(uSol), 2, 0),
                    Component(Jacobian(uSol), 2, 1), Component(Jacobian(uSol), 2, 2)),
      m_gradRec2);

    m_timing.shear += secondsSince(start);
  }

  void Fluid::computeWallShear()
  {
    Stage stage("Wall shear recovery", isRoot());
    const auto shearStart = CoronaryClock::now();
    if (isRoot() && m_step < 3)
    {
      Log("Fluid") << "Recovering the wall shear stress ..." << Alert::Raise;
      std::cout.flush();
    }

    // The operator depends on the geometry only.
    if (m_wssEpoch == m_geometryEpoch)
      m_wssProjection.assemble(AssemblyTarget::RHS);
    else
    {
      m_wssProjection.assemble();
      m_wssEpoch = m_geometryEpoch;
    }
    m_wssProjection.solve(m_wssKSP);

    // Keep the wall only: away from it the regularized projection carries a
    // sign-alternating tail that would enter the cycle accumulators. The nodal
    // indicator makes every rank write every local node consistently.
    const size_t sdim = m_mesh.getSpaceDimension();
    m_wssRec.setData(m_wssTrial.getSolution().getData());
    m_wss.project(
      VectorFunction(sdim, [this, sdim](const Point& p) -> Math::SpatialVector<Real> {
        Math::SpatialVector<Real> out(sdim);
        if (onWallNode(p))
          out = m_wssRec.getValue(p);
        else
          out.setZero();
        return out;
      }));

    m_timing.shear += secondsSince(shearStart);
  }

  void Fluid::accumulateWallShear()
  {
    // Cycle accumulators. Two are needed and they are not interchangeable: the
    // vector integral measures how much net direction survives, the scalar one
    // how much shear was applied regardless of direction.
    axpy(m_run.timeStep, m_wss.getData(), m_netShear.getData());

    // Nodal interpolation, NOT an L2 projection. |tau_w| is finite on the wall
    // and zero in the interior; the L2 projection of that onto a continuous P1
    // field undershoots, and a single negative node makes int|tau_w| dt -- and
    // therefore TAWSS, which is non-negative by definition -- come out
    // negative. Interpolating the magnitude at the nodes cannot.
    m_shearMagnitude.project(Sqrt(Dot(m_wss, m_wss)));
    axpy(m_run.timeStep, m_shearMagnitude.getData(), m_absShear.getData());
  }

  void Fluid::closeCycle(Real elapsed)
  {
    if (elapsed <= 0.0)
      return;

    // TAWSS, OSI, and the activation live on the wall only. They are nodal
    // and evaluated with the wall indicator, so partition seams carry the
    // same values on every rank.
    m_tawss.project(RealFunction([this, elapsed](const Point& p) -> Real {
      return onWallNode(p) ? m_absShear.getValue(p) / elapsed : 0.0;
    }));

    // OSI = (1/2)[1 - |int tau_w dt| / int |tau_w| dt]. Both accumulators are
    // built from the same nodal values, so OSI lies in [0, 1/2].
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

    // E = sigma_low(TAWSS) min(2 OSI, 1): a logistic in the endothelial shear
    // switch times the oscillatory index, bounded in [0, 1] so that wall
    // initiation never becomes a sink. Without the oscillatory factor,
    // unidirectional stasis activates the wall as well.
    const Real tauA = m_cfg.kinetics.activationShearStress;
    const Real width = std::max<Real>(m_cfg.kinetics.activationShearWidth, 1e-12);
    const bool useOSI = m_cfg.kinetics.activationUsesOSI;
    m_activation.project(RealFunction([this, tauA, width, useOSI](const Point& p) -> Real {
      if (!onWallNode(p))
        return 0.0;
      const Real low = 1.0 / (1.0 + std::exp((m_tawss.getValue(p) - tauA) / width));
      const Real osi =
        useOSI ? std::clamp<Real>(2.0 * m_osi.getValue(p), 0.0, 1.0) : Real(1);
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

    if (m_cfg.solveKinetics)
      extendActivation();
    m_indicesReady = true;
  }

  void Fluid::computeFluxes()
  {
    Stage stage("Fluid fluxes", isRoot());
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

  void Fluid::step()
  {
    beginStep();
    if (!solveFlow())
    {
      Alert::Exception() << "Fluid trial failed the numerical velocity guard"
                         << Alert::Raise;
    }
    recoverGradient();
    computeWallShear();
    commitStep();
  }

  void Fluid::beginStep()
  {
    if (m_finished || m_failed)
    {
      Alert::Exception() << "Cannot advance a finished or failed atrium" << Alert::Raise;
    }
    if (!m_initialized)
      initialize();
    if (m_run.timeStep != m_initialDt || m_cfg.period != m_initialPeriod ||
      m_cfg.useVentricularModel != m_initialUseVentricularModel)
    {
      Alert::Exception() << "Atrium dt, period, and ventricular mode cannot change "
                            "after initialization"
                         << Alert::Raise;
    }
    if (!m_cfg.useVentricularModel &&
      (!std::isfinite(m_cfg.prescribedOutletPressure) ||
        !std::isfinite(m_cfg.prescribedOutletResistance) ||
        m_cfg.prescribedOutletResistance < 0))
    {
      Alert::Exception() << "Invalid prescribed fluid outlet parameters" << Alert::Raise;
    }
    Stage stage("Fluid history and boundary conditions", isRoot());
    const Real PI = std::numbers::pi_v<Real>;
    // A failed solve can leave partially updated state; it cannot be retried.
    m_failed = true;
    // The species of the previous step may have changed the drag.
    ++m_coefficientEpoch;
    m_timing = Timing{};
    m_stepStart = CoronaryClock::now();
    m_linearIterations = 0;
    m_preconditionerRefreshed = false;
    m_preconditionerRetried = false;

    m_t += m_run.timeStep;

    if (m_cfg.useVentricularModel)
    {
      Stage stage("Ventricular driver", isRoot());
      const auto zeroDStart = CoronaryClock::now();
      const bool advanced = advance0D();
      m_timing.zeroD = secondsSince(zeroDStart);
      if (!advanced)
      {
        Alert::Exception() << "[0D] Newton did not converge at t = " << m_t
                           << Alert::Raise;
      }
    }
    else
    {
      m_pOut = m_cfg.prescribedOutletPressure;
    }

    m_pIn = m_cfg.inletPressureMean +
      m_cfg.inletPressureAmplitude * std::sin(2.0 * PI * m_t / m_cfg.period);
    if (m_t < m_cfg.pressureRamp)
    {
      // A coupled ventricle keeps its own pressure: the valve is closed
      // while the ramped atrium stays below it.
      const Real ramp = m_t / m_cfg.pressureRamp;
      m_pIn *= ramp;
      if (!(m_cfg.useVentricularModel && m_cfg.ventricularFeedback))
        m_pOut *= ramp;
    }

    if (m_cfg.useVentricularModel)
    {
      // Mitral valve, switched on the mean pressure of the mitral plane at
      // the previous step: the test of the ventricular model, p_v <= p_at.
      m_valveOpen = (m_outletPressure > m_pOut);
      m_mitralZ = mitralResistance(m_valveOpen) * m_outletArea;
    }
    else
    {
      m_valveOpen = true;
      m_mitralZ = m_cfg.prescribedOutletResistance * m_outletArea;
    }
    prepareSubscales();
  }

  void Fluid::commitStep()
  {
    Stage stage("Fluid diagnostics and export", isRoot());
    const int step = m_step;
    const int cycle = step / m_stepsPerCycle;
    const bool endOfCycle = step % m_stepsPerCycle == m_stepsPerCycle - 1;
    m_uOld.setData(m_u.getSolution().getData());
    ++m_coefficientEpoch;
    if (m_cfg.useVMS && !m_cfg.quasiStaticSubscales)
      m_subscaleOld = m_subscale;
    accumulateWallShear();
    m_cycleElapsed += m_run.timeStep;

    const auto fluxStart = CoronaryClock::now();
    computeFluxes();
    m_timing.fluxes = secondsSince(fluxStart);

    if (endOfCycle)
    {
      closeCycle(m_cycleElapsed);
      m_cycleElapsed = 0.0;

      // VecMax is collective: every rank must reach it, so the reductions
      // happen outside the root-only print.
      const Real maxTawss = m_tawss.max();
      const Real maxOsi = m_osi.max();
      const Real maxActivation = m_activation.max();

      if (isRoot())
      {
        Log("Fluid/Cycle") << "cycle=" << (cycle + 1) << " maxTAWSS=" << maxTawss
                           << " Pa  maxOSI=" << maxOsi
                           << "  maxActivation=" << maxActivation << Alert::Raise;
      }
    }

    // Species start once the flow cycles are complete and the first cycle
    // indices exist; they are advected by the accepted u - w on the accepted mesh.
    if (m_cfg.solveKinetics && m_indicesReady && cycle >= m_cfg.flowCycles)
      solveSpecies();

    m_output->collectDiagnostics(cycle);

    // The flow is written from the very first step, so the warm-up cycles are
    // available before the species are switched on.
    const bool output =
      endOfCycle || (m_run.outputEvery > 0 && (step + 1) % m_run.outputEvery == 0);
    if (output)
    {
      Stage stage("Fluid field export", isRoot());
      const auto outputStart = CoronaryClock::now();
      m_output->write(m_t);
      m_lastOutputStep = step + 1;
      m_timing.output = secondsSince(outputStart);
    }

    m_timing.total = secondsSince(m_stepStart);

    if (isRoot())
    {
      Log("Fluid/State") << "max_velocity=" << m_speed << " m/s  pIn=" << m_pIn
                         << " Pa  pOut=" << m_pOut << " Pa"
                         << "  qIn=" << m_qIn << "  qOut=" << m_qOut << " m^3/s"
                         << Alert::Raise;

      Log("Fluid/Timing") << "0D=" << m_timing.zeroD << "  vms=" << m_timing.vms
                          << "  asm=" << m_timing.assembly << "  ksp=" << m_timing.solve
                          << "  wss=" << m_timing.shear << "  flux=" << m_timing.fluxes
                          << "  species=" << m_timing.species
                          << "  out=" << m_timing.output
                          << "  elapsed_step=" << m_timing.total << " s (root)"
                          << Alert::Raise;
      std::cout.flush();
    }
    ++m_step;
    m_failed = false;
  }

  void Fluid::finish()
  {
    if (m_finished)
      return;
    if (m_initialized)
    {
      // Preserve the last completed state when callers stop between output events.
      if (!m_failed && m_step > 0 && m_lastOutputStep != m_step)
        m_output->write(m_t);
      m_output->close();
    }
    m_finished = true;
  }

}
