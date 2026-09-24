// LeftAtrium2DPTAG.cpp
//
// Run (from the build directory, after flattening the mesh once with
// examples/Heart/LA2D/make_la2d_mesh.py):
//   mpirun -n 4 ./examples/Heart/LeftAtrium2DPTAG -ptag_dt 1e-2
//
// Options: -ptag_mesh, -ptag_pv, -ptag_mv, -ptag_dt, -ptag_flow_cycles,
//          -ptag_species_cycles, -ptag_output_every, -ptag_vms,
//          -ptag_quasistatic, -ptag_graddiv_tau_lag, -ptag_pspg_residual,
//          -ptag_transport_iterative, -ptag_replay, -ptag_replay_flow_cycles,
//          -ptag_replay_refresh, -ptag_species_substeps, -ptag_kinetics,
//          -ptag_brinkman, -ptag_drag,
//          -ptag_fibrinogen, -ptag_wall_flux, -ptag_layer, -ptag_k_a, -ptag_th_in.
#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <sstream>
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

#include "LeftAtrium2DPTAG.h"

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
    using Clock = std::chrono::steady_clock;

    Real secondsSince(Clock::time_point start)
    {
      return std::chrono::duration<Real>(Clock::now() - start).count();
    }

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

    /// @brief Mass matrices are inverted with CG + Jacobi.
    void configureMassSolver(Rodin::Solver::KSP& ksp, const std::string& prefix)
    {
      setPrefixedDefault(prefix, "ksp_type", "cg");
      setPrefixedDefault(prefix, "pc_type", "jacobi");
      ksp.setPrefix(prefix);
    }
  }

  // ==========================================================================
  // PeriodicWaveform
  // ==========================================================================
  void PeriodicWaveform::load(const std::string& path)
  {
    std::ifstream file(path);
    if (!file)
      throw std::runtime_error("Failed to open the pressure file " + path);

    m_path = path;
    m_t.clear();
    m_p.clear();

    std::string line;
    while (std::getline(file, line))
    {
      const auto first = line.find_first_not_of(" \t\r\n");
      if (first == std::string::npos || line[first] == '#')
        continue;
      std::istringstream row(line.substr(first));
      Real t = 0.0, p = 0.0;
      if (!(row >> t >> p))
        throw std::runtime_error("Malformed line in " + path + ": " + line);
      if (!m_t.empty() && !(t > m_t.back()))
        throw std::runtime_error("Non-increasing time column in " + path);
      m_t.push_back(t);
      m_p.push_back(p);
    }

    if (m_t.size() < 2)
      throw std::runtime_error("Fewer than two samples in " + path);
    m_period = m_t.back() - m_t.front();
    if (!(m_period > 0.0))
      throw std::runtime_error("Zero period in " + path);
    m_min = *std::min_element(m_p.begin(), m_p.end());
    m_max = *std::max_element(m_p.begin(), m_p.end());
  }

  PeriodicWaveform::Real PeriodicWaveform::operator()(Real t) const
  {
    const Real t0 = m_t.front();
    Real tau = t - t0;
    tau -= m_period * std::floor(tau / m_period);
    const Real x = t0 + tau;

    const auto it = std::upper_bound(m_t.begin(), m_t.end(), x);
    if (it == m_t.begin())
      return m_p.front();
    if (it == m_t.end())
      return m_p.back();
    const size_t hi = static_cast<size_t>(it - m_t.begin());
    const size_t lo = hi - 1;
    const Real s = (x - m_t[lo]) / (m_t[hi] - m_t[lo]);
    return (1.0 - s) * m_p[lo] + s * m_p[hi];
  }

  // ==========================================================================
  // LeftAtrium2DPTAG
  // ==========================================================================
  LeftAtrium2DPTAG::MeshType LeftAtrium2DPTAG::makeMesh(
    const Context::MPI& context, const Config& cfg)
  {
    const auto& comm = context.getCommunicator();

    Rodin::MPI::Sharder sharder(context);
    if (comm.rank() == RootRank)
    {
      Geometry::Mesh<Context::Local> mesh;
      mesh.load(cfg.meshPath, IO::FileFormat::MEDIT);
      if (mesh.getSpaceDimension() != 2 || mesh.getDimension() != 2)
        throw std::runtime_error(
          "LeftAtrium2DPTAG expects a planar triangular mesh (MEDIT "
          "\"Dimension 2\"); flatten LA2D_rectLAA.mesh with "
          "examples/Heart/LA2D/make_la2d_mesh.py first.");

      const size_t D = mesh.getDimension();
      mesh.getConnectivity().compute(D, D);
      mesh.getConnectivity().compute(D, 0);
      mesh.getConnectivity().compute(D, D - 1);
      mesh.getConnectivity().compute(D - 1, D);
      mesh.getConnectivity().compute(D - 1, 0);

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
    mesh.reconcile(1);
    return mesh;
  }

  LeftAtrium2DPTAG::LeftAtrium2DPTAG(const Context::MPI& context, const Config& cfg)
    : m_cfg(cfg),
      m_mesh(makeMesh(context, m_cfg)),
      m_xdmf(context.getCommunicator(), m_cfg.xdmfBasename),
      m_inletSet(m_cfg.labels.inlets.begin(), m_cfg.labels.inlets.end()),
      m_wallSet(m_cfg.labels.wall.begin(), m_cfg.labels.wall.end()),
      m_vh(std::integral_constant<size_t, 1>{}, m_mesh, m_mesh.getSpaceDimension()),
      m_sh(std::integral_constant<size_t, 1>{}, m_mesh),
      m_u(m_vh), m_p(m_sh), m_v(m_vh), m_q(m_sh), m_uOld(m_vh), m_uOldOld(m_vh),
      m_sTrial(m_sh), m_sTest(m_sh), m_wTrial(m_vh), m_wTest(m_vh),
      m_tauFn([this](const Point& p) { return vmsTauAt(p); }),
      m_tauCFn([this](const Point& p) { return tauCAt(p); }),
      m_sqrtTauCFn([this](const Point& p) { return sqrtTauCAt(p); }),
      m_sqrtTauCOldFn([this](const Point& p) { return sqrtTauCAt(p, m_uOldOld); }),
      m_tauPFn([this](const Point& p) { return tauPAt(p); }),
      m_dragFn([this](const Point& p) { return dragAt(p); }),
      m_pspgDragFn([this](const Point& p) { return tauPAt(p) * dragAt(p); }),
      m_piTilde(m_sh), m_convProjection(m_vh), m_sub(m_vh), m_subOld(m_vh),
      m_pT(m_sh), m_tT(m_sh), m_aT(m_sh), m_gT(m_sh),
      m_vP(m_sh), m_vT(m_sh), m_vA(m_sh), m_vG(m_sh),
      m_pCur(m_sh), m_tCur(m_sh), m_aCur(m_sh), m_gCur(m_sh),
      m_pPrev(m_sh), m_tPrev(m_sh), m_aPrev(m_sh), m_gPrev(m_sh),
      m_pNext(m_sh), m_tNext(m_sh), m_aNext(m_sh), m_gNext(m_sh),
      m_fibrin(m_sh), m_tat(m_sh), m_drag(m_sh), m_drift(m_sh),
      m_eT(m_sh), m_vE(m_sh), m_layer(m_sh),
      m_wss(m_vh), m_symRec0(m_vh), m_symRec1(m_vh), m_netShear(m_vh),
      m_absShear(m_sh), m_shearMagnitude(m_sh), m_tawss(m_sh), m_osi(m_sh),
      m_activation(m_sh),
      m_qFlux(m_sh), m_one(m_sh), m_flux(m_qFlux),
      m_flow(m_u, m_p, m_v, m_q), m_flowKSP(m_flow),
      m_transportP(m_pT, m_vP), m_transportT(m_tT, m_vT),
      m_transportA(m_aT, m_vA), m_transportG(m_gT, m_vG),
      m_kspP(m_transportP), m_kspT(m_transportT), m_kspA(m_transportA),
      m_kspG(m_transportG),
      m_extension(m_eT, m_vE), m_extensionKSP(m_extension),
      m_scalarProjection(m_sTrial, m_sTest), m_scalarProjectionKSP(m_scalarProjection),
      m_vectorProjection(m_wTrial, m_wTest), m_vectorProjectionKSP(m_vectorProjection),
      m_wssTrial(m_vh), m_wssTest(m_vh),
      m_wssProjection(m_wssTrial, m_wssTest), m_wssKSP(m_wssProjection)
  {
    m_wallSet.insert(m_cfg.labels.appendage.begin(), m_cfg.labels.appendage.end());

    m_inletWave.load(m_cfg.inletPressurePath);
    m_outletWave.load(m_cfg.outletPressurePath);

    if (std::abs(m_inletWave.getPeriod() - m_outletWave.getPeriod()) > 1.0e-9)
      throw std::runtime_error("The inlet and outlet waveforms have different periods.");
    // The cycle indices are accumulated over one forcing period.
    m_cfg.period = m_inletWave.getPeriod();

    const auto cells = m_mesh.getCellCount();
    const auto vertices = m_mesh.getVertexCount();
    if (isRoot())
      Alert::Info() << "[mesh] cells=" << cells << " vertices=" << vertices
                    << "  T=" << m_cfg.period << " s" << Alert::Raise;
  }

  LeftAtrium2DPTAG::~LeftAtrium2DPTAG()
  {
    for (auto& v : m_cycleU)
      if (v)
      {
        PetscErrorCode ierr = VecDestroy(&v);
        assert(ierr == PETSC_SUCCESS);
        (void)ierr;
      }
  }

  bool LeftAtrium2DPTAG::isRoot() const
  {
    return m_mesh.getContext().getCommunicator().rank() == RootRank;
  }

  LeftAtrium2DPTAG::Real LeftAtrium2DPTAG::cellSize(const Point& p)
  {
    return std::pow(p.getPolytope().getMeasure(), 1.0 / p.getPolytope().getDimension());
  }

  void LeftAtrium2DPTAG::axpy(Real a, const ::Vec& x, ::Vec& y)
  {
    PetscErrorCode ierr = VecAXPY(y, a, x);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGhostUpdateBegin(y, INSERT_VALUES, SCATTER_FORWARD);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGhostUpdateEnd(y, INSERT_VALUES, SCATTER_FORWARD);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
  }

  // ---- Pointwise coefficients -----------------------------------------------
  LeftAtrium2DPTAG::Real LeftAtrium2DPTAG::viscosityAt(const Point& p,
    const VectorGridFunctionType& u) const
  {
    const auto& cy = m_cfg.viscosity;
    const auto sym = 0.5 * (Jacobian(u) + Transpose(Jacobian(u)));
    const Real shear = std::sqrt(cy.gammaRegularization * cy.gammaRegularization +
                                 2.0 * Dot(sym, sym).getValue(p));
    return cy.muInf + (cy.mu0 - cy.muInf) *
      std::pow(1.0 + std::pow(cy.lambda * shear, cy.yasuda), (cy.n - 1.0) / cy.yasuda);
  }

  LeftAtrium2DPTAG::Real LeftAtrium2DPTAG::fibrinAt(const Point& p) const
  {
    const Real G0 = m_cfg.kinetics.fibrinogen0;
    return std::clamp<Real>(G0 - m_gCur.getValue(p), 0.0, G0);
  }

  LeftAtrium2DPTAG::Real LeftAtrium2DPTAG::dragAt(const Point& p) const
  {
    return m_cfg.brinkman(fibrinAt(p));
  }

  // tau_1 = [ 4 nu/h^2 + 2|u|/h + sigma_B/rho ]^-1: the Brinkman drag is a
  // reaction and enters the stabilisation scale like one.
  LeftAtrium2DPTAG::Real LeftAtrium2DPTAG::tau1At(const Point& p,
    const VectorGridFunctionType& u) const
  {
    const auto uc = u.getValue(p);
    const Real h = cellSize(p);
    const Real nu = viscosityAt(p, u) / m_cfg.rho;
    return 1.0 / (4.0 * nu / (h * h) + 2.0 * std::sqrt(Math::dot(uc, uc)) / h +
                  dragAt(p) / m_cfg.rho);
  }

  LeftAtrium2DPTAG::Real LeftAtrium2DPTAG::vmsTauAt(const Point& p) const
  {
    if (!m_cfg.useVMS)
      return 0.0;
    if (m_cfg.quasiStaticSubscales)
      return m_cfg.vmsScale * tau1At(p, m_uOld) / m_cfg.rho;
    return m_cfg.vmsScale / (m_cfg.rho / m_cfg.dt + m_cfg.rho / tau1At(p, m_uOld));
  }

  LeftAtrium2DPTAG::Real LeftAtrium2DPTAG::sqrtTauCAt(const Point& p) const
  {
    return sqrtTauCAt(p, m_uOld);
  }

  LeftAtrium2DPTAG::Real LeftAtrium2DPTAG::sqrtTauCAt(const Point& p,
    const VectorGridFunctionType& u) const
  {
    if (!m_cfg.useVMS)
      return 0.0;
    const Real h = cellSize(p);
    return std::sqrt(m_cfg.gradDivScale * m_cfg.rho * h * h / (4.0 * tau1At(p, u)));
  }

  LeftAtrium2DPTAG::Real LeftAtrium2DPTAG::tauCAt(const Point& p) const
  {
    const Real s = sqrtTauCAt(p);
    return s * s;
  }

  LeftAtrium2DPTAG::Real LeftAtrium2DPTAG::tauPAt(const Point& p) const
  {
    return m_cfg.pspgScale * tau1At(p, m_uOld) / m_cfg.rho;
  }

  LeftAtrium2DPTAG::Real LeftAtrium2DPTAG::crosswind(const Point& p, Real cur,
    Real prev, const Math::SpatialVector<Real>& gradient, Real diffusivity) const
  {
    const Real gn = std::sqrt(Math::dot(gradient, gradient));
    if (gn < 1.0e-14)
      return 0.0;
    const auto uc = m_uOld.getValue(p);
    const Real residual = (cur - prev) / m_dtChem + Math::dot(uc, gradient);
    return std::max<Real>(0.0,
      m_cfg.crosswindC * std::abs(residual) * cellSize(p) / (2.0 * gn) - diffusivity);
  }

  // ---- Setup ----------------------------------------------------------------
  LeftAtrium2DPTAG& LeftAtrium2DPTAG::initialize()
  {
    if (m_cfg.speciesSubsteps < 1)
      throw std::runtime_error("speciesSubsteps must be >= 1");
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

  void LeftAtrium2DPTAG::setupSpaces()
  {
    const auto zero = Math::SpatialVector<Real>{{0.0, 0.0}};
    for (auto* g : { &m_uOld, &m_uOldOld, &m_subOld, &m_sub, &m_convProjection,
                     &m_wss, &m_netShear, &m_symRec0, &m_symRec1 })
      *g = zero;
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

    m_outletMeasure = boundaryMeasure(AttributeSet{ m_cfg.labels.outlet });
    m_inletMeasure = boundaryMeasure(m_inletSet);
    m_pIn = m_inletWave(0.0) + m_cfg.pressureOffset;
    m_pOut = m_outletWave(0.0) + m_cfg.pressureOffset;

    // sqrt(2 max|dp|/rho) is the largest velocity the forcing can account for.
    Real maxDp = 0.0;
    for (int i = 0; i <= 2000; ++i)
    {
      const Real t = m_cfg.period * static_cast<Real>(i) / 2000.0;
      maxDp = std::max<Real>(maxDp, std::abs(m_inletWave(t) - m_outletWave(t)));
    }
    m_velocityScale = std::sqrt(2.0 * maxDp / m_cfg.rho);

    const Real Tth = k.inhibitionRate * k.antithrombin0 * k.amplificationHalf /
      (k.amplificationRate * k.prothrombin0 / k.prothrombinRef -
       k.inhibitionRate * k.antithrombin0);
    if (isRoot())
      Alert::Info() << "[scale] max|dp|=" << maxDp << " Pa -> u=" << m_velocityScale
                    << " m/s;  kinetics: k_on A_0=" << (k.inhibitionRate * k.antithrombin0)
                    << " 1/s, activation threshold T_th=" << (Tth * 1e6)
                    << " nM, wall layer sqrt(D/(k_on A_0))="
                    << std::sqrt(k.diffusivityProtein / (k.inhibitionRate * k.antithrombin0))
                    << " m" << Alert::Raise;
  }

  LeftAtrium2DPTAG::Real LeftAtrium2DPTAG::boundaryMeasure(const AttributeSet& tags)
  {
    m_flux = BoundaryIntegral(m_one, m_qFlux).over(tags);
    m_flux.assemble();
    return std::max<Real>(m_flux(m_one), 1e-12);
  }

  void LeftAtrium2DPTAG::setupFlow()
  {
    const size_t dim = m_mesh.getSpaceDimension();
    const auto normal = BoundaryNormal(m_mesh);
    const auto& cy = m_cfg.viscosity;
    const Real deltaMu = cy.mu0 - cy.muInf;
    const Real rho = m_cfg.rho;
    const Real dt = m_cfg.dt;
    const Real pspgR = m_cfg.pspgResidualScale * rho;

    // Carreau-Yasuda viscosity at the lagged velocities u^n and u^{n-1}.
    const auto symLag = 0.5 * (Jacobian(m_uOld) + Transpose(Jacobian(m_uOld)));
    const auto shearLag = Sqrt(cy.gammaRegularization * cy.gammaRegularization +
                               2.0 * Dot(symLag, symLag));
    const auto muLag = cy.muInf +
      deltaMu * Pow(1.0 + Pow(cy.lambda * shearLag, cy.yasuda), (cy.n - 1.0) / cy.yasuda);
    const auto symLagOld = 0.5 * (Jacobian(m_uOldOld) + Transpose(Jacobian(m_uOldOld)));
    const auto shearLagOld = Sqrt(cy.gammaRegularization * cy.gammaRegularization +
                                  2.0 * Dot(symLagOld, symLagOld));
    const auto muLagOld = cy.muInf +
      deltaMu * Pow(1.0 + Pow(cy.lambda * shearLagOld, cy.yasuda), (cy.n - 1.0) / cy.yasuda);

    const auto convU = Mult(Jacobian(m_u), m_uOld);
    const auto temam = Div(m_uOld) * Dot(m_u, m_v);
    const auto uNormal = Dot(m_u, normal) * normal;
    const auto uTangential = m_u - uNormal;

    // Directional do-nothing: + (rho/2) max(-u^n.n, 0) u.v on every pressure
    // boundary cancels the incoming kinetic-energy flux, the only term of the
    // discrete energy balance that is not signed.
    const auto backflow = 0.5 * rho * m_cfg.backflowStabilization *
      Max(-Dot(m_uOld, normal), 0.0);

    RealFunction pInFn = [this](const Point&) { return m_pIn; };
    RealFunction pOutFn = [this](const Point&) { return m_pOut; };

    // Backward Euler, skew-symmetric convection (Temam), VMS orthogonal
    // subscales, IMEX split of the viscous term (implicit mu^{n+1} grad u,
    // explicit sqrt(mu^{n+1} mu^n) grad^T u^n: telescopic, unconditionally
    // stable), PSPG with the consistent residual, and the Brinkman drag
    //   + int sigma_B(F^n) u^{n+1}.v,
    // implicit in u and lagged in F. sigma_B >= 0 makes it a positive
    // semidefinite mass term: it can only remove kinetic energy.
    m_flow = (rho / dt) * Integral(m_u, m_v) - (rho / dt) * Integral(m_uOld, m_v)
           + rho * Integral(Dot(convU, m_v)) + 0.5 * rho * Integral(temam)

           + VMSConvectionBilinearIntegrator(m_u, m_v, m_uOld, m_tauFn, rho)
           - VMSConvectionLinearIntegrator(m_v, m_sub, m_uOld, m_convProjection,
                                           m_tauFn, rho, dt)
           + VMSGradDivBilinearIntegrator(m_u, m_v, m_tauCFn)
           - VMSGradDivLinearIntegrator(m_v, m_piTilde, m_sqrtTauCFn)

           + Integral(muLag * Jacobian(m_u), Jacobian(m_v))
           + Integral(Sqrt(muLag * muLagOld) * Transpose(Jacobian(m_uOld)), Jacobian(m_v))

           + Integral(m_dragFn * m_u, m_v)

           - Integral(m_p, Div(m_v)) + Integral(Div(m_u), m_q)
           + m_cfg.pressurePenalty * Integral(m_p, m_q)

           // PSPG: tau_p grad q . [ rho (u - u^n)/dt + rho (grad u) u^n
           //                        + sigma_B u + grad p ].
           + Integral(m_tauPFn * Grad(m_p), Grad(m_q))
           + (pspgR / dt) * Integral(m_tauPFn * m_u, Grad(m_q))
           - (pspgR / dt) * Integral(m_tauPFn * m_uOld, Grad(m_q))
           + pspgR * Integral(m_tauPFn * convU, Grad(m_q))
           + m_cfg.pspgResidualScale * Integral(m_pspgDragFn * m_u, Grad(m_q))

           + BoundaryIntegral(pInFn * Dot(m_v, normal)).over(m_inletSet)
           + BoundaryIntegral(pOutFn * Dot(m_v, normal)).over(m_cfg.labels.outlet)
           + m_cfg.inletImpedance * BoundaryIntegral(Dot(uNormal, m_v)).over(m_inletSet)
           + m_cfg.inletTangentialDamping *
               BoundaryIntegral(Dot(uTangential, m_v)).over(m_inletSet)
           + BoundaryIntegral(backflow * Dot(m_u, m_v)).over(m_inletSet)
           + BoundaryIntegral(backflow * Dot(m_u, m_v)).over(m_cfg.labels.outlet)

           + DirichletBC(m_u, Zero(dim)).on(m_wallSet);
  }

  void LeftAtrium2DPTAG::assignTransport(ScalarProblemType& problem,
    ScalarTrialFunctionType& c, ScalarTestFunctionType& v,
    ScalarGridFunctionType& cur, ScalarGridFunctionType& prev, Real D, Real inlet,
    Real robinScale, const ScalarGridFunctionType* load)
  {
    const Real dt = m_dtChem;
    const Real scale = m_cfg.kinetics.stabilizationScale;

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

    problem =
        (1.0 / dt) * Integral(c, v) - (1.0 / dt) * Integral(cur, v)
      + D * Integral(Grad(c), Grad(v))
      + Integral(Dot(m_uOld, Grad(c)), v)
      + BoundaryIntegral(robinFn * c, v).over(m_wallSet)
      - BoundaryIntegral(loadFn * v).over(m_wallSet)

      // SUPG on the transient and advective residual.
      + (1.0 / dt) * Integral(tauFn * c, pv) - (1.0 / dt) * Integral(tauFn * cur, pv)
      + Integral(tauFn * Dot(m_uOld, Grad(c)), pv)

      // Codina crosswind: SUPG leaves the crosswind direction free.
      + Integral(kdcFn * Grad(c), Grad(v))
      - Integral(kdcFn * invSpeedSqFn * Dot(m_uOld, Grad(c)), Dot(m_uOld, Grad(v)))

      // Inflow composition on the whole PV patch; the mitral outlet is free.
      + DirichletBC(c, RealFunction(Real(inlet))).on(m_inletSet);
  }

  void LeftAtrium2DPTAG::setupSpecies()
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

  void LeftAtrium2DPTAG::setupWallShear()
  {
    const auto normal = BoundaryNormal(m_mesh);
    const auto& cy = m_cfg.viscosity;
    const auto& uSol = m_u.getSolution();

    const auto symU = 0.5 * (Jacobian(uSol) + Transpose(Jacobian(uSol)));
    const auto shear = Sqrt(cy.gammaRegularization * cy.gammaRegularization +
                            2.0 * Dot(symU, symU));
    const auto mu = cy.muInf + (cy.mu0 - cy.muInf) *
      Pow(1.0 + Pow(cy.lambda * shear, cy.yasuda), (cy.n - 1.0) / cy.yasuda);

    // Traction 2 mu eps(u) n from the recovered rows of 2 eps(u), projected
    // on the wall only (a wall node belongs to two facets with different
    // normals; the projection averages them by facet measure).
    const auto traction = VectorFunction(mu * Dot(m_symRec0, normal),
                                         mu * Dot(m_symRec1, normal));
    const auto wallStress = traction - Dot(traction, normal) * normal;
    const Real reg = 1.0e-3;
    m_wssProjection = BoundaryIntegral(Dot(m_wssTrial, m_wssTest)).over(m_wallSet)
                    + reg * Integral(Dot(m_wssTrial, m_wssTest))
                    - BoundaryIntegral(Dot(wallStress, m_wssTest)).over(m_wallSet);
  }

  void LeftAtrium2DPTAG::setupLayer()
  {
    const Real delta = m_cfg.kinetics.layerThickness;
    if (delta <= 0.0)
      return;
    // Screened Poisson extension of the wall activation: E~ decays into the
    // lumen over the layer thickness; symmetric positive definite, CG.
    m_extension = (delta * delta) * Integral(Grad(m_eT), Grad(m_vE))
                + Integral(m_eT, m_vE)
                + DirichletBC(m_eT, RealFunction([this](const Point& p) -> Real {
                    return m_activation.getValue(p); })).on(m_wallSet);
  }

  void LeftAtrium2DPTAG::extendActivation()
  {
    const Real delta = m_cfg.kinetics.layerThickness;
    if (delta <= 0.0)
      return;
    m_extension.assemble();
    m_extension.solve(m_extensionKSP);
    m_layer.setData(m_eT.getSolution().getData());

    // Dose: int (j_0/delta) scale E~ dV must equal j_0 int_Gw E dA. On a flat
    // wall scale = 1; in corners and pockets narrower than delta the tails
    // overlap and scale < 1.
    m_flux = BoundaryIntegral(m_activation, m_qFlux).over(m_wallSet);
    m_flux.assemble();
    const Real wallDose = m_flux(m_one);
    m_flux = Integral(m_layer, m_qFlux);
    m_flux.assemble();
    const Real layerDose = m_flux(m_one);
    m_layerScale = (layerDose > 0.0) ? delta * wallDose / layerDose : 0.0;

    if (isRoot())
      Alert::Info() << "[layer] int E dA = " << wallDose << " m, int E~ dV = "
                    << layerDose << " m^2, dose scale = " << m_layerScale
                    << Alert::Raise;
  }

  // ---- Time stepping --------------------------------------------------------
  bool LeftAtrium2DPTAG::solveFlow()
  {
    const Real rho = m_cfg.rho;
    const Real dt = m_cfg.dt;

    if (m_cfg.useVMS)
    {
      projectVector(Mult(Jacobian(m_uOld), m_uOld), m_convProjection);
      const size_t dim = m_mesh.getSpaceDimension();
      // Quasi-static: no subscale history, m_sub stays at zero.
      if (!m_cfg.quasiStaticSubscales)
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
      if (m_cfg.lagGradDivTau)
        project(m_sqrtTauCOldFn * Div(m_uOld), m_piTilde);
      else
        project(m_sqrtTauCFn * Div(m_uOld), m_piTilde);
    }

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
    (void)ierr;

    m_uOldOld.setData(m_uOld.getData());
    m_uOld.setData(m_u.getSolution().getData());
    if (m_cfg.useVMS)
      m_subOld.setData(m_sub.getData());

    m_speed = std::max(std::abs(m_uOld.max()), std::abs(m_uOld.min()));
    return reason > 0 && std::isfinite(m_speed) && m_speed <= m_cfg.maxVelocity;
  }

  void LeftAtrium2DPTAG::storeSnapshot(int k)
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

  void LeftAtrium2DPTAG::loadSnapshot(int k)
  {
    const int n = static_cast<int>(m_cycleU.size());
    m_uOldOld.setData(m_cycleU[(k - 1 + n) % n]);
    m_uOld.setData(m_cycleU[k]);
    m_u.getSolution().setData(m_cycleU[k]);
  }

  LeftAtrium2DPTAG::Real LeftAtrium2DPTAG::fibrinMass()
  {
    m_flux = Integral(m_fibrin, m_qFlux);
    m_flux.assemble();
    return m_flux(m_one);
  }

  void LeftAtrium2DPTAG::solveSpecies()
  {
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
  }

  void LeftAtrium2DPTAG::reactionStep()
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

  void LeftAtrium2DPTAG::updateDiagnostics()
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

  void LeftAtrium2DPTAG::computeWallShear()
  {
    const auto& uSol = m_u.getSolution();

    // Rows of 2 eps(u) recovered onto the nodes (grad u_h is cellwise
    // constant on P1).
    const auto jac = Jacobian(uSol);
    const auto offDiagonal = Component(jac, 0, 1) + Component(jac, 1, 0);
    projectVector(VectorFunction(2.0 * Component(jac, 0, 0), offDiagonal), m_symRec0);
    projectVector(VectorFunction(offDiagonal, 2.0 * Component(jac, 1, 1)), m_symRec1);

    m_wssProjection.assemble();
    m_wssProjection.solve(m_wssKSP);

    // Keep the wall only: off it the regularised projection carries a
    // sign-alternating tail that would poison the cycle accumulators.
    const size_t dim = m_mesh.getSpaceDimension();
    const size_t faceDim = m_mesh.getDimension() - 1;
    const auto onWall = [this, faceDim](const Polytope& facet) {
      const auto a = m_mesh.getAttribute(faceDim, facet.getIndex());
      return a && m_wallSet.count(*a);
    };
    const auto& wssSol = m_wssTrial.getSolution();
    m_wss = Math::SpatialVector<Real>{{0.0, 0.0}};
    m_wss.project(Region::Boundary, VectorFunction(dim,
      [&wssSol, dim](const Point& p) -> Math::SpatialVector<Real> {
        const auto w = wssSol.getValue(p);
        Math::SpatialVector<Real> out(dim);
        for (Index c = 0; c < static_cast<Index>(dim); ++c)
          out(c) = w(c);
        return out;
      }), onWall);

    // Cycle accumulators: int tau_w dt (vector) and int |tau_w| dt (nodal
    // magnitude, never an L2 projection, which can go negative).
    axpy(m_cfg.dt, m_wss.getData(), m_netShear.getData());
    m_shearMagnitude.project(Sqrt(Dot(m_wss, m_wss)));
    axpy(m_cfg.dt, m_shearMagnitude.getData(), m_absShear.getData());
  }

  void LeftAtrium2DPTAG::closeCycle(Real elapsed)
  {
    if (elapsed <= 0.0)
      return;

    const size_t faceDim = m_mesh.getDimension() - 1;
    const auto onWall = [this, faceDim](const Polytope& facet) {
      const auto a = m_mesh.getAttribute(faceDim, facet.getIndex());
      return a && m_wallSet.count(*a);
    };

    // TAWSS, OSI and E live on the wall and nowhere else; all three are nodal.
    m_tawss = Real(0);
    m_osi = Real(0);
    m_activation = Real(0);

    m_tawss.project(Region::Boundary, RealFunction([this, elapsed](const Point& p) -> Real {
      return m_absShear.getValue(p) / elapsed; }), onWall);

    m_osi.project(Region::Boundary, RealFunction([this](const Point& p) -> Real {
      const Real abs = m_absShear.getValue(p);
      if (abs <= 0.0)
        return 0.0;
      const auto net = m_netShear.getValue(p);
      return std::clamp<Real>(0.5 * (1.0 - std::sqrt(Math::dot(net, net)) / abs), 0.0, 0.5);
    }), onWall);

    // E = sigma_low(TAWSS) * min(2 OSI, 1): a logistic in the measured
    // endothelial switch (~0.4 Pa) times the oscillatory index. Bounded in
    // [0,1], so the wall exchange can never become a sink of thrombin.
    const Real tauA = m_cfg.kinetics.activationShearStress;
    const Real width = std::max<Real>(m_cfg.kinetics.activationShearWidth, 1e-12);
    m_activation.project(Region::Boundary, RealFunction([this, tauA, width](const Point& p) -> Real {
      const Real low = 1.0 / (1.0 + std::exp((m_tawss.getValue(p) - tauA) / width));
      const Real osi = std::clamp<Real>(2.0 * m_osi.getValue(p), 0.0, 1.0);
      return std::clamp<Real>(low * osi, 0.0, 1.0);
    }), onWall);

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

  void LeftAtrium2DPTAG::computeFluxes()
  {
    const auto normal = BoundaryNormal(m_mesh);
    const auto& uSol = m_u.getSolution();
    // 2D: flow rates per unit depth (m^2/s); n is outward, so qIn < 0 while
    // the veins fill the atrium.
    m_flux = BoundaryIntegral(Dot(uSol, normal), m_qFlux).over(m_inletSet);
    m_flux.assemble();
    m_qIn = m_flux(m_one);
    m_flux = BoundaryIntegral(Dot(uSol, normal), m_qFlux).over(m_cfg.labels.outlet);
    m_flux.assemble();
    m_qOut = m_flux(m_one);
  }

  void LeftAtrium2DPTAG::writeCSVHeader()
  {
    m_csv << "t,cycle,pPV,pMV,qIn,qOut,maxU,maxTAWSS,maxOSI,maxE,"
             "maxT,minP,minA,maxTAT,maxF,maxDrag,invariantDrift\n";
  }

  void LeftAtrium2DPTAG::writeCSVRow(int cycle)
  {
    // Collective reductions on every rank, then a root-only write.
    const Real tawss = m_tawss.max(), osi = m_osi.max(), act = m_activation.max();
    const Real T = m_tCur.max(), P = m_pCur.min(), A = m_aCur.min();
    const Real tat = m_tat.max(), F = m_fibrin.max(), drag = m_drag.max();
    const Real drift = std::max(std::abs(m_drift.max()), std::abs(m_drift.min()));
    if (!isRoot())
      return;
    m_csv << m_t << ',' << cycle << ',' << m_pIn << ',' << m_pOut << ',' << m_qIn
          << ',' << m_qOut << ',' << m_speed << ',' << tawss << ',' << osi << ','
          << act << ',' << T << ',' << P << ',' << A << ',' << tat << ',' << F
          << ',' << drag << ',' << drift << '\n';
    m_csv.flush();
  }

  int LeftAtrium2DPTAG::run()
  {
    if (!m_initialized)
      initialize();

    const int stepsPerCycle = static_cast<int>(m_cfg.period / m_cfg.dt + 0.5);
    const int totalCycles = m_cfg.flowCycles + m_cfg.speciesCycles;
    const int totalSteps = totalCycles * stepsPerCycle;
    const int sub = m_cfg.speciesSubsteps;
    if (stepsPerCycle % sub != 0)
      throw std::runtime_error("speciesSubsteps must divide the steps per cycle");
    if (m_cfg.replay)
      m_cycleU.assign(stepsPerCycle, PETSC_NULLPTR);

    if (isRoot())
      Alert::Info() << "[run] " << totalCycles << " cycles x " << stepsPerCycle
                    << " steps; species from cycle " << (m_cfg.flowCycles + 1)
                    << " (never before the first OSI); Brinkman "
                    << (m_cfg.brinkman.enabled ? "on" : "off") << "; replay "
                    << (m_cfg.replay ? "on" : "off") << Alert::Raise;

    Real cycleElapsed = 0.0;
    Real storedFibrinMass = -1.0;
    int refreshCyclesLeft = 0;
    int flowSteps = 0;
    int chemSteps = 0;
    const auto runStart = Clock::now();

    for (int step = 0; step < totalSteps; ++step)
    {
      m_step = step;
      const auto stepStart = Clock::now();

      m_t += m_cfg.dt;
      const int cycle = step / stepsPerCycle;
      const int k = step % stepsPerCycle;
      const bool endOfCycle = (k == stepsPerCycle - 1);
      const bool warmup = cycle < m_cfg.flowCycles;
      const bool solveFlowNow = !m_cfg.replay || warmup || refreshCyclesLeft > 0;
      // Species advance on every sub-th step with dt_q = sub dt; a replayed
      // step that carries neither flow nor species costs nothing.
      const bool chemStep = ((step + 1) % sub == 0);

      m_pIn = m_inletWave(m_t) + m_cfg.pressureOffset;
      m_pOut = m_outletWave(m_t) + m_cfg.pressureOffset;

      if (solveFlowNow)
      {
        // 1. Flow with the drag of the current fibrin field; wall shear
        //    accumulated and the indices closed once per cycle.
        if (!solveFlow())
        {
          Alert::Exception() << "[flow] diverged at step " << (step + 1)
                             << ": max|u| = " << m_speed << " m/s = "
                             << (m_speed / m_velocityScale) << "x sqrt(2 max|dp|/rho)"
                             << Alert::Raise;
          return 1;
        }
        ++flowSteps;
        computeWallShear();
        computeFluxes();
        cycleElapsed += m_cfg.dt;
        if (m_cfg.replay)
          storeSnapshot(k);
        if (endOfCycle)
        {
          closeCycle(cycleElapsed);
          cycleElapsed = 0.0;
          if (refreshCyclesLeft > 0)
            --refreshCyclesLeft;
        }
      }
      else if (chemStep)
      {
        // 1'. Replay: the stored cycle stands in for the flow; E is frozen.
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
          m_subOld = Math::SpatialVector<Real>{{0.0, 0.0}};
          if (isRoot())
            Alert::Info() << "[replay] cycle " << (cycle + 1) << ": max F = " << fmax
                          << ", fibrin mass = " << mass << " -> refreshing the flow for "
                          << refreshCyclesLeft << " cycles" << Alert::Raise;
        }
      }

      if (endOfCycle || (m_cfg.outputEvery > 0 && step % m_cfg.outputEvery == 0))
        m_xdmf.write(m_t).flush();

      if (isRoot() && (step < 3 || step % 20 == 0) && (solveFlowNow || chemStep))
      {
        const Real elapsed = secondsSince(runStart);
        Alert::Info() << "step " << (step + 1) << "/" << totalSteps << "  t=" << m_t
                      << " s  " << (solveFlowNow ? "flow" : "replay")
                      << "  max|u|=" << m_speed << " m/s  dp=" << (m_pIn - m_pOut)
                      << " Pa  qIn=" << m_qIn << "  step=" << secondsSince(stepStart)
                      << " s  ETA " << (elapsed / (step + 1) * (totalSteps - step - 1) / 60.0)
                      << " min" << Alert::Raise;
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
      using Rodin::Examples::Heart::LeftAtrium2DPTAG;
      LeftAtrium2DPTAG::Config cfg;

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

      getString("-ptag_mesh", cfg.meshPath);
      getString("-ptag_pv", cfg.inletPressurePath);
      getString("-ptag_mv", cfg.outletPressurePath);
      getReal("-ptag_dt", cfg.dt);
      getInt("-ptag_flow_cycles", cfg.flowCycles);
      getInt("-ptag_species_cycles", cfg.speciesCycles);
      getInt("-ptag_output_every", cfg.outputEvery);
      getBool("-ptag_vms", cfg.useVMS);
      getBool("-ptag_quasistatic", cfg.quasiStaticSubscales);
      getBool("-ptag_graddiv_tau_lag", cfg.lagGradDivTau);
      getReal("-ptag_pspg_residual", cfg.pspgResidualScale);
      getBool("-ptag_transport_iterative", cfg.iterativeTransport);
      getBool("-ptag_replay", cfg.replay);
      getInt("-ptag_replay_flow_cycles", cfg.replayFlowCycles);
      getReal("-ptag_replay_refresh", cfg.replayRefresh);
      getInt("-ptag_species_substeps", cfg.speciesSubsteps);
      getBool("-ptag_kinetics", cfg.solveKinetics);
      getBool("-ptag_brinkman", cfg.brinkman.enabled);
      getReal("-ptag_drag", cfg.brinkman.coefficient);
      getReal("-ptag_fibrinogen", cfg.kinetics.fibrinogen0);
      getReal("-ptag_wall_flux", cfg.kinetics.wallFlux);
      getReal("-ptag_layer", cfg.kinetics.layerThickness);
      getReal("-ptag_k_a", cfg.kinetics.amplificationRate);
      getReal("-ptag_th_in", cfg.kinetics.inletThrombin);

      LeftAtrium2DPTAG simulation(context, cfg);
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
