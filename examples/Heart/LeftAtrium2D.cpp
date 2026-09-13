// LeftAtrium2D.cpp
//
// Run (from the build directory):
//   python3 ../examples/Heart/LA2D/make_la2d_mesh.py \
//           ../resources/examples/Heart/LA2D_rectLAA.mesh \
//           ../resources/examples/Heart/LA2D_rectLAA_2D.mesh
//   mpirun -n 4 ./examples/Heart/LeftAtrium2D -la2d_dt 1e-3
//
// Options: -la2d_mesh, -la2d_pv, -la2d_mv, -la2d_mesh_scale, -la2d_dt,
//          -la2d_period, -la2d_flow_cycles, -la2d_species_cycles,
//          -la2d_output_every, -la2d_vms_scale, -la2d_graddiv_scale,
//          -la2d_pspg_scale, -la2d_pspg_residual, -la2d_vms,
//          -la2d_kinetics, -la2d_th_in, -la2d_inlet_impedance.
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
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

#include "LeftAtrium2D.h"
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

  // ==========================================================================
  // PressureWaveform
  // ==========================================================================
  void PressureWaveform::load(const std::string& path)
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
      Real t = 0.0;
      Real p = 0.0;
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

    // Trapezoidal mean over one period, so the reported figure is the mean of
    // the signal and not the mean of the samples: the file is uniformly
    // sampled, but nothing here requires it to be.
    Real integral = 0.0;
    for (size_t i = 0; i + 1 < m_t.size(); ++i)
      integral += 0.5 * (m_p[i] + m_p[i + 1]) * (m_t[i + 1] - m_t[i]);
    m_mean = integral / m_period;
  }

  PressureWaveform::Real PressureWaveform::operator()(Real t) const
  {
    assert(m_t.size() >= 2);

    const Real t0 = m_t.front();
    Real tau = t - t0;
    tau -= m_period * std::floor(tau / m_period);
    const Real x = t0 + tau;

    // upper_bound, then step back: the samples are strictly increasing, so the
    // bracketing interval is [it - 1, it) and the index arithmetic cannot run
    // off either end.
    const auto it = std::upper_bound(m_t.begin(), m_t.end(), x);
    if (it == m_t.begin())
      return m_p.front();
    if (it == m_t.end())
      return m_p.back();

    const size_t hi = static_cast<size_t>(it - m_t.begin());
    const size_t lo = hi - 1;
    const Real dt = m_t[hi] - m_t[lo];
    const Real s = (x - m_t[lo]) / dt;
    return (1.0 - s) * m_p[lo] + s * m_p[hi];
  }

  // ==========================================================================
  // LeftAtrium2D
  // ==========================================================================
  LeftAtrium2D::AttributeSet LeftAtrium2D::makeInletSet(const Config& cfg)
  {
    return AttributeSet(cfg.labels.inlets.begin(), cfg.labels.inlets.end());
  }

  LeftAtrium2D::AttributeSet LeftAtrium2D::makeWallSet(const Config& cfg)
  {
    AttributeSet out(cfg.labels.wall.begin(), cfg.labels.wall.end());
    out.insert(cfg.labels.appendage.begin(), cfg.labels.appendage.end());
    return out;
  }

  LeftAtrium2D::LeftAtrium2D(const Context::MPI& context, const Config& cfg)
    : m_cfg(cfg),
      m_mesh(makeMesh(context, m_cfg)),
      m_xdmf(context.getCommunicator(), m_cfg.xdmfBasename),
      m_inletSet(makeInletSet(m_cfg)),
      m_wallSet(makeWallSet(m_cfg)),
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
      m_wss(m_vh), m_symRec0(m_vh), m_symRec1(m_vh),
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
    m_inletWave.load(m_cfg.inletPressurePath);
    m_outletWave.load(m_cfg.outletPressurePath);

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

      const auto describe = [](const char* what, const PressureWaveform& w) {
        Alert::Info() << "[" << what << "] " << w.getPath() << "  samples="
                      << w.getSampleCount() << "  T=" << w.getPeriod()
                      << " s  mean=" << w.getMean() << " Pa  range=["
                      << w.getMinimum() << ", " << w.getMaximum() << "] Pa"
                      << Alert::Raise;
      };
      describe("p_pv", m_inletWave);
      describe("p_mv", m_outletWave);
    }

    // The two files must share a period, and the run's cycle length must be
    // that period: the wall-shear indices are accumulated over one cycle, and
    // a cycle that is not a period of the forcing averages two different
    // phases of the flow into the same TAWSS.
    const Real tolerance = 1.0e-9;
    if (std::abs(m_inletWave.getPeriod() - m_outletWave.getPeriod()) > tolerance)
      throw std::runtime_error("The inlet and outlet waveforms have different periods.");

    if (std::abs(m_inletWave.getPeriod() - m_cfg.period) > 1.0e-6)
    {
      if (isRoot())
        Alert::Warning() << "[waveform] Config::period = " << m_cfg.period
                         << " s does not match the file period "
                         << m_inletWave.getPeriod()
                         << " s; taking the file's." << Alert::Raise;
      m_cfg.period = m_inletWave.getPeriod();
    }
  }

  LeftAtrium2D::~LeftAtrium2D() = default;

  bool LeftAtrium2D::isRoot() const
  {
    return m_mesh.getContext().getCommunicator().rank() == RootRank;
  }

  LeftAtrium2D::MeshType LeftAtrium2D::makeMesh(
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
          "LeftAtrium2D expects a planar triangular mesh written as MEDIT "
          "\"Dimension 2\". LA2D_rectLAA.mesh is stored as \"Dimension 3\" "
          "with z = 0; flatten it first with "
          "examples/Heart/LA2D/make_la2d_mesh.py.");

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

  LeftAtrium2D::Real LeftAtrium2D::cellSize(const Point& p)
  {
    return std::pow(p.getPolytope().getMeasure(), 1.0 / p.getPolytope().getDimension());
  }

  LeftAtrium2D::Real LeftAtrium2D::viscosityAt(const Point& p) const
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

  LeftAtrium2D::Real LeftAtrium2D::tau1At(const Point& p) const
  {
    const auto uc = m_uOld.getValue(p);
    const Real h = cellSize(p);
    const Real nu = viscosityAt(p) / m_cfg.rho;
    return 1.0 / (4.0 * nu / (h * h) + 2.0 * std::sqrt(Math::dot(uc, uc)) / h);
  }

  LeftAtrium2D::Real LeftAtrium2D::vmsTauAt(const Point& p) const
  {
    // useVMS has to be answered HERE, not only where the projections are
    // computed. The four VMS integrators are part of the form and are
    // assembled on every step; skipping the projections merely freezes
    // Pi[(grad u)u], u' and Pi[sqrt(tau_C) div u] at zero, which turns the
    // orthogonal-subscale method into a plain (non-orthogonal) SUPG/grad-div
    // one instead of removing it. Returning tau = 0 is what actually empties
    // the integrators.
    if (!m_cfg.useVMS)
      return 0.0;
    return m_cfg.vmsScale / (m_cfg.rho / m_cfg.dt + m_cfg.rho / tau1At(p));
  }

  LeftAtrium2D::Real LeftAtrium2D::sqrtTauCAt(const Point& p) const
  {
    if (!m_cfg.useVMS)
      return 0.0;
    const Real h = cellSize(p);
    return std::sqrt(m_cfg.gradDivScale * m_cfg.rho * h * h / (4.0 * tau1At(p)));
  }

  LeftAtrium2D::Real LeftAtrium2D::tauCAt(const Point& p) const
  {
    const Real s = sqrtTauCAt(p);
    return s * s;
  }

  LeftAtrium2D::Real LeftAtrium2D::tauPAt(const Point& p) const
  {
    return m_cfg.pspgScale * tau1At(p) / m_cfg.rho;
  }

  void LeftAtrium2D::axpy(Real a, const ::Vec& x, ::Vec& y)
  {
    PetscErrorCode ierr = VecAXPY(y, a, x);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGhostUpdateBegin(y, INSERT_VALUES, SCATTER_FORWARD);
    assert(ierr == PETSC_SUCCESS);
    ierr = VecGhostUpdateEnd(y, INSERT_VALUES, SCATTER_FORWARD);
    assert(ierr == PETSC_SUCCESS);
    (void)ierr;
  }

  LeftAtrium2D& LeftAtrium2D::initialize()
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

  void LeftAtrium2D::setupSpaces()
  {
    const auto zeroVector = Math::SpatialVector<Real>{{0.0, 0.0}};

    m_uOld = zeroVector;
    m_subOld = zeroVector;
    m_sub = zeroVector;
    m_convProjection = zeroVector;
    m_wss = zeroVector;
    m_netShear = zeroVector;
    m_symRec0 = zeroVector;
    m_symRec1 = zeroVector;

    m_piTilde = Real(0);
    m_absShear = Real(0);
    m_shearMagnitude = Real(0);
    m_tawss = Real(0);
    m_osi = Real(0);
    m_activation = Real(0);
    m_one = Real(1);

    const Real fg0 =
      m_cfg.thrombosis.fibrinogenSinusRhythm / m_cfg.thrombosis.fibrinogenMolarMass;
    m_thCur = Real(m_cfg.inletThrombin);
    m_thPrev = Real(m_cfg.inletThrombin);
    m_fnCur = Real(m_cfg.inletFibrin);
    m_fnPrev = Real(m_cfg.inletFibrin);
    m_fgCur = fg0;
    m_fgPrev = fg0;

    configureMassSolver(m_scalarProjectionKSP, "la2d_sproj_");
    configureMassSolver(m_vectorProjectionKSP, "la2d_vproj_");
    configureMassSolver(m_wssKSP, "la2d_wss_");

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

    // In 2D these are lengths, not areas: the "flux" through a boundary is a
    // volumetric flow per unit depth, m^2/s.
    m_outletMeasure = boundaryMeasure(AttributeSet{ m_cfg.labels.outlet });
    m_inletMeasure = boundaryMeasure(m_inletSet);

    m_pIn = m_inletWave(0.0) + m_cfg.pressureOffset;
    m_pOut = m_outletWave(0.0) + m_cfg.pressureOffset;
    m_outletPressure = m_pOut;

    // The velocity scale the forcing can account for. With rigid walls, no body
    // force and both ends on prescribed pressure, sqrt(2 max|dp| / rho) is the
    // whole budget: a peak far above it has no source in the data, and says the
    // scheme is making energy rather than that the atrium is doing something
    // interesting. Printed next to the divergence guard so the two can be read
    // against each other.
    {
      Real maxDp = 0.0;
      const int samples = 2000;
      for (int i = 0; i <= samples; ++i)
      {
        const Real t = m_cfg.period * static_cast<Real>(i) / samples;
        maxDp = std::max<Real>(maxDp, std::abs(m_inletWave(t) - m_outletWave(t)));
      }
      m_velocityScale = std::sqrt(2.0 * maxDp / m_cfg.rho);

      if (isRoot())
      {
        Alert::Info() << "[scale] max|dp| = " << maxDp
                      << " Pa  ->  sqrt(2 dp/rho) = " << m_velocityScale
                      << " m/s; mean inlet velocity at that head ~ "
                      << (m_velocityScale * m_outletMeasure / m_inletMeasure)
                      << " m/s. Divergence guard at " << m_cfg.maxVelocity
                      << " m/s = " << (m_cfg.maxVelocity / m_velocityScale)
                      << "x the scale." << Alert::Raise;

        if (m_cfg.maxVelocity > 5.0 * m_velocityScale)
          Alert::Warning() << "[scale] the guard sits more than five times "
                              "above the velocity the forcing can account for, "
                              "so a blow-up will run a long way before it trips."
                           << Alert::Raise;
      }
    }

    if (isRoot())
      Alert::Info() << "[boundary] MV length = " << m_outletMeasure
                    << " m  PV length = " << m_inletMeasure
                    << " m  |  p_pv(0) = " << m_pIn << " Pa  p_mv(0) = "
                    << m_pOut << " Pa" << Alert::Raise;
  }

  LeftAtrium2D::Real LeftAtrium2D::boundaryMeasure(const AttributeSet& tags)
  {
    m_flux = BoundaryIntegral(m_one, m_qFlux).over(tags);
    m_flux.assemble();
    return std::max<Real>(m_flux(m_one), 1e-12);
  }

  void LeftAtrium2D::setupFlow()
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

    const auto shearLag = Sqrt(cy.gammaRegularization * cy.gammaRegularization +
                               2.0 * Dot(symLag, symLag));
    const auto muLag = cy.muInf +
      deltaMu * Pow(1.0 + Pow(cy.lambda * shearLag, cy.yasuda),
                    (cy.n - 1.0) / cy.yasuda);

    const auto convU = Mult(Jacobian(m_u), m_uOld);
    const auto temam = Div(m_uOld) * Dot(m_u, m_v);

    const auto uNormal = Dot(m_u, normal) * normal;
    const auto uTangential = m_u - uNormal;

    // Incoming-kinetic-energy stabilisation, on EVERY pressure boundary and on
    // the same branch, max(-u.n, 0).
    //
    // Take v = u in the convective pair above. Integrating by parts,
    //
    //   rho (u^n.grad u, u) + (rho/2)((div u^n) u, u) = (rho/2) int_G (u^n.n)|u|^2,
    //
    // so the discrete energy balance reads
    //
    //   rho/(2 dt) d||u||^2 + 2 mu ||eps||^2
    //       = -(rho/2) int_G (u^n.n)|u|^2 - int_G p_ext (u.n).
    //
    // The normal is OUTWARD -- the run confirms it: with dp > 0 filling the
    // atrium, computeFluxes() reports qIn < 0. Wherever fluid ENTERS, u.n < 0
    // and the first term on the right is a POSITIVE, cubic, unbounded energy
    // source. That is the term that has to be cancelled, and
    // + (rho beta/2) int max(-u^n.n, 0)(u.v) is what cancels it.
    //
    // Reading max(+u^n.n, 0) at the inlets, as Atrium does, arms the branch
    // that is already dissipative and leaves the dangerous one untouched: over
    // the whole filling phase that term is identically zero. Atrium survives it
    // because its inletImpedance = 1e3 Pa s/m supplies a bound of its own;
    // setting that to 0 here, which the 30 Pa driving head demanded, removed
    // the only thing holding the inlet jets down.
    //
    // beta = 1 is the directional do-nothing condition and is what makes the
    // estimate unconditional. It also reinterprets the prescribed p: the
    // boundary now carries p + rho|u_n|^2/2, i.e. a total pressure rather than
    // a static one. Here that is not a detail -- at the target 0.24 m/s the
    // dynamic head is 30 Pa, the same size as the driving head -- so beta is
    // left exposed. See Config.
    const auto inletBackflow = 0.5 * rho * m_cfg.inletBackflowStabilization *
      Max(-Dot(m_uOld, normal), 0.0);
    const auto outletBackflow = 0.5 * rho * m_cfg.outletBackflowStabilization *
      Max(-Dot(m_uOld, normal), 0.0);

    // Time-dependent scalars enter through these, so the form is assigned once
    // and only reassembled.
    RealFunction pInFn = [this](const Point&) { return m_pIn; };
    RealFunction pOutFn = [this](const Point&) { return m_pOut; };

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

           // PSPG, tau_p grad(q) . R_M, required by the equal-order pair.
           // Written out term by term:
           //   grad p          the pressure gradient, always present;
           //   rho (u - u^n)/dt    the transient part;
           //   rho (grad u) u^n    the convective part.
           // tau_p = tau_1/rho, so the last two carry one factor of rho each
           // and the three have the same units. With pspgResidualScale = 0
           // only the first survives and the term is the Brezzi-Pitkaranta
           // penalty of the reference; at 1 the whole group vanishes on the
           // exact solution, which is what makes the method consistent.
           + Integral(m_tauPFn * Grad(m_p), Grad(m_q))
           + (pspgR / dt) * Integral(m_tauPFn * m_u, Grad(m_q))
           - (pspgR / dt) * Integral(m_tauPFn * m_uOld, Grad(m_q))
           + pspgR * Integral(m_tauPFn * convU, Grad(m_q))

           // sigma.n = -p n on both pressure boundaries. The weak form carries
           // -int_G (sigma n).v, so a prescribed p enters with a PLUS sign.
           + BoundaryIntegral(pInFn * Dot(m_v, normal)).over(m_inletSet)
           + BoundaryIntegral(pOutFn * Dot(m_v, normal)).over(m_cfg.labels.outlet)

           // Source impedance of the pressure inlets: positive semidefinite
           // and assembled implicitly. Default 0, see Config.
           + m_cfg.inletImpedance *
               BoundaryIntegral(Dot(uNormal, m_v)).over(m_inletSet)

           + m_cfg.inletTangentialDamping *
               BoundaryIntegral(Dot(uTangential, m_v)).over(m_inletSet)

           + BoundaryIntegral(inletBackflow * Dot(m_u, m_v)).over(m_inletSet)
           + BoundaryIntegral(outletBackflow * Dot(m_u, m_v))
               .over(m_cfg.labels.outlet)

           + DirichletBC(m_u, Zero(dim)).on(m_wallSet);
  }

  LeftAtrium2D::Real LeftAtrium2D::crosswind(const Point& p, Real cur, Real prev,
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

  void LeftAtrium2D::setupSpecies()
  {
    const Real Dth = m_cfg.thrombosis.diffusivityThrombin;
    const Real Dfg = m_cfg.thrombosis.diffusivityFibrinogen;
    const Real Dfn = m_cfg.thrombosis.diffusivityFibrin;
    const Real keff = m_cfg.thrombosis.reactionRate;
    const Real scale = m_cfg.thrombosis.stabilizationScale;
    const Real dt = m_cfg.dt;
    const Real fg0 =
      m_cfg.thrombosis.fibrinogenSinusRhythm / m_cfg.thrombosis.fibrinogenMolarMass;
    const Real thIn = m_cfg.inletThrombin;
    const Real fnIn = m_cfg.inletFibrin;

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
      // The activation field is identically zero until the first cycle has
      // closed, so this term contributes nothing before the OSI exists.
      - m_cfg.thrombosis.thrombinWallFlux *
          BoundaryIntegral(m_activation * m_vth).over(m_wallSet)

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

      // Pure advection needs a Dirichlet condition on the inflow boundary. The
      // veins carry plasma fibrinogen, the circulating thrombin level, and no
      // fibrin. Imposing it on the whole PV patch rather than only where
      // u.n < 0 is what makes this a *boundary* condition and not a switch:
      // the mitral outlet is left free, so the reverse-flow terms there are
      // the ones that have to hold the transport together during backflow.
      + DirichletBC(m_th, RealFunction(Real(thIn))).on(m_inletSet)
      + DirichletBC(m_fn, RealFunction(Real(fnIn))).on(m_inletSet)
      + DirichletBC(m_fg, RealFunction(Real(fg0))).on(m_inletSet);
  }

  void LeftAtrium2D::setupWallShear()
  {
    const auto normal = BoundaryNormal(m_mesh);
    const auto& cy = m_cfg.viscosity;
    const auto& uSol = m_u.getSolution();

    const auto symU = 0.5 * (Jacobian(uSol) + Transpose(Jacobian(uSol)));
    const auto shear = Sqrt(cy.gammaRegularization * cy.gammaRegularization +
                            2.0 * Dot(symU, symU));
    const auto mu = cy.muInf + (cy.mu0 - cy.muInf) *
      Pow(1.0 + Pow(cy.lambda * shear, cy.yasuda), (cy.n - 1.0) / cy.yasuda);

    // t = 2 mu eps(u) n, from the two recovered rows of 2 eps(u).
    const auto traction = VectorFunction(
      mu * Dot(m_symRec0, normal),
      mu * Dot(m_symRec1, normal));
    const auto wallStress = traction - Dot(traction, normal) * normal;

    // An L2 projection restricted to the wall, regularised in the interior so
    // the mass matrix stays invertible off it. It is used here, and not nodal
    // interpolation, for one reason: a wall node belongs to two facets with
    // different normals, so tau_w has two nodal values and the projection is
    // what averages them by facet measure. Everything downstream of this --
    // |tau_w|, TAWSS, OSI, the activation weight -- is nodal.
    const Real reg = 1.0e-3;
    m_wssProjection = BoundaryIntegral(Dot(m_wssTrial, m_wssTest)).over(m_wallSet)
                    + reg * Integral(Dot(m_wssTrial, m_wssTest))
                    - BoundaryIntegral(Dot(wallStress, m_wssTest)).over(m_wallSet);
  }

  bool LeftAtrium2D::solveFlow()
  {
    const Real rho = m_cfg.rho;
    const Real dt = m_cfg.dt;

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
      // square is the implicit coefficient.
      phase("VMS: projecting the grad-div residual");
      project(m_sqrtTauCFn * Div(m_uOld), m_piTilde);
    }

    m_timing.vms = secondsSince(vmsStart);

    if (isRoot() && m_step < 3)
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

    if (isRoot() && m_step < 3)
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

  void LeftAtrium2D::computeWallShear()
  {
    const auto shearStart = CoronaryClock::now();
    const auto& uSol = m_u.getSolution();

    if (isRoot() && m_step < 3)
    {
      ThreeDInfo() << "Recovering the wall shear stress ..." << Alert::Raise;
      std::cout.flush();
    }

    // The two rows of 2 eps(u) = grad u + grad u^T, recovered onto the nodes:
    //   row 0 = (2 du_x/dx,            du_x/dy + du_y/dx)
    //   row 1 = (du_x/dy + du_y/dx,    2 du_y/dy)
    // grad u_h is elementwise constant on P1, so a recovery is unavoidable and
    // an L2 projection is the right one here -- it is a linear functional of
    // the solution, unlike the indices built from it further down.
    const auto jac = Jacobian(uSol);
    const auto offDiagonal = Component(jac, 0, 1) + Component(jac, 1, 0);

    projectVector(
      VectorFunction(2.0 * Component(jac, 0, 0), offDiagonal), m_symRec0);
    projectVector(
      VectorFunction(offDiagonal, 2.0 * Component(jac, 1, 1)), m_symRec1);

    m_wssProjection.assemble();
    m_wssProjection.solve(m_wssKSP);

    // Keep only the wall. tau_w is an L2 projection, not a nodal quantity, and
    // the operator is (M_wall + reg M_vol). ON the wall the boundary mass
    // dominates by reg*h ~ 5e-7, so the recovered traction there is exact --
    // that part is sound. OFF the wall the only equation a node has is
    // reg M x = 0, and for a consistent P1 mass matrix that forces
    // x_i = -(1/M_ii) sum_j M_ij x_j with every M_ij > 0: the tail ALTERNATES
    // IN SIGN from node to node, decaying by about a quarter per layer. A
    // sign-alternating field is exactly what renders as dots rather than as a
    // field, and it is what put a 0/0 inside OSI. It is not a solver failure:
    // after Jacobi the operator has a condition number of about 3.
    //
    // Zeroing it leaves the wall untouched, makes the two cycle accumulators
    // identically zero off the wall, and costs one nodal pass.
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

    // Cycle accumulators. Two are needed and they are not interchangeable: the
    // vector integral measures how much net direction survives, the scalar one
    // how much shear was applied regardless of direction. Both are formed by
    // VecAXPY, i.e. degree of freedom by degree of freedom.
    axpy(m_cfg.dt, m_wss.getData(), m_netShear.getData());

    // Nodal interpolation, NOT an L2 projection. |tau_w| is finite on the wall
    // and zero in the interior; the L2 projection of that onto a continuous P1
    // field undershoots, and a single negative node makes int|tau_w| dt -- and
    // therefore TAWSS, which is non-negative by definition -- come out
    // negative. Interpolating the magnitude at the nodes cannot.
    m_shearMagnitude.project(Sqrt(Dot(m_wss, m_wss)));
    axpy(m_cfg.dt, m_shearMagnitude.getData(), m_absShear.getData());

    m_timing.shear = secondsSince(shearStart);
  }

  void LeftAtrium2D::closeCycle(Real elapsed)
  {
    if (elapsed <= 0.0)
      return;

    // The three indices live ON THE WALL and nowhere else.
    //
    // They used to be projected over the whole domain, which is not merely
    // untidy: off the wall tau_w decays to zero, so TAWSS -> 0, the logistic
    // SATURATES at its maximum 1/(1+exp(-tau_a/w)) = 0.935, and OSI becomes
    // |int tau dt| / int |tau| dt with both accumulators vanishing -- a 0/0
    // that drifts to 1/2. The product is 0.935 * 1 = 0.933, and that is
    // exactly the number the run reported as maxActivation, with maxOSI =
    // 0.499 beside it. The field maximum was sitting in the middle of the
    // cavity, where there is no endothelium at all, and the picture was of
    // the interior rather than of the wall.
    //
    // The wall flux itself was never affected -- BoundaryIntegral(...) over
    // the wall only ever reads wall nodes, whose TAWSS and OSI are genuine --
    // so this changes the diagnostics and the XDMF fields, not the physics of
    // any run already completed. What it does change is the reported maxima,
    // which were interior artefacts and are now wall values.
    const size_t faceDim = m_mesh.getDimension() - 1;
    const auto onWall = [this, faceDim](const Polytope& facet) {
      const auto a = m_mesh.getAttribute(faceDim, facet.getIndex());
      return a && m_wallSet.count(*a);
    };

    // Zeroed first, so a node that is not on the wall carries 0 rather than
    // whatever the previous cycle left there.
    m_tawss = Real(0);
    m_osi = Real(0);
    m_activation = Real(0);

    // TAWSS = (1/T) int |tau_w| dt. Still an exact scaling of the accumulator,
    // node by node: evaluating a P1 field at its own node returns the nodal
    // value, so nothing here can change its sign.
    m_tawss.project(Region::Boundary,
      RealFunction([this, elapsed](const Point& p) -> Real {
        return m_absShear.getValue(p) / elapsed;
      }), onWall);

    // OSI = (1/2)[1 - |int tau_w dt| / int |tau_w| dt]. Both accumulators are
    // built from the same nodal values, so the triangle inequality holds node
    // by node and OSI lands in [0, 1/2] on its own; the clamp is insurance.
    // Evaluated at the nodes: a ratio of two fields taken through quadrature
    // is not the ratio of the two nodal fields, and it is not bounded by 1/2
    // either.
    m_osi.project(Region::Boundary,
      RealFunction([this](const Point& p) -> Real {
        const Real abs = m_absShear.getValue(p);
        if (abs <= 0.0)
          return 0.0;
        const auto net = m_netShear.getValue(p);
        const Real mag = std::sqrt(Math::dot(net, net));
        return std::clamp<Real>(0.5 * (1.0 - mag / abs), 0.0, 0.5);
      }), onWall);

    // Smooth, bounded activation: a logistic in the measured shear threshold
    // times the oscillatory index mapped onto [0,1]. Clamped to [0,1] so the
    // endothelial thrombin flux can never turn into a sink -- a negative
    // activation is what drives thrombin, and with it fibrin, negative.
    // Reads the two fields written just above, so it must come last.
    const Real tauA = m_cfg.thrombosis.activationShearStress;
    const Real width = std::max<Real>(m_cfg.thrombosis.activationShearWidth, 1e-12);
    m_activation.project(Region::Boundary,
      RealFunction([this, tauA, width](const Point& p) -> Real {
        const Real low = 1.0 / (1.0 + std::exp((m_tawss.getValue(p) - tauA) / width));
        const Real osi = std::clamp<Real>(2.0 * m_osi.getValue(p), 0.0, 1.0);
        return std::clamp<Real>(low * osi, 0.0, 1.0);
      }), onWall);

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

    m_indicesReady = true;
  }

  void LeftAtrium2D::solveSpecies()
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

  void LeftAtrium2D::computeFluxes()
  {
    const auto normal = BoundaryNormal(m_mesh);
    const auto& uSol = m_u.getSolution();

    // In 2D these are flow rates per unit depth, m^2/s. n is outward, so qIn
    // is negative while the veins fill the atrium.
    m_flux = BoundaryIntegral(Dot(uSol, normal), m_qFlux).over(m_inletSet);
    m_flux.assemble();
    m_qIn = m_flux(m_one);

    for (size_t i = 0; i < m_cfg.labels.inlets.size(); ++i)
    {
      m_flux =
        BoundaryIntegral(Dot(uSol, normal), m_qFlux).over(m_cfg.labels.inlets[i]);
      m_flux.assemble();
      m_qInPatch[i] = m_flux(m_one);
    }

    m_flux = BoundaryIntegral(Dot(uSol, normal), m_qFlux).over(m_cfg.labels.outlet);
    m_flux.assemble();
    m_qOut = m_flux(m_one);

    m_flux = BoundaryIntegral(m_p.getSolution(), m_qFlux).over(m_cfg.labels.outlet);
    m_flux.assemble();
    m_outletPressure = m_flux(m_one) / m_outletMeasure;
  }

  void LeftAtrium2D::writeCSVHeader()
  {
    m_csv << "t,cycle,pPV,pMV,dp,qIn,qOut,";
    for (const auto tag : m_cfg.labels.inlets)
      m_csv << "qIn" << tag << ',';
    m_csv << "pOutletMean,maxU,uScale,"
          << "maxTAWSS,maxOSI,maxActivation,maxThrombin,minFibrinogen,maxFibrin\n";
  }

  void LeftAtrium2D::writeCSVRow(int cycle)
  {
    // Every one of these reduces over the communicator, so they are taken on
    // all ranks before the root-only write.
    const Real tawss = m_tawss.max();
    const Real osi = m_osi.max();
    const Real activation = m_activation.max();
    const Real th = m_thCur.max();
    const Real fg = m_fgCur.min();
    const Real fn = m_fnCur.max();

    if (!isRoot())
      return;

    m_csv << m_t << ',' << cycle << ',' << m_pIn << ',' << m_pOut << ','
          << (m_pIn - m_pOut) << ',' << m_qIn << ',' << m_qOut << ',';
    for (const auto q : m_qInPatch)
      m_csv << q << ',';
    m_csv << m_outletPressure << ',' << m_speed << ',' << m_velocityScale << ','
          << tawss << ',' << osi << ',' << activation << ',' << th << ','
          << fg << ',' << fn << '\n';
    m_csv.flush();
  }

  int LeftAtrium2D::run()
  {
    if (!m_initialized)
      initialize();

    const int stepsPerCycle = static_cast<int>(m_cfg.period / m_cfg.dt + 0.5);
    const int totalCycles = m_cfg.flowCycles + m_cfg.speciesCycles;
    const int totalSteps = totalCycles * stepsPerCycle;

    if (isRoot())
      Alert::Info() << "[run] " << totalCycles << " cycles of " << stepsPerCycle
                    << " steps (" << totalSteps << " total); species from cycle "
                    << (m_cfg.flowCycles + 1)
                    << " on, and never before the first OSI; XDMF every "
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

      // Both tractions come from the tabulated waveforms. p_pv - p_mv is
      // identically zero while the mitral valve is shut, so the pair carries
      // the valve and no diode is imposed on top of it.
      m_pIn = m_inletWave(m_t) + m_cfg.pressureOffset;
      m_pOut = m_outletWave(m_t) + m_cfg.pressureOffset;

      if (!solveFlow())
      {
        Alert::Exception() << "[flow] diverged at step " << (step + 1)
                           << ": max|u| = " << m_speed << " m/s = "
                           << (m_speed / m_velocityScale)
                           << "x sqrt(2 max|dp|/rho), i.e. a dynamic head of "
                           << (0.5 * m_cfg.rho * m_speed * m_speed)
                           << " Pa against a driving head of at most "
                           << (0.5 * m_cfg.rho * m_velocityScale * m_velocityScale)
                           << " Pa" << Alert::Raise;
        return 1;
      }

      computeWallShear();
      cycleElapsed += m_cfg.dt;

      const auto fluxStart = CoronaryClock::now();
      computeFluxes();
      m_timing.fluxes = secondsSince(fluxStart);

      // The indices are closed BEFORE the species are advanced, so that within
      // this step the activation field the wall flux reads is the one that
      // belongs to the cycle just finished.
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

      // Kinetics: only once the flow is periodic AND the first OSI exists.
      if (m_cfg.solveKinetics && m_indicesReady && cycle >= m_cfg.flowCycles)
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

      if (isRoot() && (m_step < 3 || step % 20 == 0 || endOfCycle))
      {
        const Real elapsed = secondsSince(runStart);
        const Real perStep = elapsed / static_cast<Real>(step + 1);

        Alert::Info() << "---- Step " << (step + 1) << "/" << totalSteps
                      << "  cycle " << (cycle + 1) << "/" << totalCycles
                      << "  t = " << m_t << " s"
                      << (cycle < m_cfg.flowCycles ? "  (warm-up)" : "")
                      << Alert::Raise;

        Alert::Info() << "[2D] max|u|=" << m_speed << " m/s ("
                      << (m_speed / m_velocityScale) << "x scale)  p_pv="
                      << m_pIn << " Pa  p_mv=" << m_pOut << " Pa  dp="
                      << (m_pIn - m_pOut) << " Pa  qIn=" << m_qIn
                      << "  qOut=" << m_qOut << " m^2/s" << Alert::Raise;

        Alert::Info() << "[PV] per ostium (m^2/s): "
                      << m_cfg.labels.inlets[0] << "=" << m_qInPatch[0] << "  "
                      << m_cfg.labels.inlets[1] << "=" << m_qInPatch[1] << "  "
                      << m_cfg.labels.inlets[2] << "=" << m_qInPatch[2] << "  "
                      << m_cfg.labels.inlets[3] << "=" << m_qInPatch[3]
                      << Alert::Raise;

        Alert::Info() << "[timing] vms=" << m_timing.vms
                      << "  asm=" << m_timing.assembly
                      << "  ksp=" << m_timing.solve
                      << "  wss=" << m_timing.shear
                      << "  flux=" << m_timing.fluxes
                      << "  species=" << m_timing.species
                      << "  out=" << m_timing.output
                      << "  total=" << m_timing.total
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
      Rodin::Examples::Heart::LeftAtrium2D::Config cfg;

      char buffer[512];
      PetscBool got = PETSC_FALSE;
      PetscOptionsGetString(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_mesh",
        buffer, sizeof(buffer), &got);
      if (got)
        cfg.meshPath = buffer;

      got = PETSC_FALSE;
      PetscOptionsGetString(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_pv",
        buffer, sizeof(buffer), &got);
      if (got)
        cfg.inletPressurePath = buffer;

      got = PETSC_FALSE;
      PetscOptionsGetString(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_mv",
        buffer, sizeof(buffer), &got);
      if (got)
        cfg.outletPressurePath = buffer;

      PetscReal real = 0.0;
      got = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_mesh_scale",
        &real, &got);
      if (got)
        cfg.meshScale = real;

      got = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_dt", &real, &got);
      if (got)
        cfg.dt = real;

      got = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_period", &real, &got);
      if (got)
        cfg.period = real;

      got = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_th_in", &real, &got);
      if (got)
        cfg.inletThrombin = real;

      got = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_inlet_impedance",
        &real, &got);
      if (got)
        cfg.inletImpedance = real;

      PetscInt integer = 0;
      got = PETSC_FALSE;
      PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_flow_cycles",
        &integer, &got);
      if (got)
        cfg.flowCycles = static_cast<int>(integer);

      got = PETSC_FALSE;
      PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_species_cycles",
        &integer, &got);
      if (got)
        cfg.speciesCycles = static_cast<int>(integer);

      got = PETSC_FALSE;
      PetscOptionsGetInt(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_output_every",
        &integer, &got);
      if (got)
        cfg.outputEvery = static_cast<int>(integer);

      got = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_vms_scale",
        &real, &got);
      if (got)
        cfg.vmsScale = real;

      got = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_graddiv_scale",
        &real, &got);
      if (got)
        cfg.gradDivScale = real;

      got = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_pspg_scale",
        &real, &got);
      if (got)
        cfg.pspgScale = real;

      got = PETSC_FALSE;
      PetscOptionsGetReal(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_pspg_residual",
        &real, &got);
      if (got)
        cfg.pspgResidualScale = real;

      PetscBool flag = PETSC_FALSE;
      got = PETSC_FALSE;
      PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_vms", &flag, &got);
      if (got)
        cfg.useVMS = (flag == PETSC_TRUE);

      got = PETSC_FALSE;
      PetscOptionsGetBool(PETSC_NULLPTR, PETSC_NULLPTR, "-la2d_kinetics",
        &flag, &got);
      if (got)
        cfg.solveKinetics = (flag == PETSC_TRUE);

      Rodin::Examples::Heart::LeftAtrium2D simulation(context, cfg);
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
