/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file CoupledLV0DCoronary3D.cpp
 * @brief Out-of-line implementation of @ref CoupledLV0DCoronary3D.
 *
 * Holds the parts of the coupled 0D left-ventricle / 3D coronary solver that
 * are not needed at the point of instantiation: mesh loading and partitioning,
 * the assembly of the stabilised Navier-Stokes system, the outlet closure
 * (Starling resistor plus tabulated WRMS apparent viscosity) and the time
 * loop with its XDMF/CSV output.
 */
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
#include <Rodin/Solver.h>

#ifdef RODIN_USE_SCOTCH
#include <Rodin/Scotch/MeshPartitioner.h>
#endif

#include <Rodin/Alert/Exception.h>
#include <Rodin/Alert/NewLine.h>
// Rodin/Math/RootFinding/NewtonRaphson.h and Rodin/Math/RungeKutta/RK4.h are no
// longer needed: the WRMS closure is tabulated once (buildWRMSTable) and the
// outlet update is a scalar Newton with a globally positive derivative.
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

    const char* flowModeName(CoupledLV0DCoronary3D::FlowMode mode)
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

        VecSnapshot(const VecSnapshot&) = delete;
        VecSnapshot& operator=(const VecSnapshot&) = delete;

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

  CoupledLV0DCoronary3D::CoupledLV0DCoronary3D(const Context::MPI& context)
    : CoupledLV0DCoronary3D(context, Config{})
  {}

  CoupledLV0DCoronary3D::CoupledLV0DCoronary3D(
    const Context::MPI& context, const Config& cfg)
    : m_cfg(cfg),
      m_input(makeInput(m_cfg)),
      m_model(m_input),
      m_mesh(makeMesh(context, m_cfg)),
      m_xdmf(context.getCommunicator(), m_cfg.xdmfBasename),
      m_uh(std::integral_constant<size_t, 2>{}, m_mesh, m_mesh.getSpaceDimension()),
      m_ph(std::integral_constant<size_t, 1>{}, m_mesh),
      m_uph(std::integral_constant<size_t, 2>{}, m_mesh, m_mesh.getSpaceDimension()),
      m_tauh(m_mesh),
      m_u(m_uh),
      m_p(m_ph),
      m_mu(m_ph),
      m_v(m_uh),
      m_q(m_ph),
      m_r(m_ph),
      m_uOld(m_uh),
      m_pOld(m_ph),
      m_one(m_ph),
      m_qFlux(m_ph),
      m_sub(m_uph),
      m_subOld(m_uph),
      m_up(m_uph),
      m_vp(m_uph),
      m_tau(m_tauh),
      m_t(m_tauh),
      m_tauOld(m_tauh),
      m_flux(m_qFlux),
      m_shearWall(m_uh),
      m_flow(m_u, m_p, m_v, m_q),
      m_flowKSP(m_flow),
      m_flowSolver(m_flowKSP),
      m_viscosityProjection(m_mu, m_r),
      m_viscosityProjectionKSP(m_viscosityProjection),
      m_l2ConvU(m_up, m_vp),
      m_l2ConvUSolver(m_l2ConvU),
      m_subProjection(m_sub, m_vp),
      m_subProjectionSolver(m_subProjection),
      m_tauProjection(m_tau, m_t),
      m_tauProjectionSolver(m_tauProjection),
      m_gradRecTrial(m_uh),
      m_gradRecTest(m_uh),
      m_gradRecProj(m_gradRecTrial, m_gradRecTest),
      m_gradRecKSP(m_gradRecProj),
      m_wssTrial(m_uh),
      m_wssTest(m_uh),
      m_wssProj(m_wssTrial, m_wssTest),
      m_wssKSP(m_wssProj)
  {
    auto cellCount = m_mesh.getCellCount();
    auto vertexCount = m_mesh.getVertexCount();
    auto spaceDim = m_mesh.getSpaceDimension();
    auto velocityDOFs = m_uh.getSize();
    auto pressureDOFs = m_ph.getSize();

    if (isRoot())
    {
      Alert::Info() << "---- Mesh ----" << Alert::NewLine
                    << "Number of elements in mesh: " << cellCount << Alert::NewLine
                    << "Number of vertices in mesh: " << vertexCount << Alert::NewLine
                    << "Mesh space dimension: " << spaceDim << Alert::Raise;

      Alert::Info() << "---- Function spaces ----" << Alert::NewLine
                    << "Velocity space has " << velocityDOFs << " DOFs." << Alert::NewLine
                    << "Pressure space has " << pressureDOFs << " DOFs." << Alert::Raise;
    }
  }

  CoupledLV0DCoronary3D::~CoupledLV0DCoronary3D() = default;

  CoupledLV0DCoronary3D::MeshType CoupledLV0DCoronary3D::makeMesh(
    const Context::MPI& context, const Config& cfg)
  {
    const auto& comm = context.getCommunicator();

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
        std::to_string(mesh.getContext().getCommunicator().rank()) + ".mesh",
      IO::FileFormat::MEDIT);

    return mesh;
  }

  CoupledLV0DCoronary3D::Model::Input CoupledLV0DCoronary3D::makeInput(const Config& cfg)
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
    input.pAt = [p = cfg.atrialPressure](Real t) { return atrial_pressure(p, t); };
    input.u = [a = cfg.activation](Real t) { return periodic_activation(a, t); };
    input.m0 = [low = cfg.lv.relaxationM0Low, high = cfg.lv.relaxationM0High,
                 lowEc = cfg.lv.relaxationM0LowEc,
                 highEc = cfg.lv.relaxationM0HighEc](Real ec) {
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
                  highEc = cfg.lv.relaxationM0HighEc](Real ec) {
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

  CoupledLV0DCoronary3D::WRMSTable CoupledLV0DCoronary3D::buildWRMSTable(
    const Config& cfg, const CarreauYasuda& visc)
  {
    // Universal WRMS apparent viscosity.
    //
    // The Weissenberg-Rabinowitsch-Mooney-Schofield closure gives, for a tube
    // of radius R and length L carrying a generalized-Newtonian fluid,
    //
    //   Q = pi R^3 I(tau_w) / tau_w^3,   I(tau_w) = int_0^{tau_w} tau^2 gd dtau
    //
    // with gd = gd(tau) the inverse of tau = mu(gd) gd. Writing the same flow
    // as Hagen-Poiseuille with an apparent viscosity, Q = pi R^4 dp/(8 mu_ap L),
    // and using tau_w = R dp / 2L, the radius and the length cancel identically:
    //
    //   mu_ap(tau_w) = tau_w^4 / (4 I(tau_w)),   gd_nom = 4Q/(pi R^3) = tau_w/mu_ap
    //
    // so mu_ap is a *universal* function of the wall shear stress for a given
    // rheology. It is the whole constitutive content of the closure, and it is
    // the same curve for every branch and both limbs. Tabulating it once
    // replaces, per residual evaluation, a Newton solve for gd_w and a 100-step
    // RK4 quadrature -- three orders of magnitude of arithmetic -- at a maximum
    // log-log interpolation error of 0.09 per cent over 241 nodes.
    //
    // The Newtonian check is immediate: I = tau_w^4/(4 mu) gives mu_ap = mu.
    const auto& law = cfg.outletFlowLaw;

    const Real mu0 = visc.mu0;
    const Real muInf = visc.muInf;
    const Real lambda = visc.lambda;
    const Real n = visc.n;
    const Real yasuda = visc.yasuda;
    const Real delta = mu0 - muInf;

    // Constitutive law. The WRMS closure is agnostic to which
    // generalized-Newtonian mu(gammadot) is used, so the model is selectable.
    //
    // Quemada is not offered as a better fit but because it separates what
    // Carreau-Yasuda entangles: the haematocrit sets the high-shear level and
    // k_0 the low-shear aggregation rise. A CY pair fitted to a healthy and a
    // hyperviscous condition can share mu_inf, in which case the two fluids
    // differ only below ~50 1/s -- a regime a normally perfused bed never
    // visits, because gamma_rest = r dP/(2 mu_N L) puts every compartment at
    // 10^2-10^3 1/s. Along the haematocrit axis the effect is first order.
    const bool quemada = (cfg.rheologyModel == RheologyModel::Quemada);
    const auto& qp = cfg.quemada;

    // k_0, k_inf and gamma_c are functions of the haematocrit, not constants:
    // the law has a packing limit phi_max = 2/k, so freezing k_0 at its
    // phi = 0.45 value makes the viscosity diverge above phi = 0.46, which is
    // precisely the range a polycythaemia study needs. Cokelet correlations.
    const Real phi = std::clamp<Real>(qp.hematocrit, 0.0, 0.75);
    const Real p2 = phi * phi;
    const Real p3 = p2 * phi;

    const Real qk0 =
      (qp.k0 > 0.0) ? qp.k0 : std::exp(3.874 - 10.41 * phi + 13.8 * p2 - 6.738 * p3);
    const Real qkInf = (qp.kInf > 0.0)
      ? qp.kInf
      : std::exp(1.3435 - 2.803 * phi + 2.711 * p2 - 0.6479 * p3);
    const Real qgc = (qp.gammaC > 0.0)
      ? qp.gammaC
      : std::exp(-6.1508 + 27.923 * phi - 25.6 * p2 + 3.697 * p3);

    auto muQ = [&](Real g) -> Real {
      const Real s = std::sqrt(std::max<Real>(g, 0.0) / qgc);
      const Real k = (qk0 + qkInf * s) / (1.0 + s);
      const Real b = std::max<Real>(1.0 - 0.5 * k * phi, 1e-3);
      return qp.plasmaViscosity / (b * b);
    };

    // Constant-resistance closure: freeze the law at the high-shear plateau of
    // whichever rheology is active. With mu constant the WRMS integral is
    // I = tau_w^4 / (4 mu) exactly, so mu_ap == mu at every table node and the
    // outlet becomes the linear resistance R_inf = 8 mu_inf L / (N pi R^4).
    // The 3D field keeps Config::viscosity: only the reduced closure changes.
    const Real muPlateau = quemada
      ? qp.plasmaViscosity / ((1.0 - 0.5 * qkInf * phi) * (1.0 - 0.5 * qkInf * phi))
      : muInf;

    auto mu = [&](Real g) -> Real {
      if (cfg.constantOutletResistance)
        return muPlateau;
      if (quemada)
        return muQ(g);
      return muInf +
        delta * std::pow(1.0 + std::pow(lambda * g, yasuda), (n - 1.0) / yasuda);
    };

    auto dmu = [&](Real g) -> Real {
      if (cfg.constantOutletResistance)
        return 0.0;
      if (quemada)
      {
        // Central difference in log g: the law is smooth and the table is
        // built once, so an analytic derivative buys nothing here.
        const Real h = 1e-6 * std::max<Real>(g, 1e-12);
        return (muQ(g + h) - muQ(std::max<Real>(g - h, 0.0))) / (2.0 * h);
      }

      const Real base = 1.0 + std::pow(lambda * g, yasuda);

      return delta * (n - 1.0) * std::pow(base, (n - 1.0 - yasuda) / yasuda) *
        std::pow(lambda, yasuda) * std::pow(g, yasuda - 1.0);
    };

    // gd(tau_w): tau = mu(gd) gd is strictly increasing, so a bisection-
    // safeguarded Newton on log gd converges globally. Run once per node.
    auto shearAt = [&](Real tauW) -> Real {
      Real lo = std::log(1e-14);
      Real hi = std::log(std::max<Real>(tauW / muInf, 1e-12)) + 1.0;
      Real s = 0.5 * (lo + hi);

      for (int it = 0; it < law.shearMaxIterations; ++it)
      {
        const Real g = std::exp(s);
        const Real f = mu(g) * g - tauW;

        if (f < 0.0)
          lo = s;
        else
          hi = s;

        const Real df = (mu(g) + g * dmu(g)) * g; // d f / d log g
        Real sNext = (std::abs(df) > 0.0) ? s - f / df : 0.5 * (lo + hi);

        if (!(sNext > lo && sNext < hi))
          sNext = 0.5 * (lo + hi);

        const Real step = std::abs(sNext - s);
        s = sNext;

        if (step < law.shearStepTolerance)
          break;
      }

      return std::exp(s);
    };

    // I(tau_w) = int_0^{gd_w} gd^3 mu^2 (mu + gd mu') dgd, the same integral
    // after the change of variable tau = mu gd. Composite Simpson.
    auto rheologicalIntegral = [&](Real gammaW) -> Real {
      const int m = 2 * (law.integralSteps / 2); // even
      const Real h = gammaW / static_cast<Real>(m);

      auto f = [&](Real g) -> Real {
        if (g <= 0.0)
          return 0.0;
        const Real mg = mu(g);
        return g * g * g * mg * mg * (mg + g * dmu(g));
      };

      Real sum = f(0.0) + f(gammaW);

      for (int i = 1; i < m; ++i)
        sum += ((i % 2 == 1) ? 4.0 : 2.0) * f(static_cast<Real>(i) * h);

      return sum * h / 3.0;
    };

    WRMSTable table;

    const int nodes = std::max(2, law.tableNodes);
    const Real logTauMin = std::log(law.tableTauMin);
    const Real logTauMax = std::log(law.tableTauMax);

    table.logGamma.reserve(static_cast<std::size_t>(nodes));
    table.logMu.reserve(static_cast<std::size_t>(nodes));

    for (int i = 0; i < nodes; ++i)
    {
      const Real tauW = std::exp(logTauMin +
        (logTauMax - logTauMin) * static_cast<Real>(i) / static_cast<Real>(nodes - 1));

      const Real gammaW = shearAt(tauW);
      const Real I = rheologicalIntegral(gammaW);

      if (!(I > 0.0) || !std::isfinite(I))
        continue;

      const Real muAp = std::pow(tauW, 4.0) / (4.0 * I);

      if (!(muAp > 0.0) || !std::isfinite(muAp))
        continue;

      const Real gammaNom = tauW / muAp;

      // Nodes must stay strictly increasing for the lookup to be well posed.
      if (!table.logGamma.empty() && std::log(gammaNom) <= table.logGamma.back())
        continue;

      table.logGamma.push_back(std::log(gammaNom));
      table.logMu.push_back(std::log(muAp));
    }

    if (table.logGamma.size() < 2)
    {
      // Degenerate rheology: fall back to a constant apparent viscosity.
      table.logGamma = {std::log(1e-6), std::log(1e6)};
      table.logMu = {std::log(mu0), std::log(mu0)};
    }

    return table;
  }

  void CoupledLV0DCoronary3D::updateOutlet0D(const Config& cfg, const WRMSTable& wrms,
    const Model& model, RCR& bc, Real Q, Real dt)
  {
    // ---- Coronary outlet: Starling resistor in an intramyocardial bed ------
    //
    //   3D --[R_a, assembled implicitly]--> (p_c, C) --[R_v, throat]--> P_RA
    //                                          ^
    //                                       p_im(t)
    //
    // One state, the microvascular transmural pressure p_tm = p_c - p_im:
    //
    //   C d(p_tm)/dt = q_a - q_v,           q_a = Q  (the measured 3D flux)
    //   q_v = [ p_c - max(p_im, P_RA) ] / (R_v Phi_v)   if the lumen is open
    //   q_v = 0                                          if the throat is shut
    //   p_out = p_im + p_tm            (+ R_a Phi_a Q, assembled in the 3D form)
    //
    // Nothing anywhere clamps the sign of q_a: retrograde epicardial flow in
    // early systole is a *prediction* of this outlet, governed by alpha through
    // p_out, and the outlet backflow stabilization is four orders of magnitude
    // below R_a Phi_a A, so it damps reversal without forbidding it. The
    // implicit resistance term is linear and symmetric and acts identically on
    // both signs -- it is a resistance, not a diode.
    //
    // The guiding statement is that a collapsible vessel embedded in a tissue
    // at pressure p_im drains against p_im, not against its distal reservoir
    // (Permutt-Riley; Downey & Kirk 1975). Everything else follows:
    //
    //  * In systole max(p_im, P_RA) = p_im, so q_v = p_tm/(R_v Phi_v): the
    //    tissue pressure cancels and the drainage is governed by the transmural
    //    pressure alone, which is slow. The previous formulation drained
    //    against a fixed P_RA while the compression pushed p_c up to p_im,
    //    multiplying the venous driving pressure by 7.6 and emptying 1.9 of the
    //    bed's 7.0 mL per beat.
    //  * p_im appears explicitly in p_out. That single term *is* the systolic
    //    impediment: the 3D domain discharges against ~11.2 kPa instead of
    //    ~2.0 kPa, the arteriolar driving drop falls from 7.4 to 0.9 kPa, and
    //    the inflow becomes diastole-dominant with a systolic/diastolic ratio
    //    of 0.30-0.38 over thirty times the range of C. No calibre modulation
    //    is needed to produce it.
    //  * p_tm becomes a stable relaxation variable, bounded below by 0 (the
    //    throat shuts) and above by R_v max(q_a), with time constant C R_v.
    //    Hence V = C p_tm >= 0 by construction: the collapsible-tube law, the
    //    unstressed areas, the parallel multiplicities and the non-negativity
    //    penalty all existed to bound a variable that is now bounded.
    //
    // See the RCR documentation in CoupledLV0DCoronary3D.h.
    const auto& s = model.getState();
    const auto& law = cfg.outletFlowLaw;

    // Operating (runtime) right atrial pressure. May differ from the value the
    // calibration used; see Config::operatingRightAtrialPressure.
    const Real praOp = (cfg.operatingRightAtrialPressure > 0.0)
      ? cfg.operatingRightAtrialPressure
      : cfg.rightAtrialPressure;

    const Real pim = cfg.intramyocardialFraction * s.pv;
    const Real ptmOld = bc.ptm;

    // Drainage pressure of the Starling throat. Outside the myocardium the
    // reservoir is the right atrium; inside, the collapse pressure is p_im.
    const Real pDrain = std::max<Real>(pim, praOp);

    // Is the venular lumen open?
    //
    // A Starling resistor is a check valve *only while it is collapsed*. The
    // gate therefore belongs on the state of the lumen, not on the sign of the
    // driving pressure: clamping q_v whenever p_c < p_drain would forbid
    // retrograde venous flow even with the vessel wide open, which is a real
    // and observed condition (elevated right atrial pressure, tricuspid
    // regurgitation). The distinction:
    //
    //   p_im <= P_RA : no compression anywhere along the segment, so it cannot
    //                  collapse. The limb is an ordinary resistance and q_v is
    //                  free to take either sign.
    //   p_im  > P_RA : the downstream end is compressed. The throat is a check
    //                  valve, and it is shut once p_tm <= 0.
    //
    // Backflow through the *open* branch is self-limiting and needs no floor:
    // q_v < 0 gives C dp_tm/dt > 0, which raises p_tm and closes the gradient.
    // The Jacobian stays C/dt + G_v > 0 in both branches.
    const bool waterfall = pim > praOp;

    auto venousLumenOpen = [&](Real p) { return !waterfall || p > 0.0; };

    const Real Rv = std::max<Real>(bc.Rv, 1e-300);
    const Real C = std::max<Real>(bc.C, 1e-300);
    const Real q0 = std::max<Real>(std::abs(bc.q0), 1e-300);

    // Rheological modulation. The apparent viscosity is the universal WRMS
    // curve evaluated at the nominal shear rate of the limb, which scales with
    // the flow through the calibrated operating point.
    //
    // The normalization is the *fixed* Newtonian calibration viscosity, never
    // the running rheology. That distinction is the whole point: R_a and R_v
    // are the Newtonian resistances of the pressure budget at mu_N, so a change
    // of blood properties moves Phi and therefore the flow. Normalizing by the
    // running rheology instead would make the calibration absorb the change and
    // leave the mean flow identical between a healthy and a hyperviscous run --
    // which is precisely the insensitivity the previous formulation showed.
    const Real muN = std::max<Real>(cfg.newtonianCalibrationViscosity, 1e-300);

    auto viscosityFactorV = [&](Real q) -> Real {
      const Real aq = std::abs(q);
      const Real g =
        (aq < law.zeroFlowTolerance) ? law.zeroFlowTolerance : bc.gammaV * aq / q0;
      return wrms(g) / muN;
    };

    // Implicit Euler, scalar Newton on p_tm. dq_v/dp_tm >= 0 everywhere, so
    // R' = C/dt + dq_v/dp_tm > 0 and the iteration converges globally from any
    // iterate: no line search, no bracketing, no nested solve.
    //
    // Phi_v is a function of |q_v| rather than of p_tm, so it is held fixed
    // inside the Newton step and refreshed between steps -- a fixed point on a
    // factor that moves by a few per cent per time step, converging with the
    // Newton iteration itself. It is seeded from the previous step's drainage,
    // not from zero: at zero flow the Carreau-Yasuda plateau gives mu_0, which
    // is eighty times mu_inf and would make the first step wildly off.
    Real ptm = ptmOld;
    Real qv = bc.qd;
    Real phiV = viscosityFactorV(qv);
    bool converged = false;

    // The venous lumen is regularised rather than switched. A hard shut leaves
    // J = C/dt, i.e. a pure integrator of the 3D flux Q with no restoring term,
    // so a single transient backflow drives p_tm arbitrarily negative and the
    // outlet never recovers. Keeping a small leak conductance bounds J away
    // from that degenerate case and makes the switch differentiable enough for
    // Newton to stop chattering across p_tm = 0.
    const Real leak = std::max<Real>(law.closedLumenLeak, 0.0);

    for (int it = 0; it < law.outletMaxIterations; ++it)
    {
      phiV = viscosityFactorV(qv);

      const Real drive = (ptm + pim) - pDrain;
      const Real Gv = 1.0 / (Rv * phiV);
      const Real Ge = venousLumenOpen(ptm) ? Gv : leak * Gv;

      qv = drive * Ge;

      const Real R = C * (ptm - ptmOld) / dt - Q + qv;
      const Real J = C / dt + Ge;

      const Real d = -R / J;
      ptm += d;

      // The bound the formulation assumes, imposed rather than hoped for.
      if (law.clampTransmuralPressure)
        ptm = std::max<Real>(ptm, 0.0);

      if (std::abs(d) < law.outletStepTolerance * (1.0 + std::abs(ptm)))
      {
        converged = true;
        break;
      }
    }

    if (!converged || !std::isfinite(ptm))
    {
      // Silently reusing the previous state lets a local failure propagate: the
      // next step inherits an iterate inconsistent with its own flux. Count it
      // and surface it, so a run that limps is visible instead of merely odd.
      static int failures = 0;
      ++failures;
      if (failures <= 20)
        std::cerr << "Warning: coronary outlet solve did not converge "
                  << "(failure " << failures
                  << (failures == 20 ? ", further ones suppressed" : "")
                  << "). Keeping previous state.\n";
      ptm = ptmOld;
    }

    if (law.clampTransmuralPressure)
      ptm = std::max<Real>(ptm, 0.0);

    // Final consistent evaluation at the converged state.
    const Real pc = ptm + pim;
    const Real drive = pc - pDrain;

    phiV = viscosityFactorV(qv);
    qv = (venousLumenOpen(ptm) ? 1.0 : leak) * drive / (Rv * phiV);

    bc.ptm = ptm;
    bc.pim = pim;
    bc.pc = pc;
    bc.qd = qv;
    bc.vol = C * std::max<Real>(ptm, 0.0);
    bc.muV = phiV;

    // Arteriolar modulation. R_a Phi_a is assembled on the 3D outlet boundary,
    // so the factor is exported here and consumed by the next 3D solve.
    const Real aQ = std::abs(Q);
    const Real gammaA =
      (aQ < law.zeroFlowTolerance) ? law.zeroFlowTolerance : bc.gammaA * aQ / q0;
    bc.muA = wrms(gammaA) / muN;

    // Pressure applied to the 3D outlet as a Neumann traction. The resistive
    // part, R_a Phi_a Q, is *not* included here: it is assembled implicitly as
    // R_a Phi_a A (u.n)(v.n), which is what removes the one-step lag on the
    // dominant resistance and, with it, both the need for an outlet capacitor
    // and the L_3D-C_a resonance it produced (f_0 ~ 51 Hz, zeta = 0.25).
    bc.pout = pc;
  }

  CoupledLV0DCoronary3D::Real CoupledLV0DCoronary3D::periodic_activation(
    const Activation& cfg, Real t)
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

  // Smooth return from the negative plateau back to baseline over the rest of
  // the cycle (removes the old instantaneous negativeValue -> 0 jump).
    return cfg.negativeValue *
      (1.0 - ss((tau - cfg.tNegativeEnd) / (T - cfg.tNegativeEnd)));
  }

  CoupledLV0DCoronary3D::Real CoupledLV0DCoronary3D::atrial_pressure(
    const AtrialPressure& cfg, Real t)
  {
    const Real T = cfg.period;
    const Real tau = t - T * std::floor(t / T);

  // C1-continuous (smoothstep) ramps between the same plateau values, so the
  // prescribed atrial/venous pressure has no derivative kinks.
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

  CoupledLV0DCoronary3D& CoupledLV0DCoronary3D::initialize()
  {
    setupModel();
    setupMeshAndSpaces();
    setupProjectionSolvers();
    setupDiagnostics();
    printInitialState();

    m_initialized = true;
    return *this;
  }

  void CoupledLV0DCoronary3D::setupProjectionSolvers()
  {
    // The output-path projections all invert symmetric positive-definite mass
    // matrices, for which Jacobi-preconditioned CG converges in a handful of
    // iterations with no matrix factorization at all. Give each its own PETSc
    // options prefix and default it to CG + Jacobi, so these solves no longer
    // inherit the global direct (MUMPS LU) solver used for the coupled flow
    // system. The defaults are only applied when the user has not set the
    // corresponding option, so they remain overridable from the command line.
    const auto setPrefixedDefault = [](const std::string& prefix, const char* suffix,
                                      const char* value) {
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
    };

    const auto configureMassSolver = [&](Rodin::Solver::KSP& ksp,
                                       const std::string& prefix) {
      setPrefixedDefault(prefix, "ksp_type", "cg");
      setPrefixedDefault(prefix, "pc_type", "jacobi");
      ksp.setPrefix(prefix);
    };

    configureMassSolver(m_viscosityProjectionKSP, "coronary_visc_proj_");
    configureMassSolver(m_gradRecKSP, "coronary_gradrec_proj_");
    configureMassSolver(m_wssKSP, "coronary_wss_proj_");
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
      Alert::Info() << "Setting up " << m_cfg.xdmfBasename << ".xdmf ..." << Alert::Raise;

    m_xdmf.setMesh(m_mesh);

    m_u.setName("u");
    m_p.setName("p");
    m_mu.setName("viscosity");
    m_up.setName("projected_convection");
    m_sub.setName("subscale");
    m_tau.setName("tau");
    m_shearWall.setName("shearStress");

    m_uOld = Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};
    m_pOld = 0.0;
    m_one = 1.0;
    m_mu.getSolution() = 0.0;
    m_subOld = Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};
    m_tauOld = 0.0;
    m_shearWall = Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};

    m_xdmf.add("velocity", m_u.getSolution());
    m_xdmf.add("pressure", m_p.getSolution());
    m_xdmf.add("viscosity", m_mu.getSolution());
    m_xdmf.add("subscale", m_sub.getSolution());
    m_xdmf.add("tau", m_tau.getSolution());
    m_xdmf.add("shearStress", m_shearWall);

    m_wk.clear();
    for (const Attribute tag : m_cfg.outlets)
      m_wk.emplace(tag, m_cfg.defaultRCR);

    // Universal WRMS apparent-viscosity table. Built once, shared by every
    // outlet and both limbs: mu_ap depends only on the wall shear stress, not
    // on radius, length or branch.
    m_wrms = buildWRMSTable(m_cfg, m_cfg.viscosity);

  // ---- Outlet calibration -------------------------------------------------
  //
  // Three divisions per outlet. The model needs exactly the three lumped
  // quantities that appear in the balance -- R_a, R_v and C -- so the
  // calibration produces those and nothing else. The previous scheme built an
  // explicit microvascular geometry (radii, effective lengths, parallel
  // multiplicities, unstressed areas, wall stiffnesses, transit times) from
  // which those three were emergent; that geometry never entered the physics
  // except through them, and it carried nine constants and a nested solve.
  //
  // Per branch, with dP = p_ar(0) - P_RA the resting budget and the Murray
  // split w_i = r_i^3 / sum_j r_j^3 measured on the mesh:
  //
  //   Q_i = Q_tot w_i,      C_i = C_tot w_i,
  //   R_v,i = f_v dP / Q_i, R_a,i = (1 - f_v) dP / Q_i,
  //   p_tm,i(0) = P_RA + f_v dP - alpha p_LV(0).
  //
  // The epicardial conduit is not a separate element any more: its share of
  // the head loss is under 1 per cent (its own budget table says so) and it is
  // absorbed into R_a. The characteristic impedance and the epicardial
  // compliance are gone with it -- Z_c existed only to damp C_a, and C_a
  // existed only to keep the lagged coupling stable, which the implicit
  // assembly of R_a now does unconditionally.
  //
  // The run-off time constant tau = C_i R_v,i is a *prediction*: it is what
  // the diastolic decay of the CSV should be checked against.
    if (m_cfg.autoCalibrateOutlets)
    {
      const Real PI = std::numbers::pi_v<Real>;

      // Driving pressure available at rest. In steady state the storage term
      // vanishes and the network drains to the right atrium: the tissue
      // pressure is the reference of the compartment, not a source, so it does
      // not appear in the steady pressure path.
      const Real par0 = m_model.getState().par;
      const Real dP = std::max<Real>(par0 - m_cfg.rightAtrialPressure, 1.0);

      // Anatomical outlet areas and Murray-law flow split, r = sqrt(A/pi).
      std::map<Attribute, Real> areaOf;
      std::map<Attribute, Real> rEq;
      Real sumR3 = 0.0;

      for (const Attribute tag : m_cfg.outlets)
      {
        m_flux = BoundaryIntegral(m_one, m_qFlux).over(tag);
        m_flux.assemble();
        const Real area = std::max<Real>(m_flux(m_one), 1e-12);
        areaOf[tag] = area;
        rEq[tag] = std::sqrt(area / PI);
        sumR3 += rEq[tag] * rEq[tag] * rEq[tag];
      }

      const Real fv = m_cfg.venularPressureFraction;
      const Real dPv = std::max<Real>(fv * dP, 1.0);
      const Real dPa = std::max<Real>((1.0 - fv) * dP, 1.0);

      const Real pimRest = m_cfg.intramyocardialFraction * m_model.getState().pv;

      const Real muN = std::max<Real>(m_cfg.newtonianCalibrationViscosity, 1e-300);

      // ---- Anatomical rheological operating point -------------------------
      //
      // The prescribed pair is (r, v): a calibre and a red-cell velocity, which
      // is what intravital microscopy measures. Everything else is derived,
      //
      //   g_0 = 4 v / r                     (Poiseuille wall shear rate)
      //   N   = Q_i / (pi r^2 v)            (bed multiplicity)
      //   L   = r dP_share / (2 mu_N g_0)   (pressure budget)
      //   T   = L / v                       (mean transit time)
      //
      // so the effective path length is an output. That is the right way
      // round: L is a lumped path through several generations in series and is
      // not a measurable quantity, whereas r and v are. The shear rate is where
      // the closure reads the constitutive law, so it decides how much of the
      // rheology the model can see at all -- over the plausible morphometric
      // range it moves the predicted effect of a rheology change from 4 to 27
      // per cent -- and it must therefore rest on measured quantities.
      //
      // Three measurables (r, v, T) for two degrees of freedom leaves one
      // consistency check, reported below: the derived transit time must match
      // the indicator-dilution value.
      //
      // A corollary worth stating: with a physiological pressure budget the
      // identity g_0 = r dP/(2 mu_N L) puts *every* compartment -- terminal
      // arteriole, capillary, post-capillary venule, collecting vein -- between
      // 10^2 and 10^3 1/s. The erythrocyte-aggregation regime below 50 1/s does
      // not exist at rest anywhere in the bed; it is reached only in low-flow
      // states. That is why a Carreau-Yasuda pair sharing mu_inf cannot express
      // itself here, and why the haematocrit axis (Quemada) is the one that can.
      const Real ra = std::max<Real>(m_cfg.arteriolarRadius, 1e-12);
      const Real va = std::max<Real>(m_cfg.arteriolarVelocity, 1e-12);
      const Real rv = std::max<Real>(m_cfg.venularRadius, 1e-12);
      const Real vv = std::max<Real>(m_cfg.venularVelocity, 1e-12);

      const Real gammaA0 = 4.0 * va / ra;
      const Real gammaV0 = 4.0 * vv / rv;

      const Real La = ra * dPa / (2.0 * muN * gammaA0);
      const Real Lv = rv * dPv / (2.0 * muN * gammaV0);

      const Real Ta = La / va;
      const Real Tv = Lv / vv;

      // Tone of the distal bed.  fVaso = 1 is the resting calibration; a
      // hyperaemic study (adenosine, exercise) is fVaso > 1, which divides the
      // *bed* resistance while leaving the epicardial 3D domain untouched.
      // This is the only knob that makes an epicardial stenosis flow-limiting:
      // at rest the bed carries ~99.7 per cent of the total resistance, so a
      // geometric lesion upstream is invisible; once R_a, R_v are divided by
      // 4-5 the lesion becomes the dominant term and the model traverses the
      // Gould curve.
      const Real fVaso = std::max<Real>(m_cfg.vasodilationFactor, 1e-6);

      // What dilation does to the operating shear rate.
      //
      // Dividing R by f at fixed N is, through R = 8 mu L/(pi r^4 N), a radius
      // change r -> r f^{1/4}.  The flow rises by f, so the velocity rises by
      // f/f^{1/2} = f^{1/2}, and the wall shear rate
      //
      //     g = 4 v / r  ->  g_0 f^{1/2} / f^{1/4} = g_0 f^{1/4}
      //
      // rises by only f^{1/4}.  This is the quantitative content of the
      // shear-regulation picture: a fourfold hyperaemia moves the shear rate by
      // 41 per cent, i.e. a quarter of a decade on a curve whose knee is three
      // decades away, so the *rheology* is essentially frozen through the
      // vasodilator response and the flow reserve is a resistance effect, not a
      // viscosity effect.  At f = 1 every expression below reduces exactly to
      // the resting calibration, bit for bit.
      const Real fShear = std::pow(fVaso, 0.25);
      const Real gammaA = gammaA0 * fShear;
      const Real gammaV = gammaV0 * fShear;

      // R_a and R_v remain the Newtonian resistances of the budget, so a change
      // of blood properties moves Phi and therefore the flow. The
      // normalization is mu_N and never the running rheology.
      const Real phiA0 = m_wrms(gammaA) / muN;
      const Real phiV0 = m_wrms(gammaV) / muN;

      // Initial condition: steady state of the *actual* (non-Newtonian)
      // network, so the run does not open with a spurious transient.  The
      // venular drop is R_v Phi_v Q = (dPv/(Q_i f)) Phi_v (Q_i f) = dPv Phi_v,
      // independent of f, so this line is unchanged by dilation.
      const Real ptmRest = m_cfg.rightAtrialPressure + dPv * phiV0 - pimRest;

      Real volumeTotal = 0.0;
      Real flowTotal = 0.0;
      Real countTotal = 0.0;

      for (const Attribute tag : m_cfg.outlets)
      {
        const Real w = (rEq[tag] * rEq[tag] * rEq[tag]) / std::max(sumR3, 1e-30);

        // Structural design flow of this territory.
        //
        // Two ways of fixing it, algebraically identical at the default
        // settings but epistemically different:
        //
        //   (a) prescribed  Qi = lcaTargetFlow * w.  The total bed flow is an
        //       input; the arteriolar count N_a = Qi/(pi r_a^2 v_a) is then a
        //       derived number that the model reports but nothing checks.
        //
        //   (b) morphometric  Qi = pi r_a^2 v_a * (arteriolarCount * w).  The
        //       bed is described by how many terminal arterioles it has, and
        //       the flow becomes a *prediction* that can be compared against
        //       the 1.5 mL/s of the literature.  Substituting the pressure
        //       budget dPa = the arteriolar share of the drop, the resulting
        //       Ra = dPa/Qi is exactly the Poiseuille resistance of N_a
        //       parallel tubes,  Ra = 8 mu_N L_a / (pi r_a^4 N_a),  with L_a
        //       the length computed above.  So (b) is not a different closure,
        //       it is the same closure read in the opposite direction.
        //
        // arteriolarCount is set so that (b) reproduces (a) bit for bit at the
        // default radius and velocity; the flag exists so that a morphometric
        // hypothesis (rarefaction in diabetes, capillary drop-out) can be
        // stated directly as a change of N_a instead of a change of the total
        // flow.
        const Real QiStruct = m_cfg.morphometricResistance
          ? PI * ra * ra * va * std::max<Real>(m_cfg.arteriolarCount * w, 0.0)
          : m_cfg.lcaTargetFlow * w;
        const Real Qi = std::max<Real>(QiStruct, 1e-12);

        auto& bc = m_wk.at(tag);

        bc.area = areaOf[tag];
        bc.q0 = Qi * fVaso;
        bc.Ra = dPa / (Qi * fVaso);
        bc.Rv = dPv / (Qi * fVaso);
        bc.C = std::max<Real>(m_cfg.coronaryComplianceTotal * w, 1e-300);

        // Derived operating point and bed multiplicities (predictions).
        bc.gammaA = gammaA;
        bc.gammaV = gammaV;
        bc.Na = Qi / (PI * ra * ra * va);
        bc.Nv = Qi / (PI * rv * rv * vv);

        // Steady state of the calibrated network, used as initial condition.
        bc.ptm = ptmRest;
        bc.pim = pimRest;
        bc.pc = ptmRest + pimRest;
        bc.pout = bc.pc;
        bc.qd = Qi * fVaso;
        bc.vol0 = bc.C * std::max<Real>(ptmRest, 0.0);
        bc.vol = bc.vol0;
        bc.muA = phiA0;
        bc.muV = phiV0;

        volumeTotal += bc.vol0;
        flowTotal += bc.q0;
        countTotal += bc.Na;

        const Real tau = bc.C * bc.Rv;

        if (isRoot())
          Alert::Info() << "  [calib] outlet " << tag << "  rp=" << (rEq[tag] * 1e3)
                        << " mm"
                        << "  A=" << bc.area << " m^2"
                        << "  Q=" << (Qi * 6.0e7) << " mL/min"
                        << "  Ra=" << bc.Ra << "  Rv=" << bc.Rv << " Pa s/m^3"
                        << "  C=" << bc.C << " m^3/Pa"
                        << "  tau=C*Rv=" << tau << " s"
                        << "  |  derivados: g_a=" << bc.gammaA << " g_v=" << bc.gammaV
                        << " 1/s"
                        << "  N_a=" << bc.Na << " N_v=" << bc.Nv << "  L_a=" << (La * 1e3)
                        << " L_v=" << (Lv * 1e3) << " mm"
                        << "  T_a=" << Ta << " T_v=" << Tv << " s"
                        << "  dPa/dPv=" << dPa << "/" << dPv << " Pa"
                        << "  ptm0=" << ptmRest << " Pa"
                        << "  pc0/pout0=" << bc.pc << "/" << bc.pout
                        << " Pa (par=" << par0 << ", pim=" << pimRest << ")"
                        << Alert::Raise;
      }

      if (isRoot())
      {
        // Flow budget.  With morphometricResistance the total is a prediction
        // and the ratio below is the falsifiable statement of the closure;
        // with the prescribed closure it is an identity and the ratio is 1 by
        // construction.
        Alert::Info() << "  [calib] caudal:"
                      << "  cierre="
                      << (m_cfg.morphometricResistance ? "morfometrico (N_a)"
                                                       : "prescrito (lcaTargetFlow)")
                      << "  N_a(total)=" << countTotal
                      << "  Q(reposo)=" << (flowTotal / fVaso * 6.0e7) << " mL/min"
                      << "  f_vaso=" << fVaso << "  Q(operativo)=" << (flowTotal * 6.0e7)
                      << " mL/min"
                      << "  |  referencia lcaTargetFlow=" << (m_cfg.lcaTargetFlow * 6.0e7)
                      << " mL/min"
                      << "  cociente prediccion/referencia="
                      << (flowTotal / fVaso / std::max<Real>(m_cfg.lcaTargetFlow, 1e-30))
                      << Alert::Raise;

        Alert::Info() << "  [calib] totals:"
                      << "  V0=" << (volumeTotal * 1e6) << " mL"
                      << "  C=" << m_cfg.coronaryComplianceTotal << " m^3/Pa"
                      << "  |  WRMS table nodes=" << m_wrms.logGamma.size()
                      << "  reologia="
                      << (m_cfg.rheologyModel == RheologyModel::Quemada
                             ? "Quemada"
                             : "Carreau-Yasuda")
                      << "  mu_ap(" << gammaA << ")=" << m_wrms(gammaA) << "  mu_ap("
                      << gammaV << ")=" << m_wrms(gammaV)
                      << "  (g_0 en reposo=" << gammaA0 << "/" << gammaV0
                      << ", f_shear=" << fShear << ")"
                      << "  (mu_N=" << muN << ", mu_inf=" << m_cfg.viscosity.muInf
                      << ", mu_0=" << m_cfg.viscosity.mu0 << ")"
                      << (m_cfg.constantOutletResistance
                             ? "  [CONSTANT outlet closure: mu_ap frozen at the "
                               "high-shear plateau; 3D rheology unchanged]"
                             : "")
                      << "  Phi_a0=" << phiA0 << "  Phi_v0=" << phiV0
                      << "  |  outlet resistance is assembled implicitly: the "
                      << "coupling is unconditionally stable in dt" << Alert::Raise;

        // Consistency check. (r, v) fix two degrees of freedom and the
        // pressure budget closes L, so the transit time is over-determined and
        // must agree with indicator dilution. It is a statement about the
        // lumping, not a numerical fault: if it fails, a single (r, L, N)
        // cannot carry its share of the pressure drop at that calibre and that
        // velocity, and the tree needs more than one generation.
        const Real Tref = m_cfg.referenceTransitTime;

        Alert::Info() << "  [calib] comprobacion de tiempo de transito:"
                      << "  T_a=" << Ta << " s  T_v=" << Tv << " s"
                      << "  (referencia " << Tref << " s)"
                      << "  ->  T_a/T_ref=" << (Ta / Tref)
                      << "  T_v/T_ref=" << (Tv / Tref) << Alert::Raise;

        if (Ta < 0.5 * Tref || Ta > 2.0 * Tref || Tv < 0.5 * Tref || Tv > 2.0 * Tref)
          Alert::Warning()
            << "  [calib] el tiempo de transito derivado se aparta mas de 2x "
            << "de la referencia: el calibre, la velocidad y el reparto de "
            << "presion no son mutuamente consistentes con un solo tramo "
            << "efectivo por rama." << Alert::Raise;
      }
    }
    else
    {
      // Diagnostic path: neutral, non-degenerate constants so the outlet still
      // integrates without the calibration.
      for (const Attribute tag : m_cfg.outlets)
      {
        auto& bc = m_wk.at(tag);
        const Real Qi = m_cfg.lcaTargetFlow / static_cast<Real>(m_cfg.outlets.size());
        bc.area = 1.0e-5;
        bc.q0 = Qi;
        bc.Ra = 4.5e9;
        bc.Rv = 6.8e8;
        bc.C = m_cfg.coronaryComplianceTotal / static_cast<Real>(m_cfg.outlets.size());
        bc.vol0 = bc.C * std::max<Real>(bc.ptm, 0.0);
        bc.vol = bc.vol0;
      }
    }

    m_flowSolver.setTolerances(1e-10, 1e-8, 1e-10, 50, 10000)
      .setStateUpdate([this](const PETSc::Math::Vector& x) {
        m_u.getSolution().setData(x, 0);
        m_p.getSolution().setData(x, m_uh.getSize());
      });
  }

  void CoupledLV0DCoronary3D::setupDiagnostics()
  {
    if (!isRoot())
      return;

    m_csv.open(m_cfg.csvPath);

    if (!m_csv)
      throw std::runtime_error("Failed to open coupled CSV file: " + m_cfg.csvPath);

    writeCSVHeader();
  }

  void CoupledLV0DCoronary3D::printInitialState() const
  {
    if (!isRoot())
      return;

    const auto& s = m_model.getState();

    Alert::Info() << "Initial 0D state:" << Alert::NewLine << "y     = " << s.y
                  << Alert::NewLine << "v     = " << s.v << Alert::NewLine
                  << "pv    = " << s.pv << Alert::NewLine << "par   = " << s.par
                  << Alert::NewLine << "pd    = " << s.pd << Alert::NewLine
                  << "ec    = " << s.ec << Alert::NewLine << "gamma = " << s.gamma
                  << Alert::NewLine << "beta  = " << s.beta << Alert::NewLine
                  << "w     = " << s.w << Alert::NewLine << "kc    = " << s.kc
                  << Alert::NewLine << "tauc  = " << s.tauc << Alert::NewLine
                  << "3D flow mode: " << flowModeName(m_cfg.flowMode) << Alert::Raise;
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
                  << "  iter = " << rep.iterations << "  |F| = " << rep.finalResidual
                  << "  |dx| = " << rep.finalStepNorm << Alert::Raise;
    }

    if (isRoot() && rep.converged)
    {
      const auto& s = m_model.getState();
      ZeroDInfo() << "pv  = " << s.pv << " Pa" << "  par = " << s.par << " Pa"
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

    const auto& s = m_model.getState();
    const Real pin = s.par;

    const Real mu = m_cfg.viscosity.mu0;

    const auto symDU = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));

    const auto symV = 0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));

    const auto& uLag = m_uOld;
    const auto oseenConvection = Mult(Jacobian(m_u), uLag); // (uLag·∇)du
    const auto divLag = Div(uLag);
    const auto oseenTemam = divLag * Dot(m_u, m_v); // Temam stabilization

  // Backflow damping (keeps BCs stable even in steady solve)
    const auto outletBeta = Max(-Dot(m_uOld, normal), 0.0);
    const auto outletBackflowDamping =
      0.5 * m_cfg.outletBackflowStabilization * m_cfg.rho * outletBeta;

    // Implicit outlet resistance.
    //
    // The Neumann traction at an outlet is p_out = p_c + R_a Phi_a Q, and for a
    // flat profile the resistive part is
    //
    //   int_G R_a Q (v.n) = R_a A int_G (u.n)(v.n),
    //
    // which is exact for a flat profile, symmetric and positive semidefinite
    // for any profile, and is assembled with exactly the pattern already used
    // by the inlet normal impedance. Assembling it here instead of lagging it
    // by one step is what makes the coupling unconditionally stable: the
    // amplification factor goes from |1 - R dt / L_3D| -- which diverges for
    // the per-branch arteriolar resistance, |g| = 5.4 at dt = 1e-3 -- to
    // 1/(1 + R dt / L_3D) < 1 for any dt and any R. With that, the outlet
    // capacitor C_a is no longer needed, and neither is the characteristic
    // impedance that damped it; removing them removes the L_3D-C_a resonance
    // (f_0 ~ 51 Hz, zeta = 0.25, Q = 2) that produced the observed ringing and
    // the spurious systolic peak in the outlet flux.
    const auto outletResistance = [this](const Attribute tag) {
      const auto& bc = m_wk.at(tag);
      return bc.Ra * bc.muA * bc.area;
    };

    m_flow = 2.0 * Integral(mu * symDU, symV)

      - Integral(m_p, Div(m_v)) + Integral(Div(m_u), m_q) + m_cfg.eps * Integral(m_p, m_q)

      + BoundaryIntegral(pin * Dot(m_v, normal)).over(m_cfg.inlet)

      + BoundaryIntegral(m_wk.at(outlet0).pout * Dot(m_v, normal)).over(outlet0) +
      BoundaryIntegral(m_wk.at(outlet1).pout * Dot(m_v, normal)).over(outlet1) +
      BoundaryIntegral(m_wk.at(outlet2).pout * Dot(m_v, normal)).over(outlet2) +
      BoundaryIntegral(m_wk.at(outlet3).pout * Dot(m_v, normal)).over(outlet3) +
      BoundaryIntegral(m_wk.at(outlet4).pout * Dot(m_v, normal)).over(outlet4) +
      BoundaryIntegral(m_wk.at(outlet5).pout * Dot(m_v, normal)).over(outlet5)

      + outletResistance(outlet0) *
        BoundaryIntegral(Dot(Dot(m_u, normal) * normal, m_v)).over(outlet0) +
      outletResistance(outlet1) *
        BoundaryIntegral(Dot(Dot(m_u, normal) * normal, m_v)).over(outlet1) +
      outletResistance(outlet2) *
        BoundaryIntegral(Dot(Dot(m_u, normal) * normal, m_v)).over(outlet2) +
      outletResistance(outlet3) *
        BoundaryIntegral(Dot(Dot(m_u, normal) * normal, m_v)).over(outlet3) +
      outletResistance(outlet4) *
        BoundaryIntegral(Dot(Dot(m_u, normal) * normal, m_v)).over(outlet4) +
      outletResistance(outlet5) *
        BoundaryIntegral(Dot(Dot(m_u, normal) * normal, m_v)).over(outlet5)

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

    if (isRoot())
    {
      KSPInfo() << "Static solve: " << (reason > 0 ? "Converged" : "Did NOT converge")
                << "  iterations = " << iterations << Alert::Raise;
    }
  }

  bool CoupledLV0DCoronary3D::solve3D()
  {
    const auto setup3DStart = CoronaryClock::now();

    const auto& s = m_model.getState();
    const Real pin = s.par;

    const auto normal = BoundaryNormal(m_mesh);
    const Attribute outlet0 = m_cfg.outlets[0];
    const Attribute outlet1 = m_cfg.outlets[1];
    const Attribute outlet2 = m_cfg.outlets[2];
    const Attribute outlet3 = m_cfg.outlets[3];
    const Attribute outlet4 = m_cfg.outlets[4];
    const Attribute outlet5 = m_cfg.outlets[5];

    const auto& uState = m_u.getSolution();
    const auto& pState = m_p.getSolution();
    const auto& uLag = m_uOld;

    const auto gradDU = Jacobian(m_u);
    const auto gradU = Jacobian(uState);
    const auto gradLag = Jacobian(uLag);

    const auto newtonConvection = Mult(gradDU, uState) + Mult(gradU, m_u);

    const auto stateConvection = Mult(gradU, uState);
    const auto oseenConvectionJacobian = Mult(gradDU, uLag);

    const auto divDU = Div(m_u);
    const auto divU = Div(uState);
    const auto divLag = Div(uLag);

    const auto temamJacobian1 = Dot(divDU * uState, m_v);

    const auto temamJacobian2 = divU * Dot(m_u, m_v);

    const auto temamResidual = divU * Dot(uState, m_v);

    const auto oseenTemamJacobian = divLag * Dot(m_u, m_v);

    const auto outletBeta = Max(-Dot(m_uOld, normal), 0.0);
    const auto inletBeta = Max(Dot(m_uOld, normal), 0.0);
    const auto outletBackflowDamping =
      0.5 * m_cfg.outletBackflowStabilization * m_cfg.rho * outletBeta;
    const auto inletBackflowDamping =
      0.5 * m_cfg.inletBackflowStabilization * m_cfg.rho * inletBeta;
    const auto symDU = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));

    const auto symV = 0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));

    const auto symU = 0.5 * (Jacobian(uState) + Transpose(Jacobian(uState)));

    const auto symLag = 0.5 * (gradLag + Transpose(gradLag));

    const auto duNormal = Dot(m_u, normal) * normal;

    const auto duTangential = m_u - duNormal;

    const auto uStateNormal = Dot(uState, normal) * normal;

    const auto uStateTangential = uState - uStateNormal;

    // Implicit outlet resistance coefficient, R_a Phi_a A. See solveStatic().
    const auto outletResistance = [this](const Attribute tag) {
      const auto& bc = m_wk.at(tag);
      return bc.Ra * bc.muA * bc.area;
    };

    const auto& cy = m_cfg.viscosity;
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

    const auto muLag = muInf + deltaMu * Pow(carreauBaseLag, (nCY - 1.0) / yasuda);

    const auto dgamma = 2.0 * Dot(symU, symDU) / gamma;

    const auto dmu = deltaMu * (nCY - 1.0) *
      Pow(carreauBase, (nCY - 1.0 - yasuda) / yasuda) * std::pow(lambda, yasuda) *
      Pow(gamma, yasuda - 1.0) * dgamma;

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

    m_l2ConvU = Integral(m_up, m_vp) - Integral(convectionTarget, m_vp);
    m_l2ConvU.assemble();
    m_l2ConvUSolver.solve();

    RealFunction tau = [=, this](const Point& p) -> Real {
      const auto uOld = m_uOld.getValue(p);
      const Real mu = m_mu.getSolution().getValue(p);
      const Real hK =
        std::pow(p.getPolytope().getMeasure(), 1.0 / p.getPolytope().getDimension());
      const Real order = 2;
      const Real speed = std::sqrt(dot(uOld, uOld));

      const Real Tau = 1. /
        (4.0 * std::pow(order, 4.) * mu / (m_cfg.rho * std::pow(hK, 2.0)) +
          2.0 * order * speed / hK);

      return 1. / (m_cfg.rho / m_cfg.dt + m_cfg.rho / Tau);
    };

    m_tauProjection = Integral(m_tau, m_t) - Integral(tau, m_t);
    m_tauProjection.assemble();
    m_tauProjectionSolver.solve();

    auto subUpdate = VectorFunction(
      m_mesh.getSpaceDimension(), [=, this](const Point& p) -> Math::SpatialVector<Real> {
        const auto conv = convectionTarget.getValue(p);
        const auto proj = m_up.getSolution().getValue(p);
        const auto old = m_subOld.getValue(p);
        const auto tau = m_tau.getSolution().getValue(p);

        Math::SpatialVector<Real> out(m_mesh.getSpaceDimension());

        for (size_t c = 0; c < out.size(); ++c)
        {
          out[c] = tau * m_cfg.rho * (1. / m_cfg.dt * old[c] - (conv[c] - proj[c]));
        }

        return out;
      });

  /*
   * L2 projection of the dynamic subscale into m_sub:
   *
   *   int sub_h · v = int subUpdate · v.
   */

    m_subProjection = Integral(m_sub, m_vp) - Integral(subUpdate, m_vp);
    m_subProjection.assemble();
    m_subProjectionSolver.solve();

    if (m_cfg.flowMode == FlowMode::Newton)
    {
      m_flow =
        /*
         * =========================
         * Newton Jacobian / tangent
         * =========================
         *
         * Unknowns here are Newton corrections:
         *   m_u : velocity correction
         *   m_p : pressure correction
         */

        (m_cfg.rho / m_cfg.dt) * Integral(m_u, m_v)

        /*
         * Newton linearization of convection:
         *   rho ((du · grad) uState + (uState · grad) du, v)
         *
         * In your notation this is encoded by newtonConvection.
         */
        + m_cfg.rho * Integral(Dot(newtonConvection, m_v))

        /*
         * Temam/skew-symmetric convection correction tangent.
         */
        + 0.5 * m_cfg.rho * Integral(temamJacobian1) +
        0.5 * m_cfg.rho * Integral(temamJacobian2)

        /*
         * Nonlinear viscous tangent:
         *   2 mu D(du) : D(v)
         * + 2 dmu[du] D(uState) : D(v)
         */
        + 2.0 * Integral(mu * symDU, symV) +
        2.0 * Integral(dmu * symU, symV)

        /*
         * Stokes pressure/divergence block.
         */
        - Integral(m_p, Div(m_v)) + Integral(Div(m_u), m_q) +
        m_cfg.eps * Integral(m_p, m_q)

        /*
         * Inlet normal impedance tangent.
         *
         * Boundary pressure law:
         *   p_inlet_boundary = pin + Z (u · n)
         *
         * Tangent:
         *   Z (du · n) (v · n)
         */
        + m_cfg.inletImpedance *
          BoundaryIntegral(Dot(Dot(m_u, normal) * normal, m_v)).over(m_cfg.inlet)

        /*
         * Inlet tangential damping tangent.
         *
         * Controls u_tau without directly prescribing normal inflow.
         */
        + m_cfg.inletTangentialDamping *
          BoundaryIntegral(Dot(duTangential, m_v)).over(m_cfg.inlet)

        /*
         * Backflow stabilization tangent.
         *
         * inletBeta activates when inlet behaves as an outlet:
         *   uOld · n > 0
         *
         * outletBeta activates when outlet behaves as an inlet:
         *   uOld · n < 0
         *
         * Since beta is lagged with m_uOld, the tangent is simply beta * du.
         */
        + BoundaryIntegral(inletBackflowDamping * Dot(m_u, m_v)).over(m_cfg.inlet)

        + BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet0) +
        BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet1) +
        BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet2) +
        BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet3) +
        BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet4) +
        BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet5)

        /*
         * Implicit outlet resistance tangent.
         *
         * Boundary pressure law:
         *   p_out = p_c + R_a Phi_a Q,   Q = int (u.n)
         *
         * Tangent (flat-profile form, symmetric positive semidefinite):
         *   R_a Phi_a A (du.n)(v.n)
         *
         * Assembling it rather than lagging it is what makes the 3D-0D
         * coupling unconditionally stable, and is what removes the need for the
         * outlet capacitor and, with it, the L_3D-C_a ringing.
         */
        + outletResistance(outlet0) * BoundaryIntegral(Dot(duNormal, m_v)).over(outlet0) +
        outletResistance(outlet1) * BoundaryIntegral(Dot(duNormal, m_v)).over(outlet1) +
        outletResistance(outlet2) * BoundaryIntegral(Dot(duNormal, m_v)).over(outlet2) +
        outletResistance(outlet3) * BoundaryIntegral(Dot(duNormal, m_v)).over(outlet3) +
        outletResistance(outlet4) * BoundaryIntegral(Dot(duNormal, m_v)).over(outlet4) +
        outletResistance(outlet5) * BoundaryIntegral(Dot(duNormal, m_v)).over(outlet5)

        /*
         * =========================
         * Newton residual
         * =========================
         *
         * Everything below must use the current nonlinear state:
         *   uState, pState
         *
         * Do not use m_u here, except inside DirichletBC correction terms.
         */

        + (m_cfg.rho / m_cfg.dt) * Integral(uState, m_v) -
        (m_cfg.rho / m_cfg.dt) * Integral(m_uOld, m_v)

        /*
         * Convective residual:
         *   rho ((uState · grad) uState, v)
         */
        + m_cfg.rho * Integral(Dot(stateConvection, m_v))

        /*
         * Temam/skew residual.
         */
        + 0.5 * m_cfg.rho * Integral(temamResidual)

        /*
         * Viscous residual:
         *   2 mu(uState) D(uState) : D(v)
         */
        + 2.0 * Integral(mu * symU, symV)

        /*
         * Pressure/divergence residual.
         */
        - Integral(pState, Div(m_v)) + Integral(Div(uState), m_q) +
        m_cfg.eps * Integral(pState, m_q)

        /*
         * Inlet pressure source.
         */
        + BoundaryIntegral(pin * Dot(m_v, normal)).over(m_cfg.inlet)

        /*
         * Inlet normal impedance residual.
         *
         * Must use uState, not m_u.
         */
        //+ m_cfg.inletImpedance *
        //    BoundaryIntegral(Dot(Dot(uState, normal) * normal, m_v))
        //      .over(m_cfg.inlet)

        /*
         * Inlet tangential damping residual.
         *
         * Must use uStateTangential, not duTangential.
         */
        // + m_cfg.inletTangentialDamping *
        //     BoundaryIntegral(Dot(uStateTangential, m_v)).over(m_cfg.inlet)

        /*
         * Outlet pressure Neumann residuals.
         */
        + BoundaryIntegral(m_wk.at(outlet0).pout * Dot(m_v, normal)).over(outlet0) +
        BoundaryIntegral(m_wk.at(outlet1).pout * Dot(m_v, normal)).over(outlet1) +
        BoundaryIntegral(m_wk.at(outlet2).pout * Dot(m_v, normal)).over(outlet2) +
        BoundaryIntegral(m_wk.at(outlet3).pout * Dot(m_v, normal)).over(outlet3) +
        BoundaryIntegral(m_wk.at(outlet4).pout * Dot(m_v, normal)).over(outlet4) +
        BoundaryIntegral(m_wk.at(outlet5).pout * Dot(m_v, normal)).over(outlet5)

        /*
         * Implicit outlet resistance residual.
         *
         * Must use uStateNormal, not duNormal.
         */
        + outletResistance(outlet0) *
          BoundaryIntegral(Dot(uStateNormal, m_v)).over(outlet0) +
        outletResistance(outlet1) *
          BoundaryIntegral(Dot(uStateNormal, m_v)).over(outlet1) +
        outletResistance(outlet2) *
          BoundaryIntegral(Dot(uStateNormal, m_v)).over(outlet2) +
        outletResistance(outlet3) *
          BoundaryIntegral(Dot(uStateNormal, m_v)).over(outlet3) +
        outletResistance(outlet4) *
          BoundaryIntegral(Dot(uStateNormal, m_v)).over(outlet4) +
        outletResistance(outlet5) * BoundaryIntegral(Dot(uStateNormal, m_v)).over(outlet5)

        /*
         * Backflow stabilization residual.
         *
         * Must use uState.
         */
        + BoundaryIntegral(inletBackflowDamping * Dot(uState, m_v)).over(m_cfg.inlet)

        + BoundaryIntegral(outletBackflowDamping * Dot(uState, m_v)).over(outlet0) +
        BoundaryIntegral(outletBackflowDamping * Dot(uState, m_v)).over(outlet1) +
        BoundaryIntegral(outletBackflowDamping * Dot(uState, m_v)).over(outlet2) +
        BoundaryIntegral(outletBackflowDamping * Dot(uState, m_v)).over(outlet3) +
        BoundaryIntegral(outletBackflowDamping * Dot(uState, m_v)).over(outlet4) +
        BoundaryIntegral(outletBackflowDamping * Dot(uState, m_v)).over(outlet5)

        /*
         * Variational elimination of wall Dirichlet condition.
         *
         * Since unknown is Newton correction, impose:
         *   du = -uState
         * on the wall.
         */
        + DirichletBC(m_u, -uState).on(m_cfg.wall);
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

        (m_cfg.rho / m_cfg.dt) * Integral(m_u, m_v)

        /*
         * Oseen convection with lagged convecting velocity uLag = m_uOld.
         */
        + m_cfg.rho * Integral(Dot(oseenConvectionJacobian, m_v))

        /*
         * Lagged Temam correction.
         */
        + 0.5 * m_cfg.rho * Integral(oseenTemamJacobian)

        /*
         * Lagged viscosity.
         */
        + 2.0 * Integral(muLag * symDU, symV)

        /*
         * Stokes pressure/divergence block.
         */
        - Integral(m_p, Div(m_v)) + Integral(Div(m_u), m_q) +
        m_cfg.eps * Integral(m_p, m_q)

        /*
         * VMS bilinear contribution:
         *
         *   int_K tau_K rho^2
         *     ((grad u^{n+1}) u^n)
         *     ·
         *     ((grad v) u^n).
         */

        //  + VMSConvectionBilinearIntegrator(m_u, m_v, m_uOld, m_tau.getSolution(),
        //                                  m_cfg.rho)

        /*
         * VMS linear contribution subtracted from the residual:
         *
         *   - int_K rho
         *       (tau_K rho u_proj + u'^{n+1})
         *       ·
         *       ((grad v) u^n).
         */

        //- VMSConvectionLinearIntegrator(m_v, m_sub.getSolution(), m_uOld,
        //                              m_up.getSolution(), m_tau.getSolution(),
        //                            m_cfg.rho, m_cfg.dt)

        /*
         * Inlet normal impedance.
         */
        + m_cfg.inletImpedance *
          BoundaryIntegral(Dot(Dot(m_u, normal) * normal, m_v)).over(m_cfg.inlet)

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

        + BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet0) +
        BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet1) +
        BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet2) +
        BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet3) +
        BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet4) +
        BoundaryIntegral(outletBackflowDamping * Dot(m_u, m_v)).over(outlet5)

        /*
         * Implicit outlet resistance.
         *
         * p_out = p_c + R_a Phi_a Q with Q = int (u.n); the resistive part is
         * assembled here rather than lagged, which makes the coupling
         * unconditionally stable and removes the outlet capacitor together
         * with the L_3D-C_a resonance it produced.
         */
        + outletResistance(outlet0) * BoundaryIntegral(Dot(duNormal, m_v)).over(outlet0) +
        outletResistance(outlet1) * BoundaryIntegral(Dot(duNormal, m_v)).over(outlet1) +
        outletResistance(outlet2) * BoundaryIntegral(Dot(duNormal, m_v)).over(outlet2) +
        outletResistance(outlet3) * BoundaryIntegral(Dot(duNormal, m_v)).over(outlet3) +
        outletResistance(outlet4) * BoundaryIntegral(Dot(duNormal, m_v)).over(outlet4) +
        outletResistance(outlet5) * BoundaryIntegral(Dot(duNormal, m_v)).over(outlet5)

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

        - (m_cfg.rho / m_cfg.dt) * Integral(m_uOld, m_v)

        /*
         * Inlet pressure source.
         */
        + BoundaryIntegral(pin * Dot(m_v, normal)).over(m_cfg.inlet)

        /*
         * Explicit outlet pressures.
         */
        + BoundaryIntegral(m_wk.at(outlet0).pout * Dot(m_v, normal)).over(outlet0) +
        BoundaryIntegral(m_wk.at(outlet1).pout * Dot(m_v, normal)).over(outlet1) +
        BoundaryIntegral(m_wk.at(outlet2).pout * Dot(m_v, normal)).over(outlet2) +
        BoundaryIntegral(m_wk.at(outlet3).pout * Dot(m_v, normal)).over(outlet3) +
        BoundaryIntegral(m_wk.at(outlet4).pout * Dot(m_v, normal)).over(outlet4) +
        BoundaryIntegral(m_wk.at(outlet5).pout * Dot(m_v, normal)).over(outlet5)

        /*
         * Wall no-slip condition.
         *
         * Since Oseen unknown is the new velocity itself, impose:
         *   u = 0.
         */
        + DirichletBC(m_u, Zero(m_mesh.getSpaceDimension())).on(m_cfg.wall);
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
      ThreeDInfo() << "Solving with PETSc " << (isNewtonFlow ? "SNES" : "KSP") << " ..."
                   << Alert::Raise;
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
        SNESInfo() << (m_flowSolver.converged() ? "Converged" : "Did NOT converge")
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

    if (isRoot())
    {
      KSPInfo() << (reason > 0 ? "Converged" : "Did NOT converge")
                << "  iterations = " << iterations << Alert::Raise;
    }
    return reason > 0;
  }

  void CoupledLV0DCoronary3D::computeFluxes()
  {
    const auto normal = BoundaryNormal(m_mesh);
    const auto& s = m_model.getState();

    m_flux = BoundaryIntegral(Dot(m_u.getSolution(), normal), m_qFlux).over(m_cfg.inlet);

    m_flux.assemble();
    m_stepData.qIn = m_flux(m_one);

    m_stepData.qOut.clear();
    m_stepData.qOutSum = 0.0;

    for (const Attribute tag : m_cfg.outlets)
    {
      m_flux = BoundaryIntegral(Dot(m_u.getSolution(), normal), m_qFlux).over(tag);

      m_flux.assemble();

      const Real qOut = m_flux(m_one);

      m_stepData.qOut[tag] = qOut;
      m_stepData.qOutSum += qOut;
    }

    m_stepData.flowBalance = m_stepData.qIn + m_stepData.qOutSum;
    m_stepData.t = s.t;
  }

  void CoupledLV0DCoronary3D::computeWallShear()
  {
    const auto& cy = m_cfg.viscosity;
    const auto normal = BoundaryNormal(m_mesh);
    const auto& uSol = m_u.getSolution();

    const auto jacRow0 = VectorFunction(Component(Jacobian(uSol), 0, 0),
      Component(Jacobian(uSol), 0, 1), Component(Jacobian(uSol), 0, 2));
    const auto jacRow1 = VectorFunction(Component(Jacobian(uSol), 1, 0),
      Component(Jacobian(uSol), 1, 1), Component(Jacobian(uSol), 1, 2));
    const auto jacRow2 = VectorFunction(Component(Jacobian(uSol), 2, 0),
      Component(Jacobian(uSol), 2, 1), Component(Jacobian(uSol), 2, 2));

    VelocityGridFunctionType gradRec0(m_uh);
    VelocityGridFunctionType gradRec1(m_uh);
    VelocityGridFunctionType gradRec2(m_uh);

    // Same mass matrix, three right-hand sides: the member problem + KSP keep
    // the matrix and its preconditioner across steps and across the three
    // solves here (see setupProjectionSolvers()).
    m_gradRecProj =
      Integral(m_gradRecTrial, m_gradRecTest) - Integral(jacRow0, m_gradRecTest);
    m_gradRecProj.solve(m_gradRecKSP);
    gradRec0.setData(m_gradRecTrial.getSolution().getData());

    m_gradRecProj =
      Integral(m_gradRecTrial, m_gradRecTest) - Integral(jacRow1, m_gradRecTest);
    m_gradRecProj.solve(m_gradRecKSP);
    gradRec1.setData(m_gradRecTrial.getSolution().getData());

    m_gradRecProj =
      Integral(m_gradRecTrial, m_gradRecTest) - Integral(jacRow2, m_gradRecTest);
    m_gradRecProj.solve(m_gradRecKSP);
    gradRec2.setData(m_gradRecTrial.getSolution().getData());

    // (grad u . n)_i = gradRec_i . n is the wall-normal directional derivative.
    const auto gradUn0 = Dot(gradRec0, normal);
    const auto gradUn1 = Dot(gradRec1, normal);
    const auto gradUn2 = Dot(gradRec2, normal);

    // Local Carreau-Yasuda viscosity, built from the SAME expression projected
    // into m_mu (the "viscosity" field written to XDMF). The wall shear stress
    // is tau_w = mu(gamma) * du/dn; using the constant muInf here decoupled the
    // WSS from the rheology entirely, so the WSS map could not reflect the
    // viscosity field (no correlation) and under-predicted tau_w in exactly the
    // low-shear recirculation zones where mu departs from muInf.
    const auto symU = 0.5 * (Jacobian(uSol) + Transpose(Jacobian(uSol)));
    const auto gammaDot =
      Sqrt(cy.gammaRegularization * cy.gammaRegularization + 2.0 * Dot(symU, symU));
    const auto carreauBase = 1.0 + Pow(cy.lambda * gammaDot, cy.yasuda);
    const auto muLocal =
      cy.muInf + (cy.mu0 - cy.muInf) * Pow(carreauBase, (cy.n - 1.0) / cy.yasuda);

    const auto tracRec =
      VectorFunction(muLocal * gradUn0, muLocal * gradUn1, muLocal * gradUn2);
    const auto wallStressRec = tracRec - Dot(tracRec, normal) * normal;

    const Real wssReg = 1.0e-3;
    m_wssProj = BoundaryIntegral(Dot(m_wssTrial, m_wssTest)).over(m_cfg.wall) +
      wssReg * Integral(Dot(m_wssTrial, m_wssTest)) -
      BoundaryIntegral(Dot(wallStressRec, m_wssTest)).over(m_cfg.wall);
    m_wssProj.solve(m_wssKSP);
    m_shearWall.setData(m_wssTrial.getSolution().getData());

    if (isRoot())
    {
      ::Vec fvec = m_shearWall.getData();
      PetscInt nf = 0;
      VecGetLocalSize(fvec, &nf);
      const PetscScalar* farr = nullptr;
      VecGetArrayRead(fvec, &farr);
      Real fieldMax = 0.0;
      for (PetscInt i = 0; i < nf; ++i)
        fieldMax = std::max(fieldMax, std::abs(farr[i]));
      VecRestoreArrayRead(fvec, &farr);
      Alert::Info() << "  [wss] field max = " << fieldMax << " Pa" << Alert::Raise;
    }
  }

  void CoupledLV0DCoronary3D::updateHistory()
  {
    m_uOld.setData(m_u.getSolution().getData());
    m_pOld.setData(m_p.getSolution().getData());
    m_subOld.setData(m_sub.getSolution().getData());
    m_tauOld.setData(m_tau.getSolution().getData());
  }

  void CoupledLV0DCoronary3D::writeOutputs()
  {
    const auto& cy = m_cfg.viscosity;
    const auto symU =
      0.5 * (Jacobian(m_u.getSolution()) + Transpose(Jacobian(m_u.getSolution())));
    const auto gamma =
      Sqrt(cy.gammaRegularization * cy.gammaRegularization + 2.0 * Dot(symU, symU));
    const auto carreauBase = 1.0 + Pow(cy.lambda * gamma, cy.yasuda);
    const auto mu =
      cy.muInf + (cy.mu0 - cy.muInf) * Pow(carreauBase, (cy.n - 1.0) / cy.yasuda);

    m_viscosityProjection = Integral(m_mu, m_r) - Integral(mu, m_r);
    m_viscosityProjection.solve(m_viscosityProjectionKSP);

    computeWallShear();

    m_xdmf.write(m_model.getState().t - m_cfg.dt).flush();
  }

  CoupledLV0DCoronary3D::StepData CoupledLV0DCoronary3D::collectStepData() const
  {
    StepData d;

    const auto& s = m_model.getState();

    d.t = s.t;
    d.pat = m_input.pAt(s.t);
    d.psv = m_input.pSv(s.t);

    d.y = s.y;
    d.v = s.v;
    d.radius = m_input.R0 + s.y;

    d.lvVolume = (4.0 / 3.0) * std::numbers::pi_v<Real> * d.radius * d.radius * d.radius;

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

    Real muASum = 0.0;
    Real muVSum = 0.0;
    Real ptmSum = 0.0;
    Real n = 0.0;

    for (const auto& [tag, bc] : m_wk)
    {
      const auto qOutIt = m_stepData.qOut.find(tag);
      const Real qOut = (qOutIt == m_stepData.qOut.end()) ? 0.0 : qOutIt->second;

      d.qDistal[tag] = bc.qd;
      d.qDistalSum += bc.qd;
      d.qCapChargingSum += qOut - bc.qd;
      d.pc[tag] = bc.pc;
      d.pOut[tag] = bc.pout;
      // ~identical across outlets (same pv)
      d.pim = bc.pim;

      // Mechanism diagnostics. p_tm is the state and C p_tm the stored volume
      // (the pump); the two viscosity ratios are the rheological modulation of
      // each limb, which is where a hyperviscous condition shows up. p_tm must
      // stay positive over the whole cycle: it is the model's own check that
      // the Starling throat is doing its job.
      muASum += bc.muA;
      muVSum += bc.muV;
      ptmSum += bc.ptm;
      n += 1.0;

      d.storedVolume += bc.vol;
    }

    if (n > 0.0)
    {
      d.muARatio = muASum / n;
      d.muVRatio = muVSum / n;
      d.ptm = ptmSum / n;
    }

    return d;
  }

  void CoupledLV0DCoronary3D::writeCSVHeader()
  {
    if (!isRoot())
      return;

    m_csv << "t," << "LeftAtriumPressure," << "VenaCavaPressure,"
          << "LeftVentricleDisplacement," << "LeftVentricleVelocity,"
          << "LeftVentricleRadius," << "LeftVentricleVolume," << "LeftVentriclePressure,"
          << "AortaPressure," << "DistalPressure," << "LeftVentricleFlow,"
          << "CoronaryInletFlux," << "CoronaryOutlet7Flux," << "CoronaryOutlet8Flux,"
          << "CoronaryOutlet9Flux," << "CoronaryOutlet10Flux," << "CoronaryOutlet14Flux,"
          << "CoronaryOutlet15Flux," << "CoronaryOutletFluxTotal,"
          << "CoronaryOutlet7DistalFlux," << "CoronaryOutlet8DistalFlux,"
          << "CoronaryOutlet9DistalFlux," << "CoronaryOutlet10DistalFlux,"
          << "CoronaryOutlet14DistalFlux," << "CoronaryOutlet15DistalFlux,"
          << "CoronaryDistalFluxTotal," << "CoronaryCapChargingFluxTotal,"
          << "CoronaryOutlet7CapPressure," << "CoronaryOutlet8CapPressure,"
          << "CoronaryOutlet9CapPressure," << "CoronaryOutlet10CapPressure,"
          << "CoronaryOutlet14CapPressure," << "CoronaryOutlet15CapPressure,"
          << "CoronaryOutlet7Pressure," << "CoronaryOutlet8Pressure,"
          << "CoronaryOutlet9Pressure," << "CoronaryOutlet10Pressure,"
          << "CoronaryOutlet14Pressure," << "CoronaryOutlet15Pressure," << "FlowBalance,"
          << "ec," << "gamma," << "beta," << "w," << "kc," << "tauc,"
          << "IntramyoPressure," << "TransmuralPressure,"
          << "ArteriolarViscosityRatio," << "VenularViscosityRatio,"
          << "MicrovascularVolume\n";

    m_csv.flush();
  }

  void CoupledLV0DCoronary3D::writeCSVRow()
  {
    if (!isRoot())
      return;

    const StepData d = collectStepData();

    auto get = [](const std::map<Attribute, Real>& m, Attribute a) -> Real {
      const auto it = m.find(a);
      return (it == m.end()) ? 0.0 : it->second;
    };

    m_csv << d.t << ',' << d.pat << ',' << d.psv << ',' << d.y << ',' << d.v << ','
          << d.radius << ',' << d.lvVolume << ',' << d.pv << ',' << d.par << ',' << d.pd
          << ',' << d.lvFlow << ',' << d.qIn << ',' << get(d.qOut, 7) << ','
          << get(d.qOut, 8) << ',' << get(d.qOut, 9) << ',' << get(d.qOut, 10) << ','
          << get(d.qOut, 14) << ',' << get(d.qOut, 15) << ',' << d.qOutSum << ','
          << get(d.qDistal, 7) << ',' << get(d.qDistal, 8) << ',' << get(d.qDistal, 9)
          << ',' << get(d.qDistal, 10) << ',' << get(d.qDistal, 14) << ','
          << get(d.qDistal, 15) << ',' << d.qDistalSum << ',' << d.qCapChargingSum << ','
          << get(d.pc, 7) << ',' << get(d.pc, 8) << ',' << get(d.pc, 9) << ','
          << get(d.pc, 10) << ',' << get(d.pc, 14) << ',' << get(d.pc, 15) << ','
          << get(d.pOut, 7) << ',' << get(d.pOut, 8) << ',' << get(d.pOut, 9) << ','
          << get(d.pOut, 10) << ',' << get(d.pOut, 14) << ',' << get(d.pOut, 15) << ','
          << d.flowBalance << ',' << d.ec << ',' << d.gamma << ',' << d.beta << ',' << d.w
          << ',' << d.kc << ',' << d.tauc << ',' << d.pim << ',' << d.ptm << ','
          << d.muARatio << ',' << d.muVRatio << ',' << d.storedVolume << '\n';

    m_csv.flush();
  }

  void CoupledLV0DCoronary3D::printStepTiming(int step) const
  {
    const auto& comm = m_mesh.getContext().getCommunicator();
    auto maxTime = [&](Real local) {
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
                   << "  0D = " << advance0D << " s" << "  3D-form = " << setup3DForm
                   << " s" << "  3D-assemble = " << assemble3D << " s"
                   << "  field-splits = " << fieldSplits << " s"
                   << "  3D-solve = " << solve3D << " s" << "  fluxes = " << fluxes
                   << " s" << "  RCR = " << outletRCR << " s" << "  csv = " << csv << " s"
                   << "  history = " << history << " s" << "  output = " << output << " s"
                   << Alert::Raise;
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
      Alert::Exception() << "Invalid time adaptivity maximum levels: " << maxLevels
                         << ". Expected a non-negative value." << Alert::Raise;
      return 1;
    }

    Real nextDt = baseDt;
    int acceptedStep = 0;

    {
      // solveStatic();

      // computeFluxes();
      // for (const Attribute tag : m_cfg.outlets)
      // updateOutlet0D(m_cfg, m_wrms, m_model, m_wk[tag],
      //  m_stepData.qOut.at(tag), m_cfg.dt);

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

      m_cfg.dt = physicalDt;
      m_stepTiming = StepTiming{};
      const auto stepStart = CoronaryClock::now();

      if (isRoot())
      {
        Alert::Info() << "━━━ Step " << (acceptedStep + 1) << "  t = " << t_current
                      << " s" << " / " << finalTime << " s" << "  (dt = " << physicalDt
                      << " s)" << " ━━━" << Alert::Raise;
      }

      const auto advance0DStart = CoronaryClock::now();
      const bool advanced0D = advance0D();
      m_stepTiming.advance0D = secondsSince(advance0DStart);

      if (!advanced0D)
      {
        Alert::Exception() << "0D solver failed to converge at step "
                           << (acceptedStep + 1) << ", t = " << m_model.getState().t
                           << Alert::Raise;
        return 1;
      }

      const auto advancedState = m_model.getState();
      const auto advancedHistory = m_model.getHistory();
      const auto advancedUnknowns = m_model.getUnknowns();
      const auto advancedReport = m_model.getReport();

      bool accepted = false;

      while (!accepted)
      {
        m_cfg.dt = solverDt;

        if (!solve3D())
        {
          m_model.restore(
            advancedState, advancedHistory, advancedUnknowns, advancedReport);
          m_wk = savedRCR;
          m_stepData = savedStepData;
          savedU.restore(m_u.getSolution().getData());
          savedP.restore(m_p.getSolution().getData());

          if (level >= maxLevels)
          {
            Alert::Exception() << "3D flow solver failed to converge at step "
                               << (acceptedStep + 1) << " after " << (level + 1)
                               << " attempt(s). Minimum solver dt = " << solverDt << " s."
                               << Alert::Raise;
            return 1;
          }

          solverDt *= factor;
          ++level;

          if (isRoot())
          {
            Alert::Info() << "[3D] Retrying step " << (acceptedStep + 1)
                          << " with reduced solver dt = " << solverDt
                          << " s  (adapt level = " << level << " / " << maxLevels << ")"
                          << Alert::Raise;
          }

          continue;
        }

        const auto fluxesStart = CoronaryClock::now();
        computeFluxes();
        m_stepTiming.fluxes = secondsSince(fluxesStart);

        const auto outletRCRStart = CoronaryClock::now();
        m_cfg.dt = physicalDt;
        for (const Attribute tag : m_cfg.outlets)
          updateOutlet0D(
            m_cfg, m_wrms, m_model, m_wk[tag], m_stepData.qOut.at(tag), m_cfg.dt);
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
