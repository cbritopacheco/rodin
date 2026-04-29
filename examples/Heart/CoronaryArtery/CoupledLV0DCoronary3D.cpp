#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <type_traits>

#include <Rodin/Alert.h>
#include <Rodin/Solver.h>

#include "CoupledLV0DCoronary3D.h"
#include "Rodin/Math/RungeKutta/RK4.h"
#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Math/RootFinding/NewtonRaphson.h"
#include "VMSConvectionIntegrator.h"

/**
 * @file CoupledLV0DCoronary3D.cpp
 * @brief Partitioned 0D--3D coronary haemodynamics example.
 *
 * This file implements a one-way coupled model:
 *
 *   LV 0D model
 *      -> inlet pressure for the 3D coronary domain
 *      -> semi-implicit generalized-Newtonian Navier--Stokes solve
 *      -> outlet fluxes
 *      -> nonlinear rheological RCR outlet updates.
 *
 * The 3D fluid model is a backward-Euler, Oseen-linearized,
 * incompressible Navier--Stokes problem with:
 *
 *   - Carreau--Yasuda viscosity evaluated explicitly from u^n,
 *   - dynamic projected VMS convective stabilization,
 *   - pressure regularization,
 *   - inlet/outlet pressure tractions,
 *   - outlet backflow stabilization.
 *
 * The outlet model replaces the constant RCR resistances by nonlinear
 * rheological pipe-flow laws derived from the WRMS/Rabinowitsch correction
 * for generalized Newtonian fluids.
 *
 * The coupling is partitioned and one-way at the cardiac/aortic level:
 *
 *   p_ar(t) from the 0D model drives the coronary inlet.
 *
 * The coronary outlet states are updated from the 3D outlet fluxes, but they
 * are not fed back into the LV/aortic 0D model in this driver.
 */

namespace Rodin::Examples::Heart
{
  using namespace Rodin;
  using namespace Rodin::Math;
  using namespace Rodin::Solver;
  using namespace Rodin::Geometry;
  using namespace Rodin::Variational;

  CoupledLV0DCoronary3D::CoupledLV0DCoronary3D()
    : CoupledLV0DCoronary3D(Config{})
  {
  }

  CoupledLV0DCoronary3D::CoupledLV0DCoronary3D(const Config& cfg)
    : m_cfg(cfg),
      m_input(makeInput()),
      m_model(m_input),
      m_mesh(makeMesh(m_cfg)),
      m_xdmf(m_cfg.xdmfBasename),
      m_uh(std::integral_constant<size_t, 2>{}, m_mesh, m_mesh.getSpaceDimension()),
      m_ph(std::integral_constant<size_t, 1>{}, m_mesh),
      m_muh(m_mesh),
      m_tauh(m_mesh),
      m_uph(std::integral_constant<size_t, 2>{}, m_mesh, m_mesh.getSpaceDimension()),
      m_u(m_uh),
      m_p(m_ph),
      m_tau(m_tauh),
      m_t(m_tauh),
      m_tauOld(m_tauh),
      m_mu(m_muh),
      m_up(m_uph),
      m_sub(m_uph),
      m_subOld(m_uph),
      m_v(m_uh),
      m_q(m_ph),
      m_w(m_muh),
      m_vp(m_uph),
      m_uOld(m_uh),
      m_pOld(m_ph),
      m_one(m_ph),
      m_qFlux(m_ph),
      m_flux(m_qFlux),
      m_muProjection(m_mu, m_w),
      m_l2ConvU(m_up, m_vp),
      m_subProjection(m_sub, m_vp),
      m_flow(m_u, m_p, m_v, m_q),
      m_tauProjection(m_tau, m_t),
      m_muProjectionSolver(m_muProjection),
      m_l2ConvUSolver(m_l2ConvU),
      m_subProjectionSolver(m_subProjection),
      m_flowSolver(m_flow),
      m_tauProjectionSolver(m_tauProjection)
  {
    Alert::Info() << "Number of elements in mesh: " << m_mesh.getCellCount() << '\n'
                  << "Number of vertices in mesh: " << m_mesh.getVertexCount() << '\n'
                  << "Mesh space dimension: " << m_mesh.getSpaceDimension() << '\n'
                  << Alert::Raise;

    Alert::Info() << "Velocity space has " << m_uh.getSize() << " DOFs.\n"
                  << "Pressure space has " << m_ph.getSize() << " DOFs.\n"
                  << Alert::Raise;
  }

  CoupledLV0DCoronary3D::~CoupledLV0DCoronary3D() = default;

  CoupledLV0DCoronary3D::MeshType
  CoupledLV0DCoronary3D::makeMesh(const Config& cfg)
  {
    MeshType mesh;
    mesh.load(cfg.meshPath, IO::FileFormat::MEDIT);
    mesh.scale(cfg.meshScale);

    if (mesh.getSpaceDimension() != 3)
      throw std::runtime_error("Expected a 3D coronary mesh.");

    Alert::Info() << "Computing connectivity for " << cfg.meshPath << " ..."
                  << Alert::Raise;

    mesh.getConnectivity().compute(3, 2);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    mesh.getConnectivity().compute(2, 3);

    return mesh;
  }

  CoupledLV0DCoronary3D::Model::Input
  CoupledLV0DCoronary3D::makeInput()
  {
    Model::Input input;

    input.rho = 1.0e3;
    input.R0 = 2.36e-2;
    input.d0 = 1.42e-2;

    input.Es = 3.0e5;
    input.mu = 70.0;
    input.eta = 70.0;
    input.alpha = 3.0;
    input.k0 = 1.0e5;
    input.sigma0 = 5.0e5;

    input.Rp = 8.0e6;
    input.Cp = 5.0e-9;
    input.Rd = 1.0e8;
    input.Cd = 1.0e-8;

    input.Kat = 8.0e-7;
    input.Kp  = 5.0e-10;
    input.Kar = 1.3e-5;

    input.cavityCapacity = 5.0e-12;

    input.localTolerance = 1e-12;
    input.localMaxIterations = 50;
    input.localDamping = 1.0;
    input.absRegularization = 1e-14;

    input.initFibDef = 0.0;
    input.initActiveStiffness = 0.0;
    input.initActiveStress = 0.0;

    input.pSv = [](Real) { return 1.0e3; };
    input.pAt = atrial_pressure;
    input.u = periodic_activation;

    {
      using PassiveEnergy = std::decay_t<decltype(input.passiveEnergy)>;

      typename PassiveEnergy::Parameters hp;
      hp.mu1 = 0.0;
      hp.mu2 = 0.0;
      hp.C0 = 1.9e3;
      hp.C1 = 1.1e-1;
      hp.C2 = 1.9e3;
      hp.C3 = 1.1e-1;

      input.passiveEnergy = PassiveEnergy(hp);
    }

    return input;
  }

  /**
   * Constant-resistance RCR update:
   *
   *   C (pc^{n+1} - pc^n) / dt + (pc^{n+1} - pd) / Rd = Q,
   *
   * hence
   *
   *   pc^{n+1}
   *     =
   *   ((C/dt) pc^n + Q + pd/Rd) / (C/dt + 1/Rd),
   *
   * and
   *
   *   pout^{n+1} = pc^{n+1} + Rp Q.
   *
   * Sign convention:
   *
   *   Q = int_Gamma u · n
   *
   * so Q > 0 means flow leaving the 3D coronary domain.
   */
  void CoupledLV0DCoronary3D::updateRCR(RCR& bc, Real Q, Real dt)
  {
    const Real a = bc.C / dt;

    bc.pc =
      (a * bc.pc + Q + bc.pd / bc.Rd)
      / (a + 1.0 / bc.Rd);

    bc.pout = bc.pc + bc.Rp * Q;
  }

  /**
   * Nonlinear rheological RCR update.
   *
   * This replaces the constant distal/proximal resistance laws by nonlinear
   * pipe-flow laws:
   *
   *   Q = Q_pipe(dp; L, R).
   *
   * The distal capacitor equation is solved implicitly:
   *
   *   C (pc^{n+1} - pc^n) / dt
   *   + Q_d(pc^{n+1} - pd)
   *   = Q_out.
   *
   * The proximal pressure is then recovered from:
   *
   *   Q_out = Q_p(pout - pc).
   *
   * This is still an outlet model. It is not a monolithic 0D--3D nonlinear
   * coupling because the 3D solve uses the previous outlet pressures and the
   * RCR states are updated after the 3D solve.
   */
  void CoupledLV0DCoronary3D::updateRCRNonNew(
      const Model& model, RCR& bc, Real Q, Real dt)
  {
    (void) model;

    const auto& s = m_model.getState();

    const Real cap = bc.C / dt;
    const Real pcOld = bc.pc;

    /*
     * Equivalent unresolved segment dimensions.
     *
     * These are modelling parameters for the lumped proximal and distal
     * rheological closures, not geometric data extracted from the mesh.
     */
    const Real radiusP = 0.004;
    const Real lengthP = 0.015;

    const Real radiusD = 0.0004;
    const Real lengthD = 0.002;

    /*
     * Fully developed generalized-Newtonian pipe-flow law.
     *
     * Input:
     *   dp     = pressure drop across the unresolved segment
     *   L      = segment length
     *   radius = segment radius
     *
     * Output:
     *   Q(dp)  = volumetric flow rate
     *   dQ/ddp = differential conductance
     *
     * Constitutive law:
     *
     *   mu(gamma)
     *     =
     *   mu_inf
     *   + (mu_0 - mu_inf)
     *     [1 + (lambda gamma)^a]^((n - 1) / a).
     *
     * Wall shear stress:
     *
     *   tau_w = R |dp| / (2 L).
     *
     * Wall shear rate:
     *
     *   tau_w = mu(gamma_w) gamma_w.
     *
     * WRMS/Rabinowitsch integral:
     *
     *   I =
     *   int_0^{gamma_w}
     *     gamma^3 mu(gamma)^2 d tau/d gamma d gamma,
     *
     * where
     *
     *   tau(gamma) = mu(gamma) gamma.
     *
     * Flow rate:
     *
     *   |Q| = pi R^3 I / tau_w^3,
     *
     * with sign(Q) = sign(dp).
     */
    auto flowLaw = [](Real dp, Real L, Real radius) -> std::pair<Real, Real>
    {
      const Real mu0    = 0.04868;
      const Real muInf  = 0.003605;
      const Real lambda = 3.39;
      const Real n      = 0.198;
      const Real yasuda = 1.235;
      const Real delta  = mu0 - muInf;

      const Real sgn = (dp >= 0.0) ? 1.0 : -1.0;
      const Real adp = std::abs(dp);

      /*
       * Low-pressure-drop fallback.
       *
       * At vanishing shear, Carreau--Yasuda gives mu -> mu0, so the
       * pipe law reduces to Poiseuille with viscosity mu0:
       *
       *   R0 = 8 mu0 L / (pi R^4).
       */
      const Real R0 =
        8.0 * mu0 * L /
        (std::numbers::pi_v<Real> * std::pow(radius, 4.0));

      if (adp < 1e-12)
        return {dp / R0, 1.0 / R0};

      const Real tauW = radius * adp / (2.0 * L);

      auto mu = [&](Real g) -> Real
      {
        return muInf
          + delta * std::pow(
              1.0 + std::pow(lambda * g, yasuda),
              (n - 1.0) / yasuda);
      };

      auto dmu = [&](Real g) -> Real
      {
        const Real base = 1.0 + std::pow(lambda * g, yasuda);

        return delta * (n - 1.0)
          * std::pow(base, (n - 1.0 - yasuda) / yasuda)
          * std::pow(lambda, yasuda)
          * std::pow(g, yasuda - 1.0);
      };

      /*
       * Solve:
       *
       *   F(gamma) = gamma mu(gamma) - tau_w = 0.
       *
       * Exact derivative:
       *
       *   F'(gamma) = mu(gamma) + gamma mu'(gamma).
       */
      auto tauMinusTauW = [&](Real g) -> std::pair<Real, Real>
      {
        const Real m  = mu(g);
        const Real dm = dmu(g);

        return {g * m - tauW, m + g * dm};
      };

      Math::RootFinding::NewtonRaphson<Real> rootFinder(
        1e-12, 1e-10, 1e-12, 50);

      /*
       * Safeguarded Newton needs a bracket.
       *
       * At gamma ~= 0:
       *
       *   gamma mu(gamma) - tau_w < 0.
       *
       * For large gamma, the left-hand side becomes positive. We expand the
       * upper bound until the sign changes.
       */
      Real gHi = std::max<Real>(tauW / muInf, 1e-8);

      for (int k = 0; k < 100 && tauMinusTauW(gHi).first < 0.0; ++k)
        gHi *= 2.0;

      if (tauMinusTauW(gHi).first < 0.0)
      {
        std::cerr << "Warning: failed to bracket wall shear rate. "
                  << "Using Poiseuille fallback.\n";
        return {dp / R0, 1.0 / R0};
      }

      const auto gammaRoot =
        rootFinder.solve(tauMinusTauW, 0.5 * gHi, 1e-12, gHi);

      if (!gammaRoot)
      {
        std::cerr << "Warning: failed to solve wall shear rate. "
                  << "Using Poiseuille fallback.\n";
        return {dp / R0, 1.0 / R0};
      }

      const Real gammaW = *gammaRoot;

      /*
       * Rheological integral:
       *
       *   I =
       *   int_0^{gamma_w}
       *     gamma^3 mu(gamma)^2
       *     [mu(gamma) + gamma mu'(gamma)]
       *   d gamma.
       */
      auto integrand = [&](Real g) -> Real
      {
        if (g <= 0.0)
          return 0.0;

        const Real m     = mu(g);
        const Real dm    = dmu(g);
        const Real dtau  = m + g * dm;

        return std::pow(g, 3.0) * m * m * dtau;
      };

      Math::RungeKutta::RK4 integrator;

      constexpr int steps = 100;
      const Real h = gammaW / static_cast<Real>(steps);

      Real I = 0.0;

      auto rhs = [&](Real g, Real y) -> Real
      {
        (void) y;
        return integrand(g);
      };

      /*
       * Integrate dI/dgamma = integrand(gamma), I(0) = 0.
       */
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

      const Real qAbs =
        std::numbers::pi_v<Real> * std::pow(radius, 3.0) * I /
        std::pow(tauW, 3.0);

      /*
       * Differential conductance.
       *
       * From the Rabinowitsch relation:
       *
       *   d|Q| / d|dp|
       *     =
       *   (pi R^3 gamma_w - 3 |Q|) / |dp|.
       *
       * Since Q(dp) is odd, dQ/ddp is positive for both signs of dp.
       */
      const Real dqAbs =
        (std::numbers::pi_v<Real> * std::pow(radius, 3.0) * gammaW
         - 3.0 * qAbs) / adp;

      if (!std::isfinite(qAbs) || !std::isfinite(dqAbs) || dqAbs <= 0.0)
      {
        std::cerr << "Warning: invalid WRMS flow derivative. "
                  << "Using Poiseuille fallback.\n";
        return {dp / R0, 1.0 / R0};
      }

      return {sgn * qAbs, dqAbs};
    };

    /*
     * Invert the nonlinear pipe law.
     *
     * Given targetQ, find dp such that:
     *
     *   Q_pipe(dp) = targetQ.
     */
    auto solvePressureDropForFlow =
      [&](Real targetQ, Real L, Real radius, Real guess) -> Real
      {
        if (std::abs(targetQ) < 1e-16)
          return 0.0;

        const Real sgn  = (targetQ >= 0.0) ? 1.0 : -1.0;
        const Real qAbs = std::abs(targetQ);

        /*
         * Solve in x = |dp| >= 0:
         *
         *   F(x) = sign(Q) Q_pipe(sign(Q) x) - |Q| = 0.
         *
         * Its derivative is dQ/ddp because Q_pipe is odd.
         */
        auto F = [&](Real x) -> std::pair<Real, Real>
        {
          const auto [q, dq] = flowLaw(sgn * x, L, radius);
          return {sgn * q - qAbs, dq};
        };

        Real hi = std::max<Real>(std::abs(guess), 1.0);

        for (int k = 0; k < 100 && F(hi).first < 0.0; ++k)
          hi *= 2.0;

        if (F(hi).first < 0.0)
        {
          std::cerr << "Warning: failed to bracket pressure drop for targetQ = "
                    << targetQ << ". Returning last upper bound.\n";
          return sgn * hi;
        }

        Math::RootFinding::NewtonRaphson<Real> solver(
          1e-10, 1e-9, 1e-12, 50);

        const auto root =
          solver.solve(F, std::min(std::abs(guess), hi), 0.0, hi);

        if (!root)
        {
          std::cerr << "Warning: failed to invert flow law for targetQ = "
                    << targetQ << ". Returning bracket upper bound.\n";
          return sgn * hi;
        }

        return sgn * (*root);
      };

    /*
     * Distal nonlinear capacitor equation:
     *
     *   F(pc)
     *     =
     *   C/dt (pc - pcOld)
     *   + Q_d(pc - pd)
     *   - Q_out
     *   = 0.
     *
     * Exact derivative:
     *
     *   F'(pc) = C/dt + dQ_d/d(pc - pd).
     */
    auto distalResidual = [&](Real pc) -> std::pair<Real, Real>
    {
      const Real x = pc - s.pv;
      const auto [qd, dqd] = flowLaw(x, lengthD, radiusD);

      const Real f  = cap * (pc - pcOld) + qd - Q;
      const Real df = cap + dqd;

      return {f, df};
    };

    Math::RootFinding::NewtonRaphson<Real> solver(
      1e-10, 1e-9, 1e-12, 50);

    /*
     * Bracket pc^{n+1}.
     *
     * The pressure can move above or below pd depending on the sign of Q.
     * Therefore the bracket is expanded symmetrically around pcOld and pd.
     */
    Real span = std::max<Real>(std::abs(Q) / cap + 1000.0, 1000.0);

    Real lo = std::min(pcOld, s.pv) - span;
    Real hi = std::max(pcOld, s.pv) + span;

    for (int k = 0;
         k < 100 && distalResidual(lo).first * distalResidual(hi).first > 0.0;
         ++k)
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
      const auto pcNew =
        solver.solve(distalResidual, pcOld, lo, hi);

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

    /*
     * Proximal nonlinear pressure recovery:
     *
     *   Q_out = Q_p(pout - pc).
     *
     * Recover dpP = pout - pc from the nonlinear pipe law, then set:
     *
     *   pout = pc + dpP.
     */
    const Real oldGuess = bc.pout - bc.pc;
    const Real dpP =
      solvePressureDropForFlow(Q, lengthP, radiusP, oldGuess);

    bc.pout = bc.pc + dpP;
  }

  CoupledLV0DCoronary3D::Real
  CoupledLV0DCoronary3D::periodic_activation(Real t)
  {
    const Real T = 0.85;
    const Real tau = t - T * std::floor(t / T);

    if (tau < 0.13)  return 0.0;
    if (tau < 0.141) return 35.0 * ((tau - 0.13) / 0.011);
    if (tau < 0.281) return 35.0;
    if (tau < 0.361) return 35.0 - 47.0 * ((tau - 0.281) / 0.08);
    if (tau < 0.45)  return -12.0;

    return 0.0;
  }

  CoupledLV0DCoronary3D::Real
  CoupledLV0DCoronary3D::atrial_pressure(Real t)
  {
    const Real T = 0.85;
    const Real tau = t - T * std::floor(t / T);

    const Real minValue = 500.0;
    const Real maxValue = 1000.0;
    const Real secondThreshold = 1250.0;

    const Real t1 = 0.02;
    const Real t2 = 0.15;
    const Real t3 = 0.17;
    const Real t4 = 0.56;
    const Real t5 = 0.62;
    const Real t6 = 0.85;

    Real alpha = 0.0;
    Real value = minValue;

    if (tau < t1)
    {
      alpha = -(tau - t1) / t1;
      value = alpha * minValue + (1.0 - alpha) * maxValue;
    }
    else if (tau < t2)
    {
      value = maxValue;
    }
    else if (tau < t3)
    {
      alpha = -(tau - t3) / (t3 - t2);
      value = alpha * maxValue + (1.0 - alpha) * minValue;
    }
    else if (tau < t4)
    {
      alpha = -(tau - t4) / (t4 - t3);
      value = alpha * minValue + (1.0 - alpha) * secondThreshold;
    }
    else if (tau < t5)
    {
      value = secondThreshold;
    }
    else if (tau < t6)
    {
      alpha = -(tau - t6) / (t6 - t5);
      value = alpha * secondThreshold + (1.0 - alpha) * minValue;
    }

    return value;
  }

  CoupledLV0DCoronary3D& CoupledLV0DCoronary3D::initialize()
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
    m_model.setMaxIterations(200)
           .setAbsoluteTolerance(1e-8)
           .setRelativeTolerance(1e-8)
           .setStepTolerance(1e-10)
           .setDampingFactor(1.0);

    Model::State s0;
    s0.t = 0.0;
    s0.y = 0.0;
    s0.v = 0.0;
    s0.pv = m_input.pAt(0.0) - 100.0;
    s0.par = 11000.0;
    s0.pd = 10000.0;
    s0.ec = m_input.initFibDef;
    s0.gamma = std::sqrt(std::max<Real>(m_input.initActiveStiffness, 0.0));
    s0.beta = (s0.gamma > 0.0) ? (m_input.initActiveStress / s0.gamma) : 0.0;
    s0.kc = s0.gamma * s0.gamma;
    s0.tauc = s0.gamma * s0.beta;

    m_model.initialize(s0);
  }

  void CoupledLV0DCoronary3D::setupMeshAndSpaces()
  {
    Alert::Info() << "Setting up " << m_cfg.xdmfBasename << ".xdmf ..."
                  << Alert::Raise;

    m_xdmf.setMesh(m_mesh);

    m_u.setName("u");
    m_p.setName("p");
    m_mu.setName("mu_nonNew");
    m_up.setName("projected_convection");
    m_sub.setName("subscale");
    m_tau.setName("tau");

    m_uOld = Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};
    m_pOld = 0.0;

    m_tauOld = 0.0;
    m_subOld = Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};

    m_one = 1.0;

    m_xdmf.add("velocity", m_u.getSolution());
    m_xdmf.add("pressure", m_p.getSolution());
    m_xdmf.add("viscosity", m_mu.getSolution());
    m_xdmf.add("subscale", m_sub.getSolution());
    m_xdmf.add("tau", m_tau.getSolution());

    m_wk.clear();
    for (const Attribute tag : m_cfg.outlets)
      m_wk.emplace(tag, m_cfg.defaultRCR);
  }

  void CoupledLV0DCoronary3D::setupDiagnostics()
  {
    m_csv.open(m_cfg.csvPath);

    if (!m_csv)
      throw std::runtime_error("Failed to open coupled CSV file: " + m_cfg.csvPath);

    writeCSVHeader();
  }

  void CoupledLV0DCoronary3D::printInitialState() const
  {
    const auto& s = m_model.getState();

    std::cout << "Initial 0D state:\n"
              << "  y     = " << s.y << '\n'
              << "  v     = " << s.v << '\n'
              << "  pv    = " << s.pv << '\n'
              << "  par   = " << s.par << '\n'
              << "  pd    = " << s.pd << '\n'
              << "  ec    = " << s.ec << '\n'
              << "  gamma = " << s.gamma << '\n'
              << "  beta  = " << s.beta << '\n'
              << "  kc    = " << s.kc << '\n'
              << "  tauc  = " << s.tauc << '\n';
  }

  bool CoupledLV0DCoronary3D::advance0D()
  {
    const auto rep = m_model.step(m_cfg.dt);

    std::cout << "  0D Newton step: "
              << (rep.converged ? "converged" : "not converged")
              << ", iterations = " << rep.iterations
              << ", final residual = " << rep.finalResidual
              << ", final step norm = " << rep.finalStepNorm
              << '\n';

    return rep.converged;
  }

  void CoupledLV0DCoronary3D::solve3D()
  {
    const auto& s = m_model.getState();

    /*
     * One-way coupling:
     *
     *   p_in = p_ar(t)
     *
     * from the 0D model. This pressure is imposed as a Neumann traction at
     * the coronary inlet.
     */
    const Real pin = s.par;

    const auto n = BoundaryNormal(m_mesh);

    /*
     * Oseen-linearized convection:
     *
     *   (u · grad) u
     *   approximated by
     *   (u^n · grad) u^{n+1}
     *
     * in Rodin notation:
     *
     *   Mult(Jacobian(m_u), m_uOld)
     *
     * i.e. (grad u^{n+1}) u^n.
     */
    const auto conv_u = Mult(Jacobian(m_u), m_uOld);

    /*
     * Skew-symmetric / Temam-style correction:
     *
     *   0.5 rho (div u^n) u^{n+1}.
     */
    const auto div_u_old = Div(m_uOld);

    /*
     * Backflow coefficient on outlets:
     *
     *   beta = max(-u^n · n, 0).
     */
    const auto beta = Max(-Dot(m_uOld, n), 0.0);

    auto symU =
      0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));

    auto symV =
      0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));

    auto symUOld =
      0.5 * (Jacobian(m_uOld) + Transpose(Jacobian(m_uOld)));

    /*
     * Carreau--Yasuda viscosity:
     *
     *   mu(gamma)
     *     =
     *   mu_inf
     *   + (mu_0 - mu_inf)
     *     [1 + (lambda gamma)^a]^((n - 1) / a),
     *
     * with:
     *
     *   gamma = sqrt(2 D(u^n) : D(u^n)).
     *
     * The viscosity is evaluated explicitly from u^n, making the 3D solve
     * linear in (u^{n+1}, p^{n+1}).
     */
    const Real gamma_min = 1.0e-3;
    const Real mu_0 = 0.04868;
    const Real mu_inf = 0.003605;
    const Real lambda = 3.39;
    const Real n_cy = 0.198;
    const Real a = 1.235;

    /*
     * VMS stabilization constants.
     *
     * tau_K is computed as:
     *
     *   tau_K =
     *   vmsScale /
     *   sqrt(
     *     (2 rho / dt)^2
     *     + (c2 rho |u^n| / h_K)^2
     *     + (c1 mu / h_K^2)^2
     *   ).
     */
    const Real c1 = 4.0;
    const Real c2 = 2.0;
    const Real vmsScale = 1.0;

    RealFunction muNonNew = [=, this](const Point& p) -> Real
    {
      const auto S = symUOld.getValue(p);

      const Real gamma =
        std::max(gamma_min, std::sqrt(2.0 * dot(S, S)));

      return mu_inf
        + (mu_0 - mu_inf)
        * Math::pow(
            1.0 + Math::pow(lambda * gamma, a),
            (n_cy - 1.0) / a);
    };

    Alert::Info() << "Projecting non-Newtonian viscosity ..." << Alert::Raise;

    /*
     * L2 projection of the explicit viscosity into the P0/Pcell space m_muh:
     *
     *   int mu_h w = int mu(u^n) w.
     */
    m_muProjection =
        Integral(m_mu, m_w)
      - Integral(muNonNew, m_w);

    m_muProjection.assemble();
    m_muProjectionSolver.solve();

    /*
     * Convective residual target:
     *
     *   a^n = (grad u^n) u^n.
     */
    const auto convectionTarget = Mult(Jacobian(m_uOld), m_uOld);

    RealFunction tau = [=, this](const Point& p) -> Real {
        const auto conv = convectionTarget.getValue(p);
        const auto proj = m_up.getSolution().getValue(p);
        const auto old  = m_subOld.getValue(p);
        const auto uOld = m_uOld.getValue(p);

        const Real mu = m_mu.getSolution().getValue(p);

        const Real hK =
          std::pow(
            p.getPolytope().getMeasure(),
            1.0 / p.getPolytope().getDimension());

        const Real speed = std::sqrt(dot(uOld, uOld));

        const Real invTau =
          std::sqrt(
              Math::pow2(2.0 * m_cfg.rho / m_cfg.dt)
            + Math::pow2(c2 * m_cfg.rho * speed / hK)
            + Math::pow2(c1 * mu / (hK * hK)));

         return vmsScale / invTau;
    }
    /*
     * L2 projection:
     *
     *   int u_proj · v = int a^n · v.
     *
     * Thus:
     *
     *   u_proj = Pi_h[(grad u^n) u^n].
     */
    m_l2ConvU =
        Integral(m_up, m_vp)
      - Integral(convectionTarget, m_vp);

    m_l2ConvU.assemble();
    m_l2ConvUSolver.solve();

    /*
     * Dynamic unresolved convective subscale:
     *
     *   u'^{n+1}
     *     =
     *   tau_K rho (a^n - Pi_h a^n)
     *   + tau_K (rho / dt) u'^n.
     *
     * This is explicit in the resolved velocity because a^n and Pi_h a^n
     * are both computed from u^n.
     */

    m_tauProjection =
        Integral(m_tau, m_t)
      - Integral(tau, m_t);

    m_tauProjection.assemble();
    m_tauProjectionSolver.solve();

    m_l2ConvU =
        Integral(m_up, m_vp)
      - Integral(convectionTarget, m_vp);

    m_l2ConvU.assemble();
    m_l2ConvUSolver.solve();

    auto subUpdate = VectorFunction(
        m_mesh.getSpaceDimension(),
        [=, this](const Point& p) -> Math::SpatialVector<Real>
        {
          const auto conv = convectionTarget.getValue(p);
          const auto proj = m_up.getSolution().getValue(p);
          const auto old  = m_subOld.getValue(p);
          const auto uOld = m_uOld.getValue(p);

          const Real mu = m_mu.getSolution().getValue(p);

          const Real tau = m_tau.getSolution().getValue(p);

          Math::SpatialVector<Real> out(m_mesh.getSpaceDimension());

          for (size_t c = 0; c < out.size(); ++c)
          {
            out[c] =
                tau * m_cfg.rho * (conv[c] - proj[c])
              + tau * (m_cfg.rho / m_cfg.dt) * old[c];
          }

          return out;
        });

    /*
     * L2 projection of the dynamic subscale into m_sub:
     *
     *   int sub_h · v = int subUpdate · v.
     */
    m_subProjection =
        Integral(m_sub, m_vp)
      - Integral(subUpdate, m_vp);

    m_subProjection.assemble();
    m_subProjectionSolver.solve();

    /*
     * Semi-implicit weak form:
     *
     * Find (u^{n+1}, p^{n+1}) such that, for all (v, q),
     *
     *   rho/dt (u^{n+1}, v)
     * - rho/dt (u^n, v)
     * + rho ((grad u^{n+1}) u^n, v)
     * + 0.5 rho ((div u^n) u^{n+1}, v)
     * + 2 (mu^n D(u^{n+1}), D(v))
     * - (p^{n+1}, div v)
     * + (div u^{n+1}, q)
     * + eps (p^{n+1}, q)
     *
     * plus VMS stabilization, pressure tractions, backflow stabilization,
     * and wall no-slip conditions.
     */
    m_flow =
          (m_cfg.rho / m_cfg.dt) * Integral(m_u, m_v)
        - (m_cfg.rho / m_cfg.dt) * Integral(m_uOld, m_v)

        + m_cfg.rho * Integral(Dot(conv_u, m_v))
        + 0.5 * m_cfg.rho * Integral(div_u_old * Dot(m_u, m_v))

        + 2.0 * Integral(m_mu.getSolution() * symU, symV)

        - Integral(m_p, Div(m_v))
        + Integral(Div(m_u), m_q)
        + m_cfg.eps * Integral(m_p, m_q)

        /*
         * VMS bilinear contribution:
         *
         *   int_K tau_K rho^2
         *     ((grad u^{n+1}) u^n)
         *     ·
         *     ((grad v) u^n).
         */

        + VMSConvectionBilinearIntegrator(
            m_u, m_v, m_uOld, m_tau.getSolution())

        /*
         * VMS linear contribution subtracted from the residual:
         *
         *   - int_K rho
         *       (tau_K rho u_proj + u'^{n+1})
         *       ·
         *       ((grad v) u^n).
         */
        - VMSConvectionLinearIntegrator(
            m_v, m_sub.getSolution(), m_uOld, m_up.getSolution(), m_tau.getSolution())

        /*
         * Inlet pressure traction.
         */
        + BoundaryIntegral(pin * Dot(m_v, n)).over(m_cfg.inlet)

        /*
         * Outlet pressure tractions.
         */
        + BoundaryIntegral(m_wk.at(4).pout * Dot(m_v, n)).over(4)
        + BoundaryIntegral(m_wk.at(5).pout * Dot(m_v, n)).over(5)
        + BoundaryIntegral(m_wk.at(6).pout * Dot(m_v, n)).over(6)
        + BoundaryIntegral(m_wk.at(7).pout * Dot(m_v, n)).over(7)
        + BoundaryIntegral(m_wk.at(8).pout * Dot(m_v, n)).over(8)
        + BoundaryIntegral(m_wk.at(9).pout * Dot(m_v, n)).over(9)

        /*
         * Outlet backflow stabilization:
         *
         *   - 0.5 rho int_Gamma max(-u^n · n, 0) u^{n+1} · v.
         */
        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(m_u, m_v)).over(4)
        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(m_u, m_v)).over(5)
        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(m_u, m_v)).over(6)
        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(m_u, m_v)).over(7)
        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(m_u, m_v)).over(8)
        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(m_u, m_v)).over(9)

        /*
         * No-slip wall.
         */
        + DirichletBC(m_u, Zero(m_mesh.getSpaceDimension())).on(m_cfg.wall);

    Alert::Info() << "Assembling 3D time step ..." << Alert::Raise;

    m_flow.assemble();

    if (!m_flowFieldSplitsSet)
    {
      m_flow.setFieldSplits();
      m_flowFieldSplitsSet = true;
    }

    Alert::Info() << "Solving 3D time step ..." << Alert::Raise;
    m_flowSolver.solve();
  }

  void CoupledLV0DCoronary3D::computeFluxesAndUpdateRCR()
  {
    const auto n = BoundaryNormal(m_mesh);
    const auto& s = m_model.getState();

    /*
     * Inlet flux:
     *
     *   Q_in = int_Gamma_in u · n.
     *
     * With outward normals, an entering inflow usually gives Q_in < 0.
     */
    m_flux =
      BoundaryIntegral(Dot(m_u.getSolution(), n), m_qFlux).over(m_cfg.inlet);

    m_flux.assemble();
    m_stepData.qIn = m_flux(m_one);

    m_stepData.qOut.clear();
    m_stepData.qOutSum = 0.0;

    /*
     * Outlet fluxes:
     *
     *   Q_out,k = int_Gamma_k u · n.
     *
     * With outward normals, Q_out,k > 0 means blood leaves the 3D domain.
     */
    for (const Attribute tag : m_cfg.outlets)
    {
      m_flux =
        BoundaryIntegral(Dot(m_u.getSolution(), n), m_qFlux).over(tag);

      m_flux.assemble();

      const Real qOut = m_flux(m_one);

      m_stepData.qOut[tag] = qOut;
      m_stepData.qOutSum += qOut;

      updateRCRNonNew(m_model, m_wk[tag], qOut, m_cfg.dt);
    }

    /*
     * For an exactly conservative incompressible solve:
     *
     *   Q_in + sum_k Q_out,k ~= 0.
     */
    m_stepData.flowBalance = m_stepData.qIn + m_stepData.qOutSum;
    m_stepData.t = s.t;
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
    /*
     * The 0D model has already advanced to t^{n+1}, while the current 3D
     * solution corresponds to the just-computed time slab. The existing driver
     * writes at t - dt to preserve the original output convention.
     */
    m_xdmf.write(m_model.getState().t - m_cfg.dt).flush();
  }

  CoupledLV0DCoronary3D::StepData
  CoupledLV0DCoronary3D::collectStepData() const
  {
    StepData d;

    const auto& s = m_model.getState();

    d.t = s.t;
    d.pat = m_input.pAt(s.t);
    d.psv = m_input.pSv(s.t);

    d.y = s.y;
    d.v = s.v;
    d.radius = m_input.R0 + s.y;

    d.lvVolume =
      (4.0 / 3.0) * std::numbers::pi_v<Real>
      * d.radius * d.radius * d.radius;

    d.lvFlow =
      4.0 * std::numbers::pi_v<Real>
      * d.radius * d.radius * s.v;

    d.pv = s.pv;
    d.par = s.par;
    d.pd = s.pd;

    d.ec = s.ec;
    d.gamma = s.gamma;
    d.beta = s.beta;
    d.kc = s.kc;
    d.tauc = s.tauc;

    d.qIn = m_stepData.qIn;
    d.qOutSum = m_stepData.qOutSum;
    d.flowBalance = m_stepData.flowBalance;
    d.qOut = m_stepData.qOut;

    for (const auto& [tag, bc] : m_wk)
    {
      d.pc[tag] = bc.pc;
      d.pOut[tag] = bc.pout;
    }

    return d;
  }

  void CoupledLV0DCoronary3D::writeCSVHeader()
  {
    m_csv
      << "t,"
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
      << "kc,"
      << "tauc\n";

    m_csv.flush();
  }

  void CoupledLV0DCoronary3D::writeCSVRow()
  {
    const StepData d = collectStepData();

    auto get = [](const std::map<Attribute, Real>& m, Attribute a) -> Real
    {
      const auto it = m.find(a);
      return (it == m.end()) ? 0.0 : it->second;
    };

    m_csv
      << d.t << ','
      << d.pat << ','
      << d.psv << ','
      << d.y << ','
      << d.v << ','
      << d.radius << ','
      << d.lvVolume << ','
      << d.pv << ','
      << d.par << ','
      << d.pd << ','
      << d.lvFlow << ','
      << d.qIn << ','
      << get(d.qOut, 4) << ','
      << get(d.qOut, 5) << ','
      << get(d.qOut, 6) << ','
      << get(d.qOut, 7) << ','
      << get(d.qOut, 8) << ','
      << get(d.qOut, 9) << ','
      << d.qOutSum << ','
      << get(d.pc, 4) << ','
      << get(d.pc, 5) << ','
      << get(d.pc, 6) << ','
      << get(d.pc, 7) << ','
      << get(d.pc, 8) << ','
      << get(d.pc, 9) << ','
      << get(d.pOut, 4) << ','
      << get(d.pOut, 5) << ','
      << get(d.pOut, 6) << ','
      << get(d.pOut, 7) << ','
      << get(d.pOut, 8) << ','
      << get(d.pOut, 9) << ','
      << d.flowBalance << ','
      << d.ec << ','
      << d.gamma << ','
      << d.beta << ','
      << d.kc << ','
      << d.tauc << '\n';

    m_csv.flush();
  }

  int CoupledLV0DCoronary3D::run()
  {
    if (!m_initialized)
      initialize();

    for (int i = 0; i < static_cast<int>(m_cfg.nsteps); ++i)
    {
      std::cout << "Step " << i << ": t = " << m_model.getState().t << "\n";

      if (!advance0D())
      {
        std::cerr << "0D solver failed to converge at step "
                  << i << ", t = " << m_model.getState().t << "\n";
        return 1;
      }

      solve3D();
      computeFluxesAndUpdateRCR();
      writeCSVRow();
      updateHistory();
      writeOutputs();
    }

    m_xdmf.close();
    return 0;
  }
}
