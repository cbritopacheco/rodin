#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
#include <type_traits>

#include <boost/mpi/collectives.hpp>

#include <petscvec.h>

#include <Rodin/Alert.h>
#include <Rodin/Solver.h>

#include "CoupledLV0DCoronary3D.h"
#include "Rodin/Alert/Exception.h"
#include "Rodin/Math/RungeKutta/RK4.h"
#include "Rodin/Variational/ForwardDecls.h"
#include "Rodin/Math/RootFinding/NewtonRaphson.h"

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

    struct MaxReal
    {
      Rodin::Real operator()(Rodin::Real lhs, Rodin::Real rhs) const
      {
        return std::max(lhs, rhs);
      }
    };

    Rodin::Real secondsSince(Clock::time_point start)
    {
      return std::chrono::duration<Rodin::Real>(Clock::now() - start).count();
    }
  }

  CoupledLV0DCoronary3D::CoupledLV0DCoronary3D(const Context::MPI& context)
    : CoupledLV0DCoronary3D(context, Config{})
  {
  }

  CoupledLV0DCoronary3D::CoupledLV0DCoronary3D(
      const Context::MPI& context, const Config& cfg)
    : m_cfg(cfg),
      m_input(makeInput()),
      m_model(m_input),
      m_mesh(makeMesh(context, m_cfg)),
      m_xdmf(context.getCommunicator(), m_cfg.xdmfBasename),
      m_uh(std::integral_constant<size_t, 2>{}, m_mesh, m_mesh.getSpaceDimension()),
      m_ph(std::integral_constant<size_t, 1>{}, m_mesh),
      m_u(m_uh),
      m_p(m_ph),
      m_v(m_uh),
      m_q(m_ph),
      m_uOld(m_uh),
      m_pOld(m_ph),
      m_one(m_ph),
      m_qFlux(m_ph),
      m_flux(m_qFlux),
      m_flow(m_u, m_p, m_v, m_q),
      m_flowKSP(m_flow),
      m_flowSolver(m_flowKSP)
  {
    auto cellCount = m_mesh.getCellCount();
    auto vertexCount = m_mesh.getVertexCount();
    auto spaceDim = m_mesh.getSpaceDimension();
    auto velocityDOFs = m_uh.getSize();
    auto pressureDOFs = m_ph.getSize();

    if (isRoot())
    {
      Alert::Info() << "Number of elements in mesh: " << cellCount << '\n'
                    << "Number of vertices in mesh: " << vertexCount << '\n'
                    << "Mesh space dimension: " << spaceDim << '\n'
                    << Alert::Raise;

      Alert::Info() << "Velocity space has " << velocityDOFs << " DOFs.\n"
                    << "Pressure space has " << pressureDOFs << " DOFs.\n"
                    << Alert::Raise;
    }
  }

  CoupledLV0DCoronary3D::~CoupledLV0DCoronary3D() = default;

  CoupledLV0DCoronary3D::MeshType
  CoupledLV0DCoronary3D::makeMesh(const Context::MPI& context, const Config& cfg)
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

      Alert::Info() << "Partitioning coronary mesh over "
                    << comm.size() << " MPI ranks ..."
                    << Alert::Raise;

      Geometry::BalancedCompactPartitioner partitioner(mesh);
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

    mesh.save(
        cfg.xdmfBasename + "_partitioned"
        + std::to_string(mesh.getContext().getCommunicator().rank()) + ".mesh",
        IO::FileFormat::MEDIT);

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

  void CoupledLV0DCoronary3D::updateRCR(RCR& bc, Real Q, Real dt)
  {
    const Real a = bc.C / dt;

    bc.pc =
      (a * bc.pc + Q + bc.pd / bc.Rd)
      / (a + 1.0 / bc.Rd);

    bc.pout = bc.pc + bc.Rp * Q;
  }

  void CoupledLV0DCoronary3D::updateRCRNonNew(
      const Model& model, RCR& bc, Real Q, Real dt)
  {
    const auto& s = model.getState();

    const Real cap = bc.C / dt;
    const Real pcOld = bc.pc;

    const Real radiusP = 0.004;
    const Real lengthP = 0.015;

    const Real radiusD = 0.0004;
    const Real lengthD = 0.002;

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

      auto tauMinusTauW = [&](Real g) -> std::pair<Real, Real>
      {
        const Real m  = mu(g);
        const Real dm = dmu(g);
        return {g * m - tauW, m + g * dm};
      };

      Math::RootFinding::NewtonRaphson<Real> rootFinder(
        1e-12, 1e-10, 1e-12, 50);

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

    auto solvePressureDropForFlow =
      [&](Real targetQ, Real L, Real radius, Real guess) -> Real
      {
        if (std::abs(targetQ) < 1e-16)
          return 0.0;

        const Real sgn  = (targetQ >= 0.0) ? 1.0 : -1.0;
        const Real qAbs = std::abs(targetQ);

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
    if (isRoot())
      Alert::Info() << "Setting up " << m_cfg.xdmfBasename << ".xdmf ..."
                    << Alert::Raise;

    m_xdmf.setMesh(m_mesh);

    m_u.setName("u");
    m_p.setName("p");

    m_uOld = Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};
    m_pOld = 0.0;
    m_one = 1.0;

    m_xdmf.add("velocity", m_u.getSolution());
    m_xdmf.add("pressure", m_p.getSolution());

    m_wk.clear();
    for (const Attribute tag : m_cfg.outlets)
      m_wk.emplace(tag, m_cfg.defaultRCR);

    m_flowSolver
      .setTolerances(1e-10, 1e-8, 1e-10, 50, 10000)
      .setStateUpdate([this](const PETSc::Math::Vector& x)
      {
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

  bool CoupledLV0DCoronary3D::isRoot() const
  {
    return m_mesh.getContext().getCommunicator().rank() == RootRank;
  }

  bool CoupledLV0DCoronary3D::advance0D()
  {
    if (isRoot())
      Alert::Info() << "\t[0D] Advancing LV model ..." << Alert::Raise;

    const auto rep = m_model.step(m_cfg.dt);

    if (isRoot())
    {
      Alert::Info()
        << "\t[0D] Newton: "
        << (rep.converged ? "converged" : "NOT converged")
        << "  iter = " << rep.iterations
        << "  |F| = " << rep.finalResidual
        << "  |dx| = " << rep.finalStepNorm
        << Alert::Raise;
    }

    if (isRoot() && rep.converged)
    {
      const auto& s = m_model.getState();
      Alert::Info()
        << "\t\t pv  = " << s.pv  << " Pa"
        << "  par = " << s.par << " Pa"
        << "  pd  = " << s.pd  << " Pa"
        << Alert::Raise;
    }

    return rep.converged;
  }

  bool CoupledLV0DCoronary3D::solve3D()
  {
    const auto setup3DStart = Clock::now();

    const auto& s = m_model.getState();
    const Real pin = s.par;

    const auto normal = BoundaryNormal(m_mesh);

    const auto& uState = m_u.getSolution();
    const auto& pState = m_p.getSolution();

    const auto gradDU = Jacobian(m_u);
    const auto gradU  = Jacobian(uState);

    const auto newtonConvection =
      Mult(gradDU, uState) + Mult(gradU, m_u);

    const auto stateConvection = Mult(gradU, uState);

    const auto divDU = Div(m_u);
    const auto divU  = Div(uState);

    const auto temamJacobian1 =
      Dot(divDU * uState, m_v);

    const auto temamJacobian2 =
      divU * Dot(m_u, m_v);

    const auto temamResidual =
      divU * Dot(uState, m_v);

    const auto beta = Max(-Dot(m_uOld, normal), 0.0);

    const auto symDU =
      0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));

    const auto symV =
      0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));

    const auto symU =
      0.5 * (Jacobian(uState) + Transpose(Jacobian(uState)));

    const Real gammaReg = 1.0e-3;
    const Real mu0      = 0.04868;
    const Real muInf    = 0.003605;
    const Real lambda   = 3.39;
    const Real nCY      = 0.198;
    const Real yasuda   = 1.235;
    const Real deltaMu  = mu0 - muInf;

    const auto gamma =
      Sqrt(gammaReg * gammaReg + 2.0 * Dot(symU, symU));

    const auto carreauBase =
      1.0 + Pow(lambda * gamma, yasuda);

    const auto mu =
      muInf + deltaMu * Pow(carreauBase, (nCY - 1.0) / yasuda);

    const auto dgamma =
      2.0 * Dot(symU, symDU) / gamma;

    const auto dmu =
        deltaMu
      * (nCY - 1.0)
      * Pow(carreauBase, (nCY - 1.0 - yasuda) / yasuda)
      * std::pow(lambda, yasuda)
      * Pow(gamma, yasuda - 1.0)
      * dgamma;

    m_flow =
          /* Jacobian */
          (m_cfg.rho / m_cfg.dt) * Integral(m_u, m_v)

        + m_cfg.rho * Integral(Dot(newtonConvection, m_v))

        + 0.5 * m_cfg.rho * Integral(temamJacobian1)
        + 0.5 * m_cfg.rho * Integral(temamJacobian2)

        + 2.0 * Integral(mu * symDU, symV)
        + 2.0 * Integral(dmu * symU, symV)

        - Integral(m_p, Div(m_v))
        + Integral(Div(m_u), m_q)
        + m_cfg.eps * Integral(m_p, m_q)

        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(m_u, m_v)).over(4)
        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(m_u, m_v)).over(5)
        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(m_u, m_v)).over(6)
        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(m_u, m_v)).over(7)
        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(m_u, m_v)).over(8)
        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(m_u, m_v)).over(9)

          /* Residual */
        + (m_cfg.rho / m_cfg.dt) * Integral(uState, m_v)
        - (m_cfg.rho / m_cfg.dt) * Integral(m_uOld, m_v)

        + m_cfg.rho * Integral(Dot(stateConvection, m_v))

        + 0.5 * m_cfg.rho * Integral(temamResidual)

        + 2.0 * Integral(mu * symU, symV)

        - Integral(pState, Div(m_v))
        + Integral(Div(uState), m_q)
        + m_cfg.eps * Integral(pState, m_q)

        + BoundaryIntegral(pin * Dot(m_v, normal)).over(m_cfg.inlet)

        + BoundaryIntegral(m_wk.at(4).pout * Dot(m_v, normal)).over(4)
        + BoundaryIntegral(m_wk.at(5).pout * Dot(m_v, normal)).over(5)
        + BoundaryIntegral(m_wk.at(6).pout * Dot(m_v, normal)).over(6)
        + BoundaryIntegral(m_wk.at(7).pout * Dot(m_v, normal)).over(7)
        + BoundaryIntegral(m_wk.at(8).pout * Dot(m_v, normal)).over(8)
        + BoundaryIntegral(m_wk.at(9).pout * Dot(m_v, normal)).over(9)

        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(uState, m_v)).over(4)
        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(uState, m_v)).over(5)
        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(uState, m_v)).over(6)
        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(uState, m_v)).over(7)
        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(uState, m_v)).over(8)
        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(uState, m_v)).over(9)

        + DirichletBC(m_u, -uState).on(m_cfg.wall);

    m_stepTiming.setup3DForm = secondsSince(setup3DStart);

    if (!m_flowFieldSplitsSet)
    {
      if (isRoot())
        Alert::Info() << "\t\t[3D] Initializing Newton system ..." << Alert::Raise;

      const auto assemble3DStart = Clock::now();
      m_flow.assemble();
      m_stepTiming.assemble3D = secondsSince(assemble3DStart);

      const auto fieldSplitsStart = Clock::now();
      m_flow.setFieldSplits();
      m_stepTiming.fieldSplits = secondsSince(fieldSplitsStart);
      m_flowFieldSplitsSet = true;
    }

    if (isRoot())
      Alert::Info() << "\t\t[3D] Solving with PETSc SNES ..." << Alert::Raise;

    const auto solve3DStart = Clock::now();
    m_flowSolver.solve();
    m_stepTiming.solve3D = secondsSince(solve3DStart);

    if (isRoot())
    {
      Alert::Info()
        << "\t\t\t[SNES] "
        << (m_flowSolver.converged() ? "Converged" : "Did NOT converge")
        << "  iterations = " << m_flowSolver.getIterationNumber()
        << Alert::Raise;
    }
    return m_flowSolver.converged();
  }

  void CoupledLV0DCoronary3D::computeFluxes()
  {
    const auto normal = BoundaryNormal(m_mesh);
    const auto& s = m_model.getState();

    m_flux =
      BoundaryIntegral(Dot(m_u.getSolution(), normal), m_qFlux).over(m_cfg.inlet);

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

  void CoupledLV0DCoronary3D::updateHistory()
  {
    m_uOld.setData(m_u.getSolution().getData());
    m_pOld.setData(m_p.getSolution().getData());
  }

  void CoupledLV0DCoronary3D::writeOutputs()
  {
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
    if (!isRoot())
      return;

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
    if (!isRoot())
      return;

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

  void CoupledLV0DCoronary3D::printStepTiming(int step) const
  {
    const auto& comm = m_mesh.getContext().getCommunicator();
    auto maxTime = [&](Real local)
    {
      return boost::mpi::all_reduce(comm, local, MaxReal{});
    };

    const Real total       = maxTime(m_stepTiming.total);
    const Real advance0D   = maxTime(m_stepTiming.advance0D);
    const Real setup3DForm = maxTime(m_stepTiming.setup3DForm);
    const Real assemble3D  = maxTime(m_stepTiming.assemble3D);
    const Real fieldSplits = maxTime(m_stepTiming.fieldSplits);
    const Real solve3D     = maxTime(m_stepTiming.solve3D);
    const Real fluxes      = maxTime(m_stepTiming.fluxes);
    const Real outletRCR   = maxTime(m_stepTiming.outletRCR);
    const Real csv         = maxTime(m_stepTiming.csv);
    const Real history     = maxTime(m_stepTiming.history);
    const Real output      = maxTime(m_stepTiming.output);

    if (isRoot())
    {
      Alert::Info()
        << "\t[Timing max/rank] step = " << step
        << "  total = " << total << " s"
        << "  0D = " << advance0D << " s"
        << "  3D-form = " << setup3DForm << " s"
        << "  3D-assemble = " << assemble3D << " s"
        << "  field-splits = " << fieldSplits << " s"
        << "  3D-solve = " << solve3D << " s"
        << "  fluxes = " << fluxes << " s"
        << "  RCR = " << outletRCR << " s"
        << "  csv = " << csv << " s"
        << "  history = " << history << " s"
        << "  output = " << output << " s"
        << Alert::Raise;
    }
  }

  int CoupledLV0DCoronary3D::run()
  {
    if (!m_initialized)
      initialize();

    for (int i = 0; i < static_cast<int>(m_cfg.nsteps); ++i)
    {
      m_stepTiming = StepTiming{};
      const auto stepStart = Clock::now();
      const Real t_current = m_model.getState().t;

      if (isRoot())
      {
        Alert::Info()
          << "━━━ Step " << (i + 1) << " / " << m_cfg.nsteps
          << "  t = " << t_current << " s"
          << "  (dt = " << m_cfg.dt << " s)"
          << " ━━━"
          << Alert::Raise;
      }

      const auto advance0DStart = Clock::now();
      const bool advanced0D = advance0D();
      m_stepTiming.advance0D = secondsSince(advance0DStart);

      if (!advanced0D)
      {
        Alert::Exception()
          << "0D solver failed to converge at step "
          << (i + 1) << ", t = " << m_model.getState().t
          << Alert::Raise;
        return 1;
      }

      if (!solve3D())
      {
        Alert::Exception()
          << "3D SNES solver failed to converge at step "
          << (i + 1) << ", t = " << m_model.getState().t
          << Alert::Raise;
        return 1;
      }

      const auto fluxesStart = Clock::now();
      computeFluxes();
      m_stepTiming.fluxes = secondsSince(fluxesStart);

      const auto outletRCRStart = Clock::now();
      for (const Attribute tag : m_cfg.outlets)
        updateRCRNonNew(m_model, m_wk[tag], m_stepData.qOut.at(tag), m_cfg.dt);
      m_stepTiming.outletRCR = secondsSince(outletRCRStart);

      const auto csvStart = Clock::now();
      writeCSVRow();
      m_stepTiming.csv = secondsSince(csvStart);

      const auto historyStart = Clock::now();
      updateHistory();
      m_stepTiming.history = secondsSince(historyStart);

      const auto outputStart = Clock::now();
      writeOutputs();
      m_stepTiming.output = secondsSince(outputStart);

      m_stepTiming.total = secondsSince(stepStart);
      printStepTiming(i + 1);
    }

    m_xdmf.close();
    return 0;
  }
}
