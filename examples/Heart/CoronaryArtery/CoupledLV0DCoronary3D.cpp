#include <algorithm>
#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <type_traits>

#include <Rodin/Alert.h>
#include <Rodin/Solver.h>

#include "CoupledLV0DCoronary3D.h"
#include "VMSConvectionIntegrator.h"

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
      m_uph(std::integral_constant<size_t, 2>{}, m_mesh, m_mesh.getSpaceDimension()),
      m_u(m_uh),
      m_p(m_ph),
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
      m_muProjectionSolver(m_muProjection),
      m_l2ConvUSolver(m_l2ConvU),
      m_subProjectionSolver(m_subProjection),
      m_flowSolver(m_flow)
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

  void CoupledLV0DCoronary3D::updateRCR(RCR& bc, Real Q, Real dt)
  {
    const Real a = bc.C / dt;
    bc.pc = (a * bc.pc + Q + bc.pd / bc.Rd) / (a + 1.0 / bc.Rd);
    bc.pout = bc.pc + bc.Rp * Q;
  }

  void CoupledLV0DCoronary3D::updateRCRNonNew(
      const Model& model, RCR& bc, Real Q, Real dt)
  {
    const Real a = bc.C / dt;

    const Real radiusP = 0.004;
    const Real lengthP = 0.015;

    const Real radiusD = 0.0005;
    const Real lengthD = 0.002;

    auto resistance = [](Real pp, Real pd, Real L, Real radius) -> Real
    {
      const Real m = 0.035;
      const Real n = 0.7;

      const Real dp = pp - pd;
      const Real tauW = radius * dp / (2.0 * L);

      const Real Ip =
        std::pow(1.0 / m, 1.0 / n)
        * n / (3.0 * n + 1.0)
        * std::pow(std::abs(tauW), (3.0 * n + 1.0) / n);

      if (Ip <= 0.0 || std::abs(dp) <= 0.0)
        return 1.0e-8;

      return std::pow(dp, 4.0)
        / (8.0 * std::numbers::pi_v<Real> * std::pow(L, 3.0) * Ip);
    };

    const Real Rp = resistance(bc.pc, bc.pd, lengthP, radiusP);

    const auto& s = model.getState();
    const Real pv = s.pv;

    const Real Rd = resistance(bc.pd, pv, lengthD, radiusD);

    bc.pc = (a * bc.pc + Q + bc.pd / Rd) / (a + 1.0 / Rd);
    bc.pout = bc.pc + Rp * Q;
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

    m_uOld = Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};
    m_pOld = 0.0;

    m_subOld = Math::SpatialVector<Real>{{0.0, 0.0, 0.0}};

    m_one = 1.0;

    m_xdmf.add("velocity", m_u.getSolution());
    m_xdmf.add("pressure", m_p.getSolution());
    m_xdmf.add("viscosity", m_mu.getSolution());
    m_xdmf.add("subscale", m_sub.getSolution());

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
    const Real pin = s.par;

    const auto n = BoundaryNormal(m_mesh);

    const auto conv_u = Mult(Jacobian(m_u), m_uOld);
    const auto div_u_old = Div(m_uOld);
    const auto beta = Max(-Dot(m_uOld, n), 0.0);

    auto symU = 0.5 * (Jacobian(m_u) + Transpose(Jacobian(m_u)));
    auto symV = 0.5 * (Jacobian(m_v) + Transpose(Jacobian(m_v)));
    auto symUOld = 0.5 * (Jacobian(m_uOld) + Transpose(Jacobian(m_uOld)));

    const Real n_pl = 0.7;
    const Real m_pl = 0.035;

    const Real mu_min = m_cfg.mu;
    const Real gamma_min = 1.0e-3;
    const Real mu_max = 5.0e-2;

    const Real c1 = 4.0;
    const Real c2 = 2.0;
    const Real vmsScale = 0.05;

    RealFunction muNonNew = [=, this](const Point& p) -> Real
    {
      const auto S = symUOld.getValue(p);

      const Real gamma =
        std::max(gamma_min, std::sqrt(2.0 * dot(S, S)));

      const Real mu =
        m_pl * std::pow(gamma, n_pl - 1.0);

      return std::clamp(mu, mu_min, mu_max);
    };

    Alert::Info() << "Projecting non-Newtonian viscosity ..." << Alert::Raise;

    m_muProjection =
        Integral(m_mu, m_w)
      - Integral(muNonNew, m_w);

    m_muProjection.assemble();
    m_muProjectionSolver.solve();

    const auto convectionTarget = Mult(Jacobian(m_uOld), m_uOld);

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

          const Real tau = vmsScale / invTau;

          // TODO: VERIFY SIGN CONVENTION FOR PROJECTION
          Math::SpatialVector<Real> out(3);
          for (size_t c = 0; c < 3; ++c)
            out[c] = tau * ((conv[c] - proj[c]) + old[c] / m_cfg.dt);

          return out;
        });

    m_subProjection =
        Integral(m_sub, m_vp)
      - Integral(subUpdate, m_vp);

    m_subProjection.assemble();
    m_subProjectionSolver.solve();

    m_flow =
          (m_cfg.rho / m_cfg.dt) * Integral(m_u, m_v)
        - (m_cfg.rho / m_cfg.dt) * Integral(m_uOld, m_v)
        + m_cfg.rho * Integral(Dot(conv_u, m_v))
        + 0.5 * m_cfg.rho * Integral(div_u_old * Dot(m_u, m_v))
        + 2.0 * Integral(m_mu.getSolution() * symU, symV)
        - Integral(m_p, Div(m_v))
        + Integral(Div(m_u), m_q)
        + m_cfg.eps * Integral(m_p, m_q)

        + VMSConvectionBilinearIntegrator(
            m_u, m_v, m_uOld, m_mu.getSolution(),
            m_cfg.rho, m_cfg.dt, c1, c2, vmsScale)

        - VMSConvectionLinearIntegrator(
            m_v, m_sub.getSolution(), m_uOld, m_up.getSolution(), m_mu.getSolution(),
            m_cfg.rho, m_cfg.dt, c1, c2, vmsScale)

        + BoundaryIntegral(pin * Dot(m_v, n)).over(m_cfg.inlet)

        + BoundaryIntegral(m_wk.at(4).pout * Dot(m_v, n)).over(4)
        + BoundaryIntegral(m_wk.at(5).pout * Dot(m_v, n)).over(5)
        + BoundaryIntegral(m_wk.at(6).pout * Dot(m_v, n)).over(6)
        + BoundaryIntegral(m_wk.at(7).pout * Dot(m_v, n)).over(7)
        + BoundaryIntegral(m_wk.at(8).pout * Dot(m_v, n)).over(8)
        + BoundaryIntegral(m_wk.at(9).pout * Dot(m_v, n)).over(9)

        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(m_u, m_v)).over(4)
        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(m_u, m_v)).over(5)
        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(m_u, m_v)).over(6)
        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(m_u, m_v)).over(7)
        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(m_u, m_v)).over(8)
        - BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(m_u, m_v)).over(9)

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

    m_flux = BoundaryIntegral(Dot(m_u.getSolution(), n), m_qFlux).over(m_cfg.inlet);
    m_flux.assemble();
    m_stepData.qIn = m_flux(m_one);

    m_stepData.qOut.clear();
    m_stepData.qOutSum = 0.0;

    for (const Attribute tag : m_cfg.outlets)
    {
      m_flux = BoundaryIntegral(Dot(m_u.getSolution(), n), m_qFlux).over(tag);
      m_flux.assemble();

      const Real qOut = m_flux(m_one);
      m_stepData.qOut[tag] = qOut;
      m_stepData.qOutSum += qOut;

      updateRCR(m_wk[tag], qOut, m_cfg.dt);
      // updateRCRNonNew(m_model, m_wk[tag], qOut, m_cfg.dt);
    }

    m_stepData.flowBalance = m_stepData.qIn + m_stepData.qOutSum;
    m_stepData.t = s.t;
  }

  void CoupledLV0DCoronary3D::updateHistory()
  {
    m_uOld.setData(m_u.getSolution().getData());
    m_pOld.setData(m_p.getSolution().getData());
    m_subOld.setData(m_sub.getSolution().getData());
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
