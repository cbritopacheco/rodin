#include "CoupledLV0DCoronary3D.h"

#include <cmath>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <type_traits>

#include <Rodin/Alert.h>
#include <Rodin/Solver.h>

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
      m_model(std::in_place, m_input)
  {
  }

  CoupledLV0DCoronary3D::~CoupledLV0DCoronary3D() = default;

  void CoupledLV0DCoronary3D::updateRCR(RCR& bc, Real Q, Real dt)
  {
    const Real a = bc.C / dt;
    bc.pc = (a * bc.pc + Q + bc.pd / bc.Rd) / (a + 1.0 / bc.Rd);
    bc.pout = bc.pc + bc.Rp * Q;
  }

  CoupledLV0DCoronary3D::Real CoupledLV0DCoronary3D::periodic_activation(Real t)
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

  CoupledLV0DCoronary3D::Real CoupledLV0DCoronary3D::atrial_pressure(Real t)
  {
    const Real T = 0.85;
    const Real tau = t - T * std::floor(t / T);

    const Real min_value = 500.0;
    const Real max_value = 1000.0;
    const Real second_threshold = 1250.0;

    const Real t1 = 0.02;
    const Real t2 = 0.15;
    const Real t3 = 0.17;
    const Real t4 = 0.56;
    const Real t5 = 0.62;
    const Real t6 = 0.85;

    Real alpha = 0.0;
    Real value = min_value;

    if (tau < t1)
    {
      alpha = -(tau - t1) / t1;
      value = alpha * min_value + (1.0 - alpha) * max_value;
    }
    else if (tau < t2)
    {
      value = max_value;
    }
    else if (tau < t3)
    {
      alpha = -(tau - t3) / (t3 - t2);
      value = alpha * max_value + (1.0 - alpha) * min_value;
    }
    else if (tau < t4)
    {
      alpha = -(tau - t4) / (t4 - t3);
      value = alpha * min_value + (1.0 - alpha) * second_threshold;
    }
    else if (tau < t5)
    {
      value = second_threshold;
    }
    else if (tau < t6)
    {
      alpha = -(tau - t6) / (t6 - t5);
      value = alpha * second_threshold + (1.0 - alpha) * min_value;
    }
    else
    {
      value = min_value;
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
    m_input.rho = 1.0e3;
    m_input.R0 = 2.36e-2;
    m_input.d0 = 1.42e-2;

    m_input.Es = 3.0e5;
    m_input.mu = 70.0;
    m_input.eta = 70.0;
    m_input.alpha = 3.0;
    m_input.k0 = 1.0e5;
    m_input.sigma0 = 5.0e5;

    m_input.Rp = 8.0e6;
    m_input.Cp = 5.0e-9;
    m_input.Rd = 1.0e8;
    m_input.Cd = 1.0e-8;

    m_input.Kat = 8.0e-7;
    m_input.Kp  = 5.0e-10;
    m_input.Kar = 1.3e-5;

    m_input.cavityCapacity = 5.0e-12;

    m_input.localTolerance = 1e-12;
    m_input.localMaxIterations = 50;
    m_input.localDamping = 1.0;
    m_input.absRegularization = 1e-14;

    m_input.initFibDef = 0.0;
    m_input.initActiveStiffness = 0.0;
    m_input.initActiveStress = 0.0;

    m_input.pSv = [](Real) { return 1.0e3; };
    m_input.pAt = atrial_pressure;
    m_input.u = periodic_activation;

    {
      using PassiveEnergy = std::decay_t<decltype(m_input.passiveEnergy)>;
      typename PassiveEnergy::Parameters hp;
      hp.mu1 = 0.0;
      hp.mu2 = 0.0;
      hp.C0 = 1.9e3;
      hp.C1 = 1.1e-1;
      hp.C2 = 1.9e3;
      hp.C3 = 1.1e-1;
      m_input.passiveEnergy = PassiveEnergy(hp);
    }

    m_model.emplace(m_input);
    m_model->setMaxIterations(200)
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

    m_model->initialize(s0);
  }

  void CoupledLV0DCoronary3D::setupMeshAndSpaces()
  {
    m_mesh.load(m_cfg.meshPath, IO::FileFormat::MEDIT);
    m_mesh.scale(m_cfg.meshScale);

    Alert::Info() << "Computing connectivity for " << m_cfg.meshPath << " ..."
                  << Alert::Raise;

    m_mesh.getConnectivity().compute(3, 2);
    m_mesh.getConnectivity().compute(2, 1);
    m_mesh.getConnectivity().compute(1, 0);
    m_mesh.getConnectivity().compute(2, 3);

    Alert::Info() << "Setting up " << m_cfg.xdmfBasename << ".xdmf ..."
                  << Alert::Raise;

    m_xdmf = std::make_unique<IO::XDMF>(m_cfg.xdmfBasename);
    m_xdmf->setMesh(m_mesh);

    const size_t dim = m_mesh.getSpaceDimension();
    if (dim != 3)
      throw std::runtime_error("Expected a 3D coronary mesh.");

    m_uh = std::make_unique<VelocityFESType>(std::integral_constant<size_t, 2>{}, m_mesh, dim);
    m_ph = std::make_unique<PressureFESType>(std::integral_constant<size_t, 1>{}, m_mesh);

    m_u = std::make_unique<VelocityTrialFunctionType>(*m_uh);
    m_p = std::make_unique<PressureTrialFunctionType>(*m_ph);
    m_v = std::make_unique<VelocityTestFunctionType>(*m_uh);
    m_q = std::make_unique<PressureTestFunctionType>(*m_ph);
    m_u->setName("u");
    m_p->setName("p");

    m_uOld = std::make_unique<VelocityGridFunctionType>(*m_uh);
    m_pOld = std::make_unique<PressureGridFunctionType>(*m_ph);
    *m_uOld = Math::Vector<Real>{{0.0, 0.0, 0.0}};
    *m_pOld = 0.0;

    m_one = std::make_unique<PressureGridFunctionType>(*m_ph);
    *m_one = 1.0;
    m_qFlux = std::make_unique<PressureTestFunctionType>(*m_ph);
    m_flux = std::make_unique<FluxLinearFormType>(*m_qFlux);

    m_xdmf->add("velocity", m_u->getSolution());
    m_xdmf->add("pressure", m_p->getSolution());

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
    const auto& s = m_model->getState();
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
    const auto rep = m_model->step(m_cfg.dt);
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
    const auto& s = m_model->getState();
    const Real pin = s.par;

    const auto n = BoundaryNormal(m_mesh);
    const auto conv_u = Mult(Jacobian(*m_u), *m_uOld);
    const auto div_u_old = Div(*m_uOld);
    const auto beta = Max(-Dot(*m_uOld, n), 0.0);
    auto symU = 0.5 * (Jacobian(*m_u) + Transpose(Jacobian(*m_u)));
    auto symV = 0.5 * (Jacobian(*m_v) + Transpose(Jacobian(*m_v)));

    Problem flow(*m_u, *m_p, *m_v, *m_q);
    flow =
          (m_cfg.rho / m_cfg.dt) * Integral(*m_u, *m_v)
        - (m_cfg.rho / m_cfg.dt) * Integral(*m_uOld, *m_v)
        + m_cfg.rho * Integral(Dot(conv_u, *m_v))
        + 0.5 * m_cfg.rho * Integral(div_u_old * Dot(*m_u, *m_v))
        + 2 * m_cfg.mu * Integral(symU, symV)
        - Integral(*m_p, Div(*m_v))
        + Integral(Div(*m_u), *m_q)
        + m_cfg.eps * Integral(*m_p, *m_q)
        + BoundaryIntegral(pin * Dot(*m_v, n)).over(m_cfg.inlet)
        + BoundaryIntegral(m_wk.at(4).pout * Dot(*m_v, n)).over(4)
        + BoundaryIntegral(m_wk.at(5).pout * Dot(*m_v, n)).over(5)
        + BoundaryIntegral(m_wk.at(6).pout * Dot(*m_v, n)).over(6)
        + BoundaryIntegral(m_wk.at(7).pout * Dot(*m_v, n)).over(7)
        + BoundaryIntegral(m_wk.at(8).pout * Dot(*m_v, n)).over(8)
        + BoundaryIntegral(m_wk.at(9).pout * Dot(*m_v, n)).over(9)
        + BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(*m_u, *m_v)).over(4)
        + BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(*m_u, *m_v)).over(5)
        + BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(*m_u, *m_v)).over(6)
        + BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(*m_u, *m_v)).over(7)
        + BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(*m_u, *m_v)).over(8)
        + BoundaryIntegral(0.5 * m_cfg.rho * beta * Dot(*m_u, *m_v)).over(9)
        + DirichletBC(*m_u, Zero(m_mesh.getSpaceDimension())).on(m_cfg.wall);

    Alert::Info() << "Assembling 3D time step ..." << Alert::Raise;
    flow.assemble().setFieldSplits();

    Alert::Info() << "Solving 3D time step ..." << Alert::Raise;
    Solver::KSP(flow).solve();
  }

  void CoupledLV0DCoronary3D::computeFluxesAndUpdateRCR()
  {
    const auto n = BoundaryNormal(m_mesh);
    const auto& s = m_model->getState();

    (*m_flux) = BoundaryIntegral(Dot(m_u->getSolution(), n), *m_qFlux).over(m_cfg.inlet);
    m_flux->assemble();
    m_stepData.qIn = (*m_flux)(*m_one);

    m_stepData.qOut.clear();
    m_stepData.qOutSum = 0.0;
    for (const Attribute tag : m_cfg.outlets)
    {
      (*m_flux) = BoundaryIntegral(Dot(m_u->getSolution(), n), *m_qFlux).over(tag);
      m_flux->assemble();
      const Real qOut = (*m_flux)(*m_one);
      m_stepData.qOut[tag] = qOut;
      m_stepData.qOutSum += qOut;
      updateRCR(m_wk[tag], qOut, m_cfg.dt);
    }

    m_stepData.flowBalance = m_stepData.qIn + m_stepData.qOutSum;
    m_stepData.t = s.t;
  }

  void CoupledLV0DCoronary3D::updateHistory()
  {
    m_uOld->setData(m_u->getSolution().getData());
    m_pOld->setData(m_p->getSolution().getData());
  }

  void CoupledLV0DCoronary3D::writeOutputs() const
  {
    m_xdmf->write(m_model->getState().t).flush();
  }

  CoupledLV0DCoronary3D::StepData CoupledLV0DCoronary3D::collectStepData() const
  {
    StepData d;
    const auto& s = m_model->getState();
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

    for (int i = 0; i < m_cfg.nsteps; ++i)
    {
      std::cout << "Step " << i << ": t = " << m_model->getState().t << "\n";
      if (!advance0D())
      {
        std::cerr << "0D solver failed to converge at step "
                  << i << ", t = " << m_model->getState().t << "\n";
        return 1;
      }

      solve3D();
      computeFluxesAndUpdateRCR();
      writeCSVRow();
      updateHistory();
      writeOutputs();
    }

    if (m_xdmf)
      m_xdmf->close();
    return 0;
  }
}
