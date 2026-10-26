/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Solid.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <Rodin/Alert.h>
#include "Output.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Examples::Heart::AtriumContraction
{
  Solid::Material Solid::makeMaterial(const Parameters& cfg)
  {
    Rodin::Solid::ActiveFiberLaw::Parameters active;
    active.stiffness = cfg.activeStiffness;
    active.damping = cfg.activeDamping;
    active.destructionRate = cfg.destructionRate;
    active.crossBridgeStiffness = cfg.crossBridgeStiffness;
    active.contractility = cfg.contractility;
    return Material(
      PassiveMaterial(Rodin::Solid::DemirayGasser(cfg.passive), cfg.viscosity),
      Rodin::Solid::ActiveFiberLaw(active));
  }

  Solid::Solid(const Parameters& parameters)
    : m_parameters(validateOptions(parameters)),
      m_mesh(loadWall(m_parameters)),
      m_data(m_mesh, m_parameters),
      m_space(m_mesh, 3),
      m_activeSpace(m_mesh),
      m_displacement(m_space),
      m_velocity(m_space),
      m_acceleration(m_space),
      m_predictedDisplacement(m_space),
      m_predictedVelocity(m_space),
      m_previousDisplacement(m_space),
      m_extension(m_activeSpace),
      m_gamma(m_activeSpace),
      m_beta(m_activeSpace),
      m_previousExtension(m_activeSpace),
      m_previousGamma(m_activeSpace),
      m_previousBeta(m_activeSpace),
      m_timeStep(m_parameters.timeStep > 0 ? m_parameters.timeStep
                                           : m_data.endTime / m_parameters.steps),
      m_density(m_parameters.density * DensityConversion),
      m_accelerationCoefficient(1 / (NewmarkBeta * m_timeStep * m_timeStep)),
      m_massCoefficient(m_density * m_accelerationCoefficient),
      m_velocityCoefficient(NewmarkGamma / (NewmarkBeta * m_timeStep)),
      m_law(makeMaterial(m_parameters)),
      m_increment(m_space),
      m_extensionIncrement(m_activeSpace),
      m_gammaIncrement(m_activeSpace),
      m_betaIncrement(m_activeSpace),
      m_test(m_space),
      m_extensionTest(m_activeSpace),
      m_gammaTest(m_activeSpace),
      m_betaTest(m_activeSpace),
      m_work(m_law, m_displacement),
      m_mechanics(m_increment, m_extensionIncrement, m_gammaIncrement, m_betaIncrement,
        m_test, m_extensionTest, m_gammaTest, m_betaTest),
      m_linearSolver(m_mechanics),
      m_newton(m_linearSolver)
  {
    const VectorFunction zero(3, [](const Point&) {
      Math::SpatialVector<Real> value(3);
      value.setZero();
      return value;
    });
    m_displacement = zero;
    m_velocity = zero;
    m_acceleration = zero;
    m_predictedDisplacement = zero;
    m_predictedVelocity = zero;
    m_previousDisplacement = zero;
    m_extension = Real(0);
    m_beta = Real(0);
    m_gamma = m_law.getActiveLaw().getMinimumGamma();
    m_previousExtension = m_extension;
    m_previousGamma = m_gamma;
    m_previousBeta = m_beta;
    m_displacement.setName("Displacement");
    m_velocity.setName("Velocity");
    m_acceleration.setName("Acceleration");
    m_extension.setName("ActiveExtension");
    m_gamma.setName("ActiveGamma");
    m_beta.setName("ActiveBeta");

    m_work.setInput([this](Rodin::Solid::ConstitutivePoint& point) { setInputs(point); });
    m_work.setQuadratureOrder(QuadratureOrder);
    m_mechanics = m_work(m_increment, m_test) +
      m_massCoefficient * Integral(m_increment, m_test) +
      m_massCoefficient * Integral(m_displacement, m_test) -
      m_massCoefficient * Integral(m_predictedDisplacement, m_test) +
      DirichletBC(m_increment, zero)
        .on(Attribute{1}, Attribute{2}, Attribute{3}, Attribute{5}, Attribute{6});

    const auto coupling = [this](bool momentum) {
      return MatrixFunction(3, 3, [this, momentum](const Point& point) {
        const auto prepared = prepare(point);
        const Real factor = momentum ? Real(0.5) * prepared.getActiveResponse().dStressDc
                                     : prepared.getActiveResidual().strainDerivative[0];
        return Math::SpatialMatrix<Real>(factor *
          prepared.getKinematicState().getDeformationGradient() *
          prepared.getFiberKinematics().tensor());
      });
    };
    auto momentumCoupling =
      Integral(m_extensionIncrement * coupling(true), Jacobian(m_test));
    auto extensionCoupling =
      Integral(Dot(coupling(false), Jacobian(m_increment)), m_extensionTest);
    momentumCoupling.setOrder(QuadratureOrder);
    extensionCoupling.setOrder(QuadratureOrder);
    m_mechanics += momentumCoupling;
    m_mechanics += extensionCoupling;
    const std::array increments{std::ref(m_extensionIncrement),
      std::ref(m_gammaIncrement), std::ref(m_betaIncrement)};
    const std::array tests{
      std::ref(m_extensionTest), std::ref(m_gammaTest), std::ref(m_betaTest)};
    for (size_t row = 0; row < 3; ++row)
    {
      auto residual = RealFunction([this, row](const Point& point) {
        return prepare(point).getActiveResidual().value[row];
      });
      auto load = Integral(residual, tests[row].get());
      load.setOrder(QuadratureOrder);
      m_mechanics += load;
      for (size_t column = 0; column < 3; ++column)
      {
        auto derivative = RealFunction([this, row, column](const Point& point) {
          return prepare(point).getActiveResidual().derivative(row, column);
        });
        auto block = Integral(derivative * increments[column].get(), tests[row].get());
        block.setOrder(QuadratureOrder);
        m_mechanics += block;
      }
    }
    m_newton.setMaxIterations(MaximumNewtonIterations)
      .setAbsoluteTolerance(AbsoluteResidualTolerance)
      .setRelativeTolerance(0);
    m_newton.setStepPolicy([this](auto& values, const auto& system, const auto&) {
      const Math::Vector<Real> saved = values, increment = system.getSolution();
      const Real residual = system.getVector().norm();
      Real scale = 1;
      for (size_t attempt = 0; attempt < MaximumBacktracks; ++attempt)
      {
        values = saved + scale * increment;
        scatter();
        if (isAdmissible())
        {
          m_mechanics.assemble();
          const Real trialResidual = m_mechanics.getLinearSystem().getVector().norm();
          if (std::isfinite(trialResidual) && trialResidual < residual)
          {
            return Newton::StepResult{
              .accepted = true, .stepNorm = scale * increment.norm()};
          }
        }
        scale *= Real(0.5);
      }
      values = saved;
      scatter();
      return Newton::StepResult{.accepted = false, .stepNorm = 0};
    });
    m_newton.setMonitor([](const auto& report) {
      std::cout << "  Newton " << report.iterations
                << " residual=" << report.finalResidual
                << " damping=" << report.dampingFactor << std::endl;
    });
  }

  Solid::~Solid() = default;

  void Solid::setInputs(Rodin::Solid::ConstitutivePoint& point) const
  {
    using namespace Rodin::Solid;
    const Point& reference = point.getPoint()->get();
    const auto fiber = m_data.fiberAt(reference);
    point.set<Tags::FiberDirection>(fiber);
    point.set<Tags::VelocityGradient>(
      Math::SpatialMatrix<Real>(Jacobian(m_predictedVelocity).getValue(reference) +
        m_velocityCoefficient *
          (point.getKinematicState().getDisplacementGradient() -
            Jacobian(m_predictedDisplacement).getValue(reference))));
    point.set<Tags::VelocityGradientDerivative>(m_velocityCoefficient);
    point.set<Tags::TimeStep>(m_timeStep);
    point.set<Tags::ActiveExtension>(m_extension(reference));
    point.set<Tags::ActiveGamma>(m_gamma(reference));
    point.set<Tags::ActiveBeta>(m_beta(reference));
    point.set<Tags::PreviousActiveExtension>(m_previousExtension(reference));
    point.set<Tags::PreviousActiveGamma>(m_previousGamma(reference));
    point.set<Tags::PreviousActiveBeta>(m_previousBeta(reference));
    KinematicState previous(3);
    previous.setDisplacementGradient(
      Jacobian(m_previousDisplacement).getValue(reference));
    point.set<Tags::PreviousFiberStrain>(FiberKinematics(previous, fiber).strain());
    point.set<Tags::ElectricalActivation>(m_data.electricalAt(reference, m_time));
  }

  Solid::Material Solid::prepare(const Point& point) const
  {
    Rodin::Solid::KinematicState state(3);
    state.setDisplacementGradient(Jacobian(m_displacement).getValue(point));
    Rodin::Solid::ConstitutivePoint cp(point, state);
    setInputs(cp);
    auto law = m_law;
    law.setConstitutivePoint(cp);
    return law;
  }

  void Solid::pack()
  {
    m_iterate.resize(m_displacement.getData().size() + 3 * m_extension.getData().size());
    m_iterate << m_displacement.getData(), m_extension.getData(), m_gamma.getData(),
      m_beta.getData();
  }

  void Solid::scatter()
  {
    const auto mechanics = m_displacement.getData().size(),
               active = m_extension.getData().size();
    m_displacement.getData() = m_iterate.head(mechanics);
    m_extension.getData() = m_iterate.segment(mechanics, active);
    m_gamma.getData() = m_iterate.segment(mechanics + active, active);
    m_beta.getData() = m_iterate.tail(active);
  }

  bool Solid::isAdmissible() const
  {
    if (!m_iterate.allFinite() || m_gamma.getData().minCoeff() <= 0 ||
      m_extension.getData().minCoeff() <= -Real(0.5))
      return false;
    const auto minimum = m_data.minimumDeformedJacobian(m_displacement);
    return minimum && *minimum > 0;
  }

  Solid& Solid::initialize()
  {
    if (m_finished || m_failed)
    {
      Alert::Exception() << "Cannot initialize a finished or failed atrial solid"
                         << Alert::Raise;
    }
    if (m_initialized)
      return *this;
    m_failed = true;
    Real residual = 0;
    if (!m_parameters.fibersOnly)
    {
      // With no fluid pressure, passive equilibrium has traction-free wall faces.
      auto passive = m_law.getPassiveLaw();
      auto work = Rodin::Solid::InternalVirtualWork(passive, m_displacement);
      work.setInput([this](Rodin::Solid::ConstitutivePoint& cp) {
        const auto& point = cp.getPoint()->get();
        Math::SpatialMatrix<Real> zero(3, 3);
        zero.setZero();
        cp.set<Rodin::Solid::Tags::FiberDirection>(m_data.fiberAt(point));
        cp.set<Rodin::Solid::Tags::VelocityGradient>(zero);
        cp.set<Rodin::Solid::Tags::VelocityGradientDerivative>(0);
      });
      work.setQuadratureOrder(QuadratureOrder);
      const VectorFunction zero(3, [](const Point&) {
        Math::SpatialVector<Real> value(3);
        value.setZero();
        return value;
      });
      Problem equilibrium(m_increment, m_test);
      equilibrium = work(m_increment, m_test) +
        DirichletBC(m_increment, zero)
          .on(Attribute{1}, Attribute{2}, Attribute{3}, Attribute{5}, Attribute{6});
      LinearSolver linear(equilibrium);
      Newton newton(linear);
      newton.setMaxIterations(MaximumNewtonIterations)
        .setAbsoluteTolerance(AbsoluteResidualTolerance)
        .setRelativeTolerance(0);
      newton.solve(m_displacement);
      if (!newton.getReport().converged)
      {
        Alert::Exception() << "Passive atrial equilibrium did not converge"
                           << Alert::Raise;
      }
      residual = newton.getReport().finalResidual;
    }
    m_extension = [this](const Point& point) {
      Rodin::Solid::KinematicState state(3);
      state.setDisplacementGradient(Jacobian(m_displacement).getValue(point));
      return Rodin::Solid::FiberKinematics(state, m_data.fiberAt(point)).strain();
    };
    m_previousDisplacement = m_displacement;
    m_previousExtension = m_extension;
    m_previousGamma = m_gamma;
    m_previousBeta = m_beta;
    m_predictedDisplacement = m_displacement;
    m_output = std::make_unique<Output>(m_parameters, m_data, m_displacement, m_velocity,
      m_acceleration, m_extension, m_gamma, m_beta, m_timeStep, m_density, NewmarkBeta,
      NewmarkGamma);
    m_output->write(0, 0, residual);
    m_initialized = true;
    m_failed = false;
    return *this;
  }

  void Solid::step()
  {
    if (m_finished || m_failed)
    {
      Alert::Exception() << "Cannot advance a finished or failed atrial solid"
                         << Alert::Raise;
    }
    if (!m_initialized)
      initialize();
    if (m_parameters.fibersOnly)
    {
      Alert::Exception() << "Fiber inspection mode does not advance mechanics"
                         << Alert::Raise;
    }
    m_failed = true;
    m_previousDisplacement = m_displacement;
    m_previousExtension = m_extension;
    m_previousGamma = m_gamma;
    m_previousBeta = m_beta;
    m_time = (m_step + 1) * m_timeStep;
    std::cout << "Step " << (m_step + 1) << '/' << m_parameters.steps
              << " time=" << m_time << std::endl;
    m_predictedDisplacement = m_displacement + m_timeStep * m_velocity +
      (m_timeStep * m_timeStep * (Real(0.5) - NewmarkBeta)) * m_acceleration;
    m_predictedVelocity = m_velocity + (m_timeStep * (1 - NewmarkGamma)) * m_acceleration;
    for (auto cell = m_mesh.getCell(); cell; ++cell)
    {
      const Point point(*cell, Polytope::Traits(cell->getGeometry()).getCentroid());
      Rodin::Solid::ActiveFiberLaw::State previous;
      previous.gamma = m_previousGamma(point);
      previous.beta = m_previousBeta(point);
      const auto predicted =
        m_law.getActiveLaw().update(m_timeStep, previous, m_previousExtension(point),
          m_extension(point), m_data.maximumElectrical(*cell, m_time));
      const Index index =
        m_activeSpace.getGlobalIndex({m_mesh.getDimension(), cell->getIndex()}, 0);
      m_gamma[index] = predicted.gamma;
      m_beta[index] = predicted.beta;
    }
    pack();
    m_newton.solve(m_iterate);
    if (!m_newton.getReport().converged)
    {
      Alert::Exception() << "Coupled atrial solid Newton solve did not converge"
                         << Alert::Raise;
    }
    m_acceleration =
      m_accelerationCoefficient * (m_displacement - m_predictedDisplacement);
    m_velocity = m_predictedVelocity + (NewmarkGamma * m_timeStep) * m_acceleration;
    m_output->write(
      m_time, m_newton.getReport().iterations, m_newton.getReport().finalResidual);
    ++m_step;
    m_failed = false;
  }

  void Solid::finish()
  {
    if (m_finished)
      return;
    if (m_output)
      m_output->close();
    m_finished = true;
  }

  void Solid::run()
  {
    initialize();
    const size_t count = m_parameters.stopAfter
      ? std::min(m_parameters.steps, m_parameters.stopAfter)
      : m_parameters.steps;
    if (!m_parameters.fibersOnly)
    {
      while (m_step < count)
        step();
    }
    finish();
  }

  Real Solid::getTime() const noexcept
  {
    return m_time;
  }
  Real Solid::getTimeStep() const noexcept
  {
    return m_timeStep;
  }
  size_t Solid::getStepCount() const noexcept
  {
    return m_step;
  }
  const Solid::Parameters& Solid::getParameters() const noexcept
  {
    return m_parameters;
  }
  const LocalMesh& Solid::getMesh() const noexcept
  {
    return m_mesh;
  }
  const VectorField& Solid::getDisplacement() const noexcept
  {
    return m_displacement;
  }
  const VectorField& Solid::getVelocity() const noexcept
  {
    return m_velocity;
  }
  const VectorField& Solid::getAcceleration() const noexcept
  {
    return m_acceleration;
  }
  const CellScalarField& Solid::getActiveExtension() const noexcept
  {
    return m_extension;
  }
  const CellScalarField& Solid::getActiveGamma() const noexcept
  {
    return m_gamma;
  }
  const CellScalarField& Solid::getActiveBeta() const noexcept
  {
    return m_beta;
  }
}
