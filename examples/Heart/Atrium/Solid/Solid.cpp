/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Solid.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <Rodin/Alert.h>
#include <Rodin/Solid.h>
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Solver/SparseLU.h>
#include "Output.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Examples::Heart::AtriumContraction
{
  Solid::Material Solid::makeMaterial(const Parameters& cfg)
  {
    const Real mu = cfg.young / (2 * (1 + cfg.poisson));
    const Real lambda =
      cfg.young * cfg.poisson / ((1 + cfg.poisson) * (1 - 2 * cfg.poisson));
    Rodin::Solid::ActiveFiberLaw::Parameters parameters;
    parameters.stiffness = cfg.activeStiffness;
    return Material(
      Rodin::Solid::NeoHookean(lambda, mu), Rodin::Solid::ActiveFiberLaw(parameters));
  }

  Solid::Solid(const Parameters& parameters)
    : m_parameters(validateOptions(parameters)),
      m_mesh(loadWall(m_parameters)),
      m_data(m_mesh, m_parameters),
      m_space(m_mesh, 3),
      m_displacement(m_space),
      m_velocity(m_space),
      m_acceleration(m_space),
      m_predictedDisplacement(m_space),
      m_predictedVelocity(m_space),
      m_timeStep(m_parameters.timeStep > 0 ? m_parameters.timeStep
                                           : m_data.endTime / m_parameters.steps),
      m_density(m_parameters.density * DensityConversion),
      m_accelerationCoefficient(1 / (NewmarkBeta * m_timeStep * m_timeStep)),
      m_massCoefficient(m_density * m_accelerationCoefficient),
      m_law(makeMaterial(m_parameters)),
      m_increment(m_space),
      m_test(m_space),
      m_work(m_law, m_displacement),
      m_mechanics(m_increment, m_test),
      m_linearSolver(m_mechanics),
      m_newton(m_linearSolver)
  {
    const VectorFunction zero(3, [](const Point&) {
      Math::SpatialVector<Real> result(3);
      result.setZero();
      return result;
    });
    m_displacement = zero;
    m_velocity = zero;
    m_acceleration = zero;
    m_predictedDisplacement = zero;
    m_predictedVelocity = zero;
    m_displacement.setName("Displacement");
    m_velocity.setName("Velocity");
    m_acceleration.setName("Acceleration");

    m_work.setInput([this](Rodin::Solid::ConstitutivePoint& cp) {
      const Point& point = cp.getPoint()->get();
      cp.set<Rodin::Solid::Tags::FiberDirection>(m_data.fiberAt(point));
      cp.set<Rodin::Solid::Tags::ActiveExtension>(
        -m_parameters.contraction * m_data.activationAt(point, m_time));
    });
    m_work.setQuadratureOrder(QuadratureOrder);
    m_mechanics = m_work(m_increment, m_test) +
      m_massCoefficient * Integral(m_increment, m_test) +
      m_massCoefficient * Integral(m_displacement, m_test) -
      m_massCoefficient * Integral(m_predictedDisplacement, m_test) +
      DirichletBC(m_increment, zero)
        .on(Attribute{1}, Attribute{2}, Attribute{3}, Attribute{5}, Attribute{6});
    const Real velocityCoefficient = NewmarkGamma / (NewmarkBeta * m_timeStep);
    if (m_parameters.rayleighAlpha > 0)
    {
      const Real dampingMass = m_parameters.rayleighAlpha * m_density;
      m_mechanics += dampingMass * velocityCoefficient * Integral(m_increment, m_test);
      m_mechanics += dampingMass * Integral(m_predictedVelocity, m_test);
      m_mechanics += dampingMass * velocityCoefficient * Integral(m_displacement, m_test);
      m_mechanics -=
        dampingMass * velocityCoefficient * Integral(m_predictedDisplacement, m_test);
    }
    if (m_parameters.rayleighBeta > 0)
    {
      const Real mu = m_parameters.young / (2 * (1 + m_parameters.poisson));
      const Real lambda = m_parameters.young * m_parameters.poisson /
        ((1 + m_parameters.poisson) * (1 - 2 * m_parameters.poisson));
      const auto symmetricTest = 0.5 * (Jacobian(m_test) + Transpose(Jacobian(m_test)));
      const auto strainRate =
        0.5 * (Jacobian(m_predictedVelocity) + Transpose(Jacobian(m_predictedVelocity))) +
        (0.5 * velocityCoefficient) *
          (Jacobian(m_displacement) + Transpose(Jacobian(m_displacement)) -
            Jacobian(m_predictedDisplacement) -
            Transpose(Jacobian(m_predictedDisplacement)));
      const auto divergenceRate = Div(m_predictedVelocity) +
        velocityCoefficient * (Div(m_displacement) - Div(m_predictedDisplacement));
      const auto symmetricIncrement =
        0.5 * (Jacobian(m_increment) + Transpose(Jacobian(m_increment)));
      auto bulkTangent = Integral(Div(m_increment), Div(m_test));
      auto shearTangent = Integral(symmetricIncrement, symmetricTest);
      auto bulkResidual = Integral(divergenceRate, Div(m_test));
      auto shearResidual = Integral(strainRate, symmetricTest);
      bulkTangent.setOrder(QuadratureOrder);
      shearTangent.setOrder(QuadratureOrder);
      bulkResidual.setOrder(QuadratureOrder);
      shearResidual.setOrder(QuadratureOrder);
      m_mechanics +=
        (m_parameters.rayleighBeta * velocityCoefficient * lambda) * bulkTangent;
      m_mechanics +=
        (2 * m_parameters.rayleighBeta * velocityCoefficient * mu) * shearTangent;
      m_mechanics += (m_parameters.rayleighBeta * lambda) * bulkResidual;
      m_mechanics += (2 * m_parameters.rayleighBeta * mu) * shearResidual;
    }
    m_newton.setMaxIterations(30).setAbsoluteTolerance(1e-8).setRelativeTolerance(1e-7);
    m_newton.setMonitor([](const auto& report) {
      std::cout << "  Newton " << report.iterations
                << " residual=" << report.finalResidual << std::endl;
    });
  }

  Solid::~Solid() = default;

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
    m_output = std::make_unique<Output>(m_parameters, m_data, m_displacement, m_velocity,
      m_acceleration, m_timeStep, m_density, NewmarkBeta, NewmarkGamma);
    Real initialResidual = 0;
    if (!m_parameters.fibersOnly)
    {
      // M a_0 = -f_int(u_0, 0) with zero displacement and velocity.
      const VectorFunction zero(3, [](const Point&) {
        Math::SpatialVector<Real> result(3);
        result.setZero();
        return result;
      });
      Trial initialAcceleration(m_space);
      MechanicsProblem initial(initialAcceleration, m_test);
      initial = m_density * Integral(initialAcceleration, m_test) +
        m_work.Residual(m_test) +
        DirichletBC(initialAcceleration, zero)
          .on(Attribute{1}, Attribute{2}, Attribute{3}, Attribute{5}, Attribute{6});
      LinearSolver initialSolver(initial);
      initialSolver.solve();
      m_acceleration = initialAcceleration.getSolution();
      const auto& system = initial.getLinearSystem();
      initialResidual =
        (system.getOperator() * system.getSolution() - system.getVector()).norm();
      if (!std::isfinite(initialResidual) ||
        initialResidual > 1e-8 * std::max<Real>(1, system.getVector().norm()))
      {
        Alert::Exception() << "Initial acceleration solve failed" << Alert::Raise;
      }
    }
    m_output->write(0, 0, initialResidual);
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
    const size_t nextStep = m_step + 1;
    m_time = nextStep * m_timeStep;
    std::cout << "Step " << nextStep << '/' << m_parameters.steps << " time=" << m_time
              << std::endl;
    m_predictedDisplacement = m_displacement + m_timeStep * m_velocity +
      (m_timeStep * m_timeStep * (0.5 - NewmarkBeta)) * m_acceleration;
    m_predictedVelocity = m_velocity + (m_timeStep * (1 - NewmarkGamma)) * m_acceleration;
    m_newton.solve(m_displacement);
    if (!m_newton.getReport().converged)
    {
      Alert::Exception() << "Transient contraction Newton solve did not converge"
                         << Alert::Raise;
    }
    m_acceleration =
      m_accelerationCoefficient * (m_displacement - m_predictedDisplacement);
    m_velocity = m_predictedVelocity + (NewmarkGamma * m_timeStep) * m_acceleration;
    m_output->write(
      m_time, m_newton.getReport().iterations, m_newton.getReport().finalResidual);
    m_step = nextStep;
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
}
