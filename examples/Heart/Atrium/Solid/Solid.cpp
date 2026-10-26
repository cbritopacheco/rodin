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
  namespace
  {
    constexpr Real NewmarkBeta = 0.25;
    constexpr Real NewmarkGamma = 0.5;
    // 1 kPa s^2/mm^2 = 10^9 kg/m^3.
    constexpr Real DensityConversion = 1e-9;
    constexpr size_t QuadratureOrder = 4;

    using Material = Rodin::Solid::ActiveContraction<Rodin::Solid::NeoHookean,
      Rodin::Solid::ActiveFiberLaw>;
    using Trial = TrialFunction<VectorField, VectorSpace>;
    using Test = TestFunction<VectorSpace>;
    using LinearSystem = Math::LinearSystem<Math::SparseMatrix<Real>, Math::Vector<Real>>;
    using MechanicsProblem = Problem<LinearSystem, Trial, Test>;
    using LinearSolver = Solver::SparseLU<LinearSystem>;
    using Newton = Solver::NewtonSolver<LinearSolver>;
    using InternalWork = Rodin::Solid::InternalVirtualWork<Material, VectorField>;

    Material makeMaterial(const Options& cfg)
    {
      const Real mu = cfg.young / (2 * (1 + cfg.poisson));
      const Real lambda =
        cfg.young * cfg.poisson / ((1 + cfg.poisson) * (1 - 2 * cfg.poisson));
      Rodin::Solid::ActiveFiberLaw::Parameters parameters;
      parameters.stiffness = cfg.activeStiffness;
      return Material(
        Rodin::Solid::NeoHookean(lambda, mu), Rodin::Solid::ActiveFiberLaw(parameters));
    }
  }

  struct Solid::Impl
  {
    Options cfg;
    LocalMesh mesh;
    ProblemData data;
    VectorSpace space;
    VectorField displacement, velocity, acceleration;
    VectorField predictedDisplacement, predictedVelocity;
    Real dt, density, accelerationCoefficient, massCoefficient;
    Real time = 0;
    size_t step = 0;
    Material law;
    Trial increment;
    Test test;
    InternalWork work;
    MechanicsProblem mechanics;
    LinearSolver linearSolver;
    Newton newton;
    std::unique_ptr<Output> output;
    bool initialized = false, finished = false, failed = false;

    explicit Impl(const Options& parameters)
      : cfg(validateOptions(parameters)),
        mesh(loadWall(cfg)),
        data(mesh, cfg),
        space(mesh, 3),
        displacement(space),
        velocity(space),
        acceleration(space),
        predictedDisplacement(space),
        predictedVelocity(space),
        dt(cfg.timeStep > 0 ? cfg.timeStep : data.endTime / cfg.steps),
        density(cfg.density * DensityConversion),
        accelerationCoefficient(1 / (NewmarkBeta * dt * dt)),
        massCoefficient(density * accelerationCoefficient),
        law(makeMaterial(cfg)),
        increment(space),
        test(space),
        work(law, displacement),
        mechanics(increment, test),
        linearSolver(mechanics),
        newton(linearSolver)
    {
      const VectorFunction zero(3, [](const Point&) {
        Math::SpatialVector<Real> result(3);
        result.setZero();
        return result;
      });
      displacement = zero;
      velocity = zero;
      acceleration = zero;
      predictedDisplacement = zero;
      predictedVelocity = zero;
      displacement.setName("Displacement");
      velocity.setName("Velocity");
      acceleration.setName("Acceleration");

      work.setInput([this](Rodin::Solid::ConstitutivePoint& cp) {
        const Point& point = cp.getPoint()->get();
        cp.set<Rodin::Solid::Tags::FiberDirection>(data.fiberAt(point));
        cp.set<Rodin::Solid::Tags::ActiveExtension>(
          -cfg.contraction * data.activationAt(point, time));
      });
      work.setQuadratureOrder(QuadratureOrder);
      mechanics = work(increment, test) + massCoefficient * Integral(increment, test) +
        massCoefficient * Integral(displacement, test) -
        massCoefficient * Integral(predictedDisplacement, test) +
        DirichletBC(increment, zero)
          .on(Attribute{1}, Attribute{2}, Attribute{3}, Attribute{5}, Attribute{6});
      newton.setMaxIterations(30).setAbsoluteTolerance(1e-8).setRelativeTolerance(1e-7);
      newton.setMonitor([](const auto& report) {
        std::cout << "  Newton " << report.iterations
                  << " residual=" << report.finalResidual << std::endl;
      });
    }

    void initialize()
    {
      if (finished || failed)
      {
        Alert::Exception() << "Cannot initialize a finished or failed atrial solid"
                           << Alert::Raise;
      }
      if (initialized)
        return;
      failed = true;
      output = std::make_unique<Output>(cfg, data, displacement, velocity, acceleration,
        dt, density, NewmarkBeta, NewmarkGamma);
      Real initialResidual = 0;
      if (!cfg.fibersOnly)
      {
        // M a_0 = -f_int(u_0, 0) with zero displacement and velocity.
        const VectorFunction zero(3, [](const Point&) {
          Math::SpatialVector<Real> result(3);
          result.setZero();
          return result;
        });
        Trial initialAcceleration(space);
        MechanicsProblem initial(initialAcceleration, test);
        initial = density * Integral(initialAcceleration, test) + work.Residual(test) +
          DirichletBC(initialAcceleration, zero)
            .on(Attribute{1}, Attribute{2}, Attribute{3}, Attribute{5}, Attribute{6});
        LinearSolver initialSolver(initial);
        initialSolver.solve();
        acceleration = initialAcceleration.getSolution();
        const auto& system = initial.getLinearSystem();
        initialResidual =
          (system.getOperator() * system.getSolution() - system.getVector()).norm();
        if (!std::isfinite(initialResidual) ||
          initialResidual > 1e-8 * std::max<Real>(1, system.getVector().norm()))
        {
          Alert::Exception() << "Initial acceleration solve failed" << Alert::Raise;
        }
      }
      output->write(0, 0, initialResidual);
      initialized = true;
      failed = false;
    }

    void advance()
    {
      if (finished || failed)
      {
        Alert::Exception() << "Cannot advance a finished or failed atrial solid"
                           << Alert::Raise;
      }
      if (!initialized)
        initialize();
      if (cfg.fibersOnly)
      {
        Alert::Exception() << "Fiber inspection mode does not advance mechanics"
                           << Alert::Raise;
      }
      failed = true;
      const size_t nextStep = step + 1;
      time = nextStep * dt;
      std::cout << "Step " << nextStep << '/' << cfg.steps << " time=" << time
                << std::endl;
      predictedDisplacement =
        displacement + dt * velocity + (dt * dt * (0.5 - NewmarkBeta)) * acceleration;
      predictedVelocity = velocity + (dt * (1 - NewmarkGamma)) * acceleration;
      newton.solve(displacement);
      if (!newton.getReport().converged)
      {
        Alert::Exception() << "Transient contraction Newton solve did not converge"
                           << Alert::Raise;
      }
      acceleration = accelerationCoefficient * (displacement - predictedDisplacement);
      velocity = predictedVelocity + (NewmarkGamma * dt) * acceleration;
      output->write(
        time, newton.getReport().iterations, newton.getReport().finalResidual);
      step = nextStep;
      failed = false;
    }

    void finish()
    {
      if (finished)
        return;
      if (output)
        output->close();
      finished = true;
    }
  };

  Solid::Solid(const Parameters& parameters)
    : m_impl(std::make_unique<Impl>(parameters))
  {}
  Solid::~Solid() = default;
  Solid& Solid::initialize()
  {
    m_impl->initialize();
    return *this;
  }
  void Solid::step()
  {
    m_impl->advance();
  }
  void Solid::finish()
  {
    m_impl->finish();
  }
  void Solid::run()
  {
    initialize();
    const size_t count = m_impl->cfg.stopAfter
      ? std::min(m_impl->cfg.steps, m_impl->cfg.stopAfter)
      : m_impl->cfg.steps;
    if (!m_impl->cfg.fibersOnly)
    {
      while (m_impl->step < count)
        step();
    }
    finish();
  }
  Real Solid::getTime() const noexcept
  {
    return m_impl->time;
  }
  Real Solid::getTimeStep() const noexcept
  {
    return m_impl->dt;
  }
  size_t Solid::getStepCount() const noexcept
  {
    return m_impl->step;
  }
  const Solid::Parameters& Solid::getParameters() const noexcept
  {
    return m_impl->cfg;
  }
  const LocalMesh& Solid::getMesh() const noexcept
  {
    return m_impl->mesh;
  }
  const VectorField& Solid::getDisplacement() const noexcept
  {
    return m_impl->displacement;
  }
  const VectorField& Solid::getVelocity() const noexcept
  {
    return m_impl->velocity;
  }
  const VectorField& Solid::getAcceleration() const noexcept
  {
    return m_impl->acceleration;
  }
}
