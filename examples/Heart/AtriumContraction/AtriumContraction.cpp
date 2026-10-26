/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @brief Transient atrial solid activation with prescribed material inputs.
 *
 * The active constitutive law, inertial weak form, and Newmark updates define
 * the example. Geometry, prescribed fields, and output are supporting modules.
 * Coordinates are in mm, time in s, and stresses in kPa.
 */
#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <Rodin/Alert.h>
#include <Rodin/Solid.h>
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Solver/SparseLU.h>
#include "ProblemData.h"
#include "Output.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Examples::Heart::AtriumContraction;

int main(int argc, char** argv)
{
  try
  {
    const auto cfg = readOptions(argc, argv);
    auto mesh = loadWall(cfg);
    ProblemData data(mesh, cfg);
    P1 vectorSpace(mesh, 3);
    const VectorFunction zero(3, [](const Point&) {
      Math::SpatialVector<Real> result(3);
      result.setZero();
      return result;
    });
    GridFunction displacement(vectorSpace);
    GridFunction velocity(vectorSpace), acceleration(vectorSpace);
    GridFunction predictedDisplacement(vectorSpace), predictedVelocity(vectorSpace);
    displacement = zero;
    displacement.setName("Displacement");
    velocity = zero;
    acceleration = zero;
    predictedDisplacement = zero;
    predictedVelocity = zero;
    velocity.setName("Velocity");
    acceleration.setName("Acceleration");
    const Real dt = data.endTime / cfg.steps;
    // Newmark average acceleration, as in the coronary solid dynamics example.
    constexpr Real NewmarkBeta = 0.25;
    constexpr Real NewmarkGamma = 0.5;
    // 1 kPa s^2/mm^2 = 10^9 kg/m^3.
    constexpr Real DensityConversion = 1e-9;
    const Real density = cfg.density * DensityConversion;
    const Real accelerationCoefficient = 1 / (NewmarkBeta * dt * dt);
    const Real massCoefficient = density * accelerationCoefficient;
    Real time = 0;
    const Real mu = cfg.young / (2 * (1 + cfg.poisson));
    const Real lambda =
      cfg.young * cfg.poisson / ((1 + cfg.poisson) * (1 - 2 * cfg.poisson));
    Solid::ActiveFiberLaw::Parameters parameters;
    parameters.stiffness = cfg.activeStiffness;
    Solid::ActiveContraction<Solid::NeoHookean, Solid::ActiveFiberLaw> law(
      Solid::NeoHookean(lambda, mu), Solid::ActiveFiberLaw(parameters));
    const auto input = [&](Solid::ConstitutivePoint& cp) {
      const Point& point = cp.getPoint()->get();
      cp.set<Solid::Tags::FiberDirection>(data.fiberAt(point));
      cp.set<Solid::Tags::ActiveExtension>(
        -cfg.contraction * data.activationAt(point, time));
    };
    TrialFunction increment(vectorSpace);
    TestFunction test(vectorSpace);
    auto work = Solid::InternalVirtualWork(law, displacement).setInput(input);
    // A quadrature policy for this non-polynomial hyperelastic example.
    constexpr size_t QuadratureOrder = 4;
    work.setQuadratureOrder(QuadratureOrder);
    Problem mechanics(increment, test);
    mechanics = work(increment, test) + massCoefficient * Integral(increment, test) +
      massCoefficient * Integral(displacement, test) -
      massCoefficient * Integral(predictedDisplacement, test) +
      DirichletBC(increment, zero)
        .on(Attribute{1}, Attribute{2}, Attribute{3}, Attribute{5}, Attribute{6});
    Solver::SparseLU linearSolver(mechanics);
    Solver::NewtonSolver newton(linearSolver);
    newton.setMaxIterations(30).setAbsoluteTolerance(1e-8).setRelativeTolerance(1e-7);
    newton.setMonitor([](const auto& r) {
      std::cout << "  Newton " << r.iterations << " residual=" << r.finalResidual
                << std::endl;
    });
    Output output(cfg, data, displacement, velocity, acceleration, dt, density,
      NewmarkBeta, NewmarkGamma);
    if (cfg.fibersOnly)
    {
      output.write(0, 0, 0);
      output.close();
      return 0;
    }
    // Initial rest with a consistent acceleration: M a_0 = -f_int(u_0, t_0).
    TrialFunction initialAcceleration(vectorSpace);
    Problem initial(initialAcceleration, test);
    initial = density * Integral(initialAcceleration, test) + work.Residual(test) +
      DirichletBC(initialAcceleration, zero)
        .on(Attribute{1}, Attribute{2}, Attribute{3}, Attribute{5}, Attribute{6});
    Solver::SparseLU initialSolver(initial);
    initialSolver.solve();
    acceleration = initialAcceleration.getSolution();
    const auto& initialSystem = initial.getLinearSystem();
    const Real initialResidual =
      (initialSystem.getOperator() * initialSystem.getSolution() -
        initialSystem.getVector())
        .norm();
    if (!std::isfinite(initialResidual) ||
      initialResidual > 1e-8 * std::max<Real>(1, initialSystem.getVector().norm()))
    {
      Alert::Exception() << "Initial acceleration solve failed" << Alert::Raise;
    }
    const size_t count = cfg.stopAfter ? std::min(cfg.steps, cfg.stopAfter) : cfg.steps;
    for (size_t step = 0; step <= count; ++step)
    {
      time = step * dt;
      std::cout << "Step " << step << '/' << cfg.steps << " time=" << time << std::endl;
      if (step > 0)
      {
        predictedDisplacement =
          displacement + dt * velocity + (dt * dt * (0.5 - NewmarkBeta)) * acceleration;
        predictedVelocity = velocity + (dt * (1 - NewmarkGamma)) * acceleration;
        // Keep the last valid displacement as Newton's initial iterate. The
        // kinematic predictor may fold cells under a coarsely sampled pulse.
        newton.solve(displacement);
        if (!newton.getReport().converged)
        {
          Alert::Exception() << "Transient contraction Newton solve did not converge"
                             << Alert::Raise;
        }
        acceleration = accelerationCoefficient * (displacement - predictedDisplacement);
        velocity = predictedVelocity + (NewmarkGamma * dt) * acceleration;
      }
      output.write(time, step == 0 ? 0 : newton.getReport().iterations,
        step == 0 ? initialResidual : newton.getReport().finalResidual);
    }
    output.close();
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << std::endl;
    return 1;
  }
  return 0;
}
