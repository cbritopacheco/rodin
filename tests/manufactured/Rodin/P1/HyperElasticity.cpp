/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Hyperelasticity manufactured solution tests.
 *
 * These tests assemble Rodin variational forms for a hyperelasticity manufactured solution, solve the problem on the configured mesh, and compare against analytic fields or expected residual/error behavior. They protect the P1 finite-element and solver path, including boundary-condition handling, geometry coverage, and numerical accuracy of the manufactured workflow.
 */

#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Solid.h"
#include "Rodin/Solver/NewtonSolver.h"
#include "Rodin/Solver/SparseLU.h"
#include "Rodin/Test/FDProbe.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

namespace Rodin::Tests::Manufactured::HyperElasticity
{
  namespace
  {
    enum class ResidualSign
    {
      Correct,
      Wrong
    };

    struct SolveResult
    {
      Real l2ErrorSquared;
      bool converged;
      Real initialResidual;
      Real finalResidual;
    };

    auto makeUnitSquareMesh()
    {
      constexpr Real meshElementSize = 1.0 / 16.0;
      Mesh mesh;
      mesh = mesh.UniformGrid(Polytope::Type::Quadrilateral, { 17, 17 });
      mesh.scale(meshElementSize);
      mesh.getConnectivity().compute(1, 2);
      return mesh;
    }

    auto makeUnitCubeMesh()
    {
      constexpr Real meshElementSize = 1.0 / 4.0;
      Mesh mesh;
      mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {5, 5, 5});
      mesh.scale(meshElementSize);
      mesh.getConnectivity().compute(2, 3);
      return mesh;
    }

    template <class FES, class State>
    Rodin::Test::FDProbeReport checkNeoHookeanInternalVirtualWorkFD(
      FES& Vh, State& uCurrent)
    {
      const Real lambda = 2.0;
      const Real mu = 1.0;
      Solid::NeoHookean law(lambda, mu);

      TrialFunction du(Vh);
      TestFunction v(Vh);
      du.getSolution() = uCurrent;

      auto ivw = Solid::InternalVirtualWork(law, du.getSolution());
      Problem problem(du, v);
      problem = ivw(du, v);

      Rodin::Test::FDProbe probe(problem);
      return probe.test(1e-6);
    }

    template <class FES, class State>
    Rodin::Test::FDProbeReport checkActiveContractionInternalVirtualWorkFD(
      FES& Vh, State& uCurrent)
    {
      const size_t dim = Vh.getMesh().getSpaceDimension();

      Solid::NeoHookean passive(0.0, 0.0);
      Solid::ActiveFiberLaw::Parameters activeInput;
      activeInput.stiffness = 80.0;
      activeInput.damping = 0.4;
      activeInput.destructionRate = 0.5;
      activeInput.crossBridgeStiffness = 60.0;
      activeInput.contractility = 40.0;
      Solid::ActiveFiberLaw active(activeInput);
      Solid::ActiveContraction law(passive, active);
      law.setLocalTolerance(1e-13).setLocalMaxIterations(80);

      auto setActiveInput = [dim](Solid::ConstitutivePoint& cp) {
        Math::SpatialVector<Real> fiber(static_cast<std::uint8_t>(dim));
        fiber.setZero();
        fiber[0] = 0.8;
        fiber[1] = 0.6;
        if (dim == 3)
          fiber[2] = 0.3;
        fiber *= 1.0 / fiber.norm();

        cp.set<Solid::Tags::FiberDirection>(fiber);
        cp.set<Solid::Tags::TimeStep>(0.01);
        cp.set<Solid::Tags::PreviousActiveExtension>(-0.035);
        cp.set<Solid::Tags::PreviousActiveGamma>(1.7);
        cp.set<Solid::Tags::PreviousActiveBeta>(0.9);
        cp.set<Solid::Tags::ElectricalActivation>(1.1);
      };

      TrialFunction du(Vh);
      TestFunction v(Vh);
      du.getSolution() = uCurrent;

      auto ivw =
        Solid::InternalVirtualWork(law, du.getSolution()).setInput(setActiveInput);
      Problem problem(du, v);
      problem = ivw(du, v);

      Rodin::Test::FDProbe probe(problem);
      return probe.test(1e-7);
    }

    SolveResult solveAffineHyperElasticity(ResidualSign residualSign)
    {
      Mesh mesh = makeUnitSquareMesh();
      const size_t dim = mesh.getSpaceDimension();
      P1 Vh(mesh, dim);

      const Real lambda = 2.0;
      const Real mu = 1.0;
      Solid::NeoHookean law(lambda, mu);

      GridFunction uCurrent(Vh);
      const Real xDisplacementScale = 0.10;
      const Real xDisplacementOffset = 0.05;
      const Real yDisplacementScale = -0.08;
      const Real yDisplacementOffset = 0.10;
      const Real perturbationScale = 0.03;
      const Real pi = Math::Constants::pi();
      // Affine displacement is represented exactly in P1, making this a strict
      // manufactured test for the nonlinear solver/integrator plumbing.
      auto uExact = VectorFunction{
        xDisplacementScale * F::x + xDisplacementOffset,
        yDisplacementScale * F::y + yDisplacementOffset
      };
      auto zero = VectorFunction{ Zero(), Zero() };
      auto perturbation = perturbationScale * sin(pi * F::x) * sin(pi * F::y);
      uCurrent = VectorFunction{
        xDisplacementScale * F::x + xDisplacementOffset + perturbation,
        yDisplacementScale * F::y + yDisplacementOffset - perturbation
      };

      TrialFunction du(Vh);
      TestFunction v(Vh);

      // Both tangent and residual are linearized at the same iterate.
      auto ivw = Solid::InternalVirtualWork(law, uCurrent);

      Problem newtonProblem(du, v);
      if (residualSign == ResidualSign::Correct)
      {
        newtonProblem = ivw(du, v)
                      + DirichletBC(du, zero);
      }
      else
      {
        newtonProblem = ivw.Tangent(du, v)
                      - ivw.Residual(v)
                      + DirichletBC(du, zero);
      }

      SparseLU linearSolver(newtonProblem);
      NewtonSolver newton(linearSolver);
      if (residualSign == ResidualSign::Correct)
      {
        newton.setMaxIterations(12)
          .setAbsoluteTolerance(1e-12)
          .setRelativeTolerance(1e-10);
      }
      else
      {
        // A wrong residual sign turns Newton's descent direction into an
        // ascent direction for the residual norm, so the residual grows
        // instead of shrinking. Take a couple of small damped steps so every
        // assembled iterate stays physical (det F > 0, see
        // KinematicState::setDisplacementGradient) while the regression is
        // still clearly observable. Running undamped to divergence would drive
        // an element through inversion and trip the det F > 0 invariant.
        newton.setMaxIterations(2)
          .setDampingFactor(0.25)
          .setAbsoluteTolerance(1e-12)
          .setRelativeTolerance(1e-10);
      }
      newton.solve(uCurrent);

      P1 scalar(mesh);
      GridFunction err2(scalar);
      err2 = Pow(Frobenius(uCurrent - uExact), 2);
      return {Integral(err2).compute(), newton.converged(),
        newton.getReport().initialResidual, newton.getReport().finalResidual};
    }
  }

  /// @brief Verifies hyper elasticity neo hookean affine newton solver for manufactured P1 by checking tolerance-based numerical results, true predicates, manufactured-solution convergence.
  TEST(Rodin_Manufactured_P1, HyperElasticity_NeoHookean_Affine_NewtonSolver)
  {
    const auto result = solveAffineHyperElasticity(ResidualSign::Correct);
    EXPECT_TRUE(result.converged);
    EXPECT_LT(result.finalResidual, result.initialResidual);
    EXPECT_NEAR(result.l2ErrorSquared, 0.0, RODIN_FUZZY_CONSTANT);
  }

  /// @brief Verifies hyper elasticity neo hookean wrong residual sign regresses for manufactured P1 by checking false predicates, manufactured-solution convergence.
  TEST(Rodin_Manufactured_P1, HyperElasticity_NeoHookean_WrongResidualSignRegresses)
  {
    const auto result = solveAffineHyperElasticity(ResidualSign::Wrong);
    EXPECT_FALSE(result.converged);
    // The wrong sign makes the residual grow rather than decay.
    EXPECT_GT(result.finalResidual, result.initialResidual);
    EXPECT_GT(result.l2ErrorSquared, 1e-6);
  }

  /// @brief Verifies the NeoHookean internal virtual work tangent against a central finite difference of the residual in 2D.
  TEST(
    Rodin_Manufactured_P1, HyperElasticity_NeoHookean_InternalVirtualWork_FDConsistency2D)
  {
    Mesh mesh = makeUnitSquareMesh();
    P1 Vh(mesh, mesh.getSpaceDimension());

    GridFunction uCurrent(Vh);
    const Real pi = Math::Constants::pi();
    uCurrent = VectorFunction{0.07 * F::x + 0.02 * sin(pi * F::x) * sin(pi * F::y),
      -0.05 * F::y + 0.01 * cos(pi * F::x) * sin(pi * F::y)};

    const auto result = checkNeoHookeanInternalVirtualWorkFD(Vh, uCurrent);
    EXPECT_LT(result.relativeError, 1e-6)
      << "absoluteError = " << result.absoluteError
      << ", tangentNorm = " << result.tangentNorm
      << ", finiteDifferenceNorm = " << result.finiteDifferenceNorm;
  }

  /// @brief Verifies the NeoHookean internal virtual work tangent against a central finite difference of the residual in 3D.
  TEST(
    Rodin_Manufactured_P1, HyperElasticity_NeoHookean_InternalVirtualWork_FDConsistency3D)
  {
    Mesh mesh = makeUnitCubeMesh();
    P1 Vh(mesh, mesh.getSpaceDimension());

    GridFunction uCurrent(Vh);
    const Real pi = Math::Constants::pi();
    uCurrent = VectorFunction{0.04 * F::x + 0.01 * sin(pi * F::x) * sin(pi * F::y),
      -0.03 * F::y + 0.01 * sin(pi * F::y) * sin(pi * F::z),
      0.02 * F::z + 0.01 * sin(pi * F::x) * sin(pi * F::z)};

    const auto result = checkNeoHookeanInternalVirtualWorkFD(Vh, uCurrent);
    EXPECT_LT(result.relativeError, 1e-6)
      << "absoluteError = " << result.absoluteError
      << ", tangentNorm = " << result.tangentNorm
      << ", finiteDifferenceNorm = " << result.finiteDifferenceNorm;
  }

  /// @brief Verifies the ActiveContraction internal virtual work tangent against a central finite difference of the residual in 2D.
  TEST(Rodin_Manufactured_P1,
    HyperElasticity_ActiveContraction_InternalVirtualWork_FDConsistency2D)
  {
    Mesh mesh = makeUnitSquareMesh();
    P1 Vh(mesh, mesh.getSpaceDimension());

    GridFunction uCurrent(Vh);
    const Real pi = Math::Constants::pi();
    uCurrent = VectorFunction{0.05 * F::x + 0.015 * sin(pi * F::x) * sin(pi * F::y),
      -0.04 * F::y + 0.012 * cos(pi * F::x) * sin(pi * F::y)};

    const auto result = checkActiveContractionInternalVirtualWorkFD(Vh, uCurrent);
    EXPECT_GT(result.finiteDifferenceNorm, 1e-10);
    EXPECT_LT(result.relativeError, 2e-5)
      << "absoluteError = " << result.absoluteError
      << ", tangentNorm = " << result.tangentNorm
      << ", finiteDifferenceNorm = " << result.finiteDifferenceNorm;
  }

  /// @brief Verifies the ActiveContraction internal virtual work tangent against a central finite difference of the residual in 3D.
  TEST(Rodin_Manufactured_P1,
    HyperElasticity_ActiveContraction_InternalVirtualWork_FDConsistency3D)
  {
    Mesh mesh = makeUnitCubeMesh();
    P1 Vh(mesh, mesh.getSpaceDimension());

    GridFunction uCurrent(Vh);
    const Real pi = Math::Constants::pi();
    uCurrent = VectorFunction{0.035 * F::x + 0.008 * sin(pi * F::x) * sin(pi * F::y),
      -0.025 * F::y + 0.007 * sin(pi * F::y) * sin(pi * F::z),
      0.020 * F::z + 0.006 * sin(pi * F::x) * sin(pi * F::z)};

    const auto result = checkActiveContractionInternalVirtualWorkFD(Vh, uCurrent);
    EXPECT_GT(result.finiteDifferenceNorm, 1e-10);
    EXPECT_LT(result.relativeError, 2e-5)
      << "absoluteError = " << result.absoluteError
      << ", tangentNorm = " << result.tangentNorm
      << ", finiteDifferenceNorm = " << result.finiteDifferenceNorm;
  }
}
