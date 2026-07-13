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
}
