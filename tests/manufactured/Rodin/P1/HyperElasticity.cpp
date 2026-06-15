/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
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
      newton.setMaxIterations(12)
        .setAbsoluteTolerance(1e-12)
        .setRelativeTolerance(1e-10);
      newton.solve(uCurrent);

      P1 scalar(mesh);
      GridFunction err2(scalar);
      err2 = Pow(Frobenius(uCurrent - uExact), 2);
      return {
        Integral(err2).compute(),
        newton.converged(),
        newton.getReport().initial_residual,
        newton.getReport().final_residual
      };
    }
  }

  TEST(Rodin_Manufactured_P1, HyperElasticity_NeoHookean_Affine_NewtonSolver)
  {
    const auto result = solveAffineHyperElasticity(ResidualSign::Correct);
    EXPECT_TRUE(result.converged);
    EXPECT_LT(result.finalResidual, result.initialResidual);
    EXPECT_NEAR(result.l2ErrorSquared, 0.0, RODIN_FUZZY_CONSTANT);
  }

  TEST(Rodin_Manufactured_P1, HyperElasticity_NeoHookean_WrongResidualSignRegresses)
  {
    const auto result = solveAffineHyperElasticity(ResidualSign::Wrong);
    EXPECT_FALSE(result.converged);
    EXPECT_GT(result.l2ErrorSquared, 1e-6);
  }
}
