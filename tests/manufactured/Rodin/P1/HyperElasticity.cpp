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
    auto makeUnitSquareMesh()
    {
      constexpr Real meshElementSize = 1.0 / 16.0;
      Mesh mesh;
      mesh = mesh.UniformGrid(Polytope::Type::Quadrilateral, { 17, 17 });
      mesh.scale(meshElementSize);
      mesh.getConnectivity().compute(1, 2);
      return mesh;
    }
  }

  TEST(Rodin_Manufactured_P1, HyperElasticity_NeoHookean_Affine_NewtonSolver)
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
    // Affine displacement is represented exactly in P1, making this a strict
    // manufactured test for the nonlinear solver/integrator plumbing.
    auto uExact = VectorFunction{
      xDisplacementScale * F::x + xDisplacementOffset,
      yDisplacementScale * F::y + yDisplacementOffset
    };
    auto zero = VectorFunction{ Zero(), Zero() };
    uCurrent = uExact;

    TrialFunction du(Vh);
    TestFunction v(Vh);

    // Both tangent and residual are linearized at the same iterate.
    Solid::MaterialTangent tangent(law, du, v);
    tangent.setLinearizationPoint(uCurrent);

    Solid::InternalForce residual(law, v);
    residual.setLinearizationPoint(uCurrent);

    Problem newtonProblem(du, v);
    newtonProblem = tangent
                  - residual
                  + DirichletBC(du, zero);

    SparseLU linearSolver(newtonProblem);
    NewtonSolver newton(linearSolver);
    newton.setMaxIterations(20)
      .setAbsoluteTolerance(1e-12)
      .setRelativeTolerance(1e-10);
    newton.solve(uCurrent);

    P1 scalar(mesh);
    GridFunction err2(scalar);
    err2 = Pow(Frobenius(uCurrent - uExact), 2);
    const Real l2ErrorSquared = Integral(err2).compute();
    EXPECT_NEAR(l2ErrorSquared, 0.0, RODIN_FUZZY_CONSTANT);
  }
}
