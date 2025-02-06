/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Test/Random.h"
#include "Rodin/Variational.h"
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;
using namespace Rodin::Test::Random;

/**
 * @brief Manufactured solution for Poisson equation
 *
 * Miaow miaow miaow
 */
TEST(Rodin_Manufactured_P1, Poisson_1)
{
  auto pi = Rodin::Math::Constants::pi();

  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Quadrilateral, { 16, 16 });
  mesh.scale(1.0 / 15);
  mesh.getConnectivity().compute(1, 2);

  P1 vh(mesh);
  auto f = 2 * pi * sin(pi * F::x) * sin(pi * F::y);

  TrialFunction u(vh);
  TestFunction  v(vh);

  // Define problem
  Problem poisson(u, v);
  poisson = Integral(Grad(u), Grad(v))
          - Integral(f, v)
          + DirichletBC(u, Zero());

  // Solve
  CG(poisson).solve();

  auto solution = sin(pi * F::x) * sin(pi * F::y);

  GridFunction diff(vh);
  diff = Pow(u.getSolution() - solution, 2);
  diff.setWeights();

  Real error = Integral(diff).compute();

  EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
}


/**
 * @brief Manufactured solution for Poisson equation
 *
 * woof woof woof
 */
TEST(Rodin_Manufactured_P1, Poisson_2)
{
  auto pi = Rodin::Math::Constants::pi();

  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Quadrilateral, { 16, 16 });
  mesh.scale(1.0 / 15);
  mesh.getConnectivity().compute(1, 2);

  P1 vh(mesh);
  auto f = 2 * pi * sin(pi * F::x) * sin(pi * F::y);

  TrialFunction u(vh);
  TestFunction  v(vh);

  // Define problem
  Problem poisson(u, v);
  poisson = Integral(Grad(u), Grad(v))
          - Integral(f, v)
          + DirichletBC(u, Zero());

  // Solve
  CG(poisson).solve();

  auto solution = sin(pi * F::x) * sin(pi * F::y);

  GridFunction diff(vh);
  diff = Pow(u.getSolution() - solution, 2);
  diff.setWeights();

  Real error = Integral(diff).compute();

  EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
}
