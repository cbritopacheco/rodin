/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

namespace Rodin::Tests::Manufactured::H1Poisson
{
  template <size_t M>
  class Manufactured_Poisson_H1_Test : public ::testing::TestWithParam<Polytope::Type>
  {
    protected:
      Mesh<Context::Local> getMesh()
      {
        Mesh mesh;
        mesh = mesh.UniformGrid(GetParam(), { M, M });
        mesh.scale(1.0 / (M - 1));
        mesh.getConnectivity().compute(1, 2);
        return mesh;
      }
  };

  using Manufactured_Poisson_H1_Test_16x16 =
    Rodin::Tests::Manufactured::H1Poisson::Manufactured_Poisson_H1_Test<16>;

  constexpr auto kQuadraticOrder = std::integral_constant<size_t, 2>{};

  TEST_P(Manufactured_Poisson_H1_Test_16x16, Poisson_SimpleSine_H1)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    H1 vh(kQuadraticOrder, mesh);

    TrialFunction u(vh);
    TestFunction  v(vh);

    auto f = 2 * pi * pi * sin(pi * F::x) * sin(pi * F::y);

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, Zero());
    CG(poisson).solve();

    auto solution = sin(pi * F::x) * sin(pi * F::y);

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  TEST_P(Manufactured_Poisson_H1_Test_16x16, Poisson_NonhomogeneousDirichlet_H1)
  {
    auto pi = Rodin::Math::Constants::pi();

    Mesh mesh = this->getMesh();

    H1 vh(kQuadraticOrder, mesh);

    TrialFunction u(vh);
    TestFunction  v(vh);

    auto solution = cos(pi * F::x) * cos(pi * F::y);
    auto f = 2 * pi * pi * solution;

    Problem poisson(u, v);
    poisson = Integral(Grad(u), Grad(v))
            - Integral(f, v)
            + DirichletBC(u, solution);
    CG(poisson).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - solution, 2);

    Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  INSTANTIATE_TEST_SUITE_P(
    MeshParams16x16,
    Manufactured_Poisson_H1_Test_16x16,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );
}
