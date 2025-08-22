/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

namespace Rodin::Tests::Manufactured::Helmholtz
{
  template <size_t M>
  class Manufactured_Helmholtz_Test : public ::testing::TestWithParam<Polytope::Type>
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

  using Manufactured_Helmholtz_Test_32x32 = Manufactured_Helmholtz_Test<32>;
  using Manufactured_Helmholtz_Test_64x64 = Manufactured_Helmholtz_Test<64>;

  TEST_P(Manufactured_Helmholtz_Test_32x32, Helmholtz_SimpleSine)
  {
    auto pi = Math::Constants::pi();
    const Real kappa = 3.0; // wavenumber

    Mesh mesh = this->getMesh();
    P1 vh(mesh);

    auto u_expr = sin(pi * F::x) * sin(pi * F::y);
    auto f = (2 * pi * pi - kappa * kappa) * u_expr;

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem helmholtz(u, v);
    helmholtz = Integral(Grad(u), Grad(v))
              - kappa * kappa * Integral(u, v)
              - Integral(f, v)
              + DirichletBC(u, Zero());
    CG(helmholtz).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - u_expr, 2);
    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  TEST_P(Manufactured_Helmholtz_Test_64x64, Helmholtz_VariableFrequency)
  {
    auto pi = Math::Constants::pi();
    const Real kappa = 2.5;

    const Real omega1 = 3;
    const Real omega2 = 2;

    Mesh mesh = this->getMesh();
    P1 vh(mesh);

    auto u_expr = sin(omega1 * pi * F::x) * sin(omega2 * pi * F::y);
    auto f = ( (omega1 * omega1 + omega2 * omega2) * pi * pi - kappa * kappa ) * u_expr;

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem helmholtz(u, v);
    helmholtz = Integral(Grad(u), Grad(v))
              - kappa * kappa * Integral(u, v)
              - Integral(f, v)
              + DirichletBC(u, Zero());
    CG(helmholtz).solve();

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - u_expr, 2);
    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  INSTANTIATE_TEST_SUITE_P(
    MeshParams32x32,
    Manufactured_Helmholtz_Test_32x32,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );

  INSTANTIATE_TEST_SUITE_P(
    MeshParams64x64,
    Manufactured_Helmholtz_Test_64x64,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );
}

