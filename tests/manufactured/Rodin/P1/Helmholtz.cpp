/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Helmholtz manufactured solution tests.
 *
 * These tests assemble Rodin variational forms for a Helmholtz manufactured solution, solve the problem on the configured mesh, and compare against analytic fields or expected residual/error behavior. They protect the P1 finite-element and solver path, including boundary-condition handling, geometry coverage, and numerical accuracy of the manufactured workflow.
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

  /// @brief Helper used by the manufactured tests to Manufactured Helmholtz Test 32 x 32.
  using Manufactured_Helmholtz_Test_32x32 = Manufactured_Helmholtz_Test<32>;
  /// @brief Helper used by the manufactured tests to Manufactured Helmholtz Test 64 x 64.
  using Manufactured_Helmholtz_Test_64x64 = Manufactured_Helmholtz_Test<64>;

  /// @brief Verifies helmholtz P1 exact residual for manufactured helmholtz test 32 x 32 by checking tolerance-based numerical results, solver behavior, manufactured-solution convergence.
  TEST_P(Manufactured_Helmholtz_Test_32x32, Helmholtz_P1ExactResidual)
  {
    const Real kappa = 1.5;

    Mesh mesh = this->getMesh();
    P1 vh(mesh);

    auto u_expr = 2 * F::x - F::y + 1;
    auto f = -kappa * kappa * u_expr;

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem helmholtz(u, v);
    helmholtz = Integral(Grad(u), Grad(v))
              - kappa * kappa * Integral(u, v)
              - Integral(f, v)
              + DirichletBC(u, u_expr);
    CG(helmholtz).solve();

    GridFunction u_exact(vh);
    u_exact = u_expr;

    auto& A = helmholtz.getLinearSystem().getOperator();
    auto& b = helmholtz.getLinearSystem().getVector();
    auto& x = helmholtz.getLinearSystem().getSolution();

    auto r = A * x - b;
    auto re = A * u_exact.getData() - b;

    const Real scale = std::max<Real>(b.norm(), 1);
    EXPECT_NEAR(r.norm() / scale, 0, 1e-10);
    EXPECT_NEAR(re.norm() / scale, 0, 1e-12);

    GridFunction diff(vh);
    diff = Pow(u.getSolution() - u_expr, 2);
    EXPECT_NEAR(Integral(diff).compute(), 0, 1e-12);
  }

  /// @brief Verifies helmholtz simple sine for manufactured helmholtz test 32 x 32 by checking tolerance-based numerical results, solver behavior, manufactured-solution convergence.
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

  /// @brief Verifies helmholtz variable frequency for manufactured helmholtz test 64 x 64 by checking tolerance-based numerical results, solver behavior, manufactured-solution convergence.
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

  /// @brief Instantiates Manufactured Helmholtz Test 32 x 32 over the Mesh Params 32 x 32 parameter coverage.
  INSTANTIATE_TEST_SUITE_P(
    MeshParams32x32,
    Manufactured_Helmholtz_Test_32x32,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );

  /// @brief Instantiates Manufactured Helmholtz Test 64 x 64 over the Mesh Params 64 x 64 parameter coverage.
  INSTANTIATE_TEST_SUITE_P(
    MeshParams64x64,
    Manufactured_Helmholtz_Test_64x64,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );
}
