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

/**
 * @brief Manufactured solutions for the divergence Poisson problem using P1 vector spaces.
 *
 * We solve the first-order problem:
 * @f[
 *   \Div(u) = f \quad 	ext{in } \Omega,
 *   \quad u = g \quad \text{on } \partial\Omega.
 * @f]
 *
 * The least-squares weak form is: Find @f$u\in V^2@f$ such that
 * @f[
 *   \int_\Omega \Div(u)\,\Div(v)\,dx = \int_\Omega f\,\Div(v)\,dx,
 * @f]
 * for all @f$v\in V^2@f$, with the essential boundary condition @f$u=g@f$.
 */
namespace Rodin::Tests::Manufactured::DivPoisson
{
  template <size_t M>
  class Manufactured_DivPoisson_Test : public ::testing::TestWithParam<Polytope::Type>
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

  using Manufactured_DivPoisson_Test_32x32 =
    Rodin::Tests::Manufactured::DivPoisson::Manufactured_DivPoisson_Test<32>;
  using Manufactured_DivPoisson_Test_64x64 =
    Rodin::Tests::Manufactured::DivPoisson::Manufactured_DivPoisson_Test<64>;

  /**
   * u(x,y) = (sin(pi x), sin(pi y))
   * Div(u) = pi cos(pi x) + pi cos(pi y)
   */
  TEST_P(Manufactured_DivPoisson_Test_32x32, Divergence_SimpleSin)
  {
    constexpr auto pi = Math::Constants::pi();
    Mesh mesh = this->getMesh();

    P1 vh(mesh, mesh.getSpaceDimension());
    auto ux = sin(pi * F::x);
    auto uy = sin(pi * F::y);
    VectorFunction u_exact{ ux, uy };
    auto f = pi * cos(pi * F::x) + pi * cos(pi * F::y);

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem divProb(u, v);
    divProb = Integral(Div(u), Div(v))
            - Integral(f, Div(v))
            + DirichletBC(u, u_exact);
    CG(divProb).solve();

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution() - u_exact), 2);

    auto error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * u(x,y) = (x, y)
   * Div(u) = 2
   */
  TEST_P(Manufactured_DivPoisson_Test_32x32, Divergence_Linear)
  {
    Mesh mesh = this->getMesh();

    P1 vh(mesh, mesh.getSpaceDimension());
    auto ux = F::x;
    auto uy = F::y;
    VectorFunction u_exact{ ux, uy };
    auto f = RealFunction(2);

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem divProb(u, v);
    divProb = Integral(Div(u), Div(v))
            - Integral(f, Div(v))
            + DirichletBC(u, u_exact);
    CG(divProb).solve();

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution() - u_exact), 2);

    auto error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /**
   * u(x,y) = (x(1-x)y, y(1-y)x)
   * Div(u) = (1-2x)y + (1-2y)x
   */
  TEST_P(Manufactured_DivPoisson_Test_32x32, Divergence_Polynomial)
  {
    Mesh mesh = this->getMesh();

    P1 vh(mesh, mesh.getSpaceDimension());
    auto ux = F::x * (1.0 - F::x) * F::y;
    auto uy = F::y * (1.0 - F::y) * F::x;
    VectorFunction u_exact{ ux, uy };
    auto f = (1.0 - 2.0 * F::x) * F::y + (1.0 - 2.0 * F::y) * F::x;

    TrialFunction u(vh);
    TestFunction  v(vh);

    Problem divProb(u, v);
    divProb = Integral(Div(u), Div(v))
            - Integral(f, Div(v))
            + DirichletBC(u, u_exact);
    CG(divProb).solve();

    P1 sh(mesh);
    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution() - u_exact), 2);

    auto error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  INSTANTIATE_TEST_SUITE_P(
    MeshParams32x32,
    Manufactured_DivPoisson_Test_32x32,
    ::testing::Values(Polytope::Type::Quadrilateral, Polytope::Type::Triangle)
  );
}

