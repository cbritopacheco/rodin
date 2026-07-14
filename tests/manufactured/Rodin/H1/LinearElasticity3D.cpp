/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Linear-elasticity manufactured solution tests.
 *
 * These tests assemble Rodin variational forms for a linear-elasticity manufactured solution, solve the problem on the configured mesh, and compare against analytic fields or expected residual/error behavior. They protect the H1 finite-element and solver path, including boundary-condition handling, geometry coverage, and numerical accuracy of the manufactured workflow.
 */

#include <algorithm>
#include <gtest/gtest.h>

#include "Rodin/Assembly.h"
#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"
#include "Rodin/Solver/CG.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

namespace Rodin::Tests::Manufactured::H1LinearElasticity3D
{
  template <size_t M>
  class Manufactured_LinearElasticity3D_H1_Test
    : public ::testing::TestWithParam<Polytope::Type>
  {
  protected:
    void SetUp() override
    {
      m_mesh = Mesh().UniformGrid(GetParam(), { M, M, M });
      m_mesh.scale(1.0 / (M - 1));
      m_mesh.getConnectivity().compute(2, 3);
      m_mesh.getConnectivity().compute(3, 2);
      m_mesh.getConnectivity().compute(2, 1);
      m_mesh.getConnectivity().compute(1, 0);
    }

    const auto& getMesh() const { return m_mesh; }

    static auto polynomialSolution()
    {
      return VectorFunction{
        F::x * (1 - F::x),
        F::y * (1 - F::y),
        F::z * (1 - F::z)
      };
    }

    static auto polynomialForcing(Real lambda, Real mu)
    {
      // For u_i = coord(1 - coord):
      // Δu_i = -2, div u = 3 - 2(x+y+z), grad(div u) = (-2,-2,-2)
      // f = -mu Δu - (lambda + mu) grad(div u)
      return VectorFunction{
        2 * mu + 2 * (lambda + mu),
        2 * mu + 2 * (lambda + mu),
        2 * mu + 2 * (lambda + mu)
      };
    }

  private:
    Mesh<Context::Local> m_mesh;
  };

  /// @brief Helper used by the manufactured tests to Manufactured Linear Elasticity 3 D H1 Test 8.
  using Manufactured_LinearElasticity3D_H1_Test_8 =
    Manufactured_LinearElasticity3D_H1_Test<8>;

  /// @brief Helper used by the manufactured tests to Manufactured Linear Elasticity 3 D H1 Test 16.
  using Manufactured_LinearElasticity3D_H1_Test_16 =
    Manufactured_LinearElasticity3D_H1_Test<16>;

  /// @brief Helper used by the manufactured tests to Manufactured Linear Elasticity 3 D H1 Test 32.
  using Manufactured_LinearElasticity3D_H1_Test_32 =
    Manufactured_LinearElasticity3D_H1_Test<32>;

  /// @brief Verifies linear elasticity 3 D P1 exact residual for manufactured linear elasticity 3 D H1 test 8 by checking tolerance-based numerical results, solver behavior, manufactured-solution convergence.
  TEST_P(Manufactured_LinearElasticity3D_H1_Test_8, LinearElasticity3D_P1ExactResidual)
  {
    constexpr auto order = std::integral_constant<size_t, 1>{};
    const Real lambda = 1.0;
    const Real mu = 1.0;

    const auto& mesh = this->getMesh();

    H1 sh(order, mesh);
    H1 vh(order, mesh, mesh.getSpaceDimension());

    TrialFunction u(vh);
    TestFunction  v(vh);

    VectorFunction u_exact{ F::x, F::y, F::z };
    VectorFunction f_body{ Zero(), Zero(), Zero() };

    Problem elasticity(u, v);
    elasticity = Integral(lambda * Div(u), Div(v))
               + Integral(
                   mu * (Jacobian(u) + Jacobian(u).T()),
                   0.5 * (Jacobian(v) + Jacobian(v).T()))
               - Integral(f_body, v)
               + DirichletBC(u, u_exact);

    CG(elasticity).solve();

    GridFunction u_exact_coeffs(vh);
    u_exact_coeffs = u_exact;

    auto& A = elasticity.getLinearSystem().getOperator();
    auto& b = elasticity.getLinearSystem().getVector();
    auto& x = elasticity.getLinearSystem().getSolution();

    auto r = A * x - b;
    auto re = A * u_exact_coeffs.getData() - b;

    const Real scale = std::max<Real>(b.norm(), 1);
    EXPECT_NEAR(r.norm() / scale, 0, 1e-10);
    EXPECT_NEAR(re.norm() / scale, 0, 1e-12);

    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution() - u_exact), 2);
    EXPECT_NEAR(Integral(diff).compute(), 0, 1e-12);
  }

  /// @brief Verifies manufactured linear elasticity 3 D H1 2 for manufactured linear elasticity 3 D H1 test 8 by checking tolerance-based numerical results, solver behavior, manufactured-solution convergence.
  TEST_P(Manufactured_LinearElasticity3D_H1_Test_8, Manufactured_LinearElasticity3D_H1_2)
  {
    constexpr auto order = std::integral_constant<size_t, 2>{};
    const Real lambda = 1.5;
    const Real mu = 0.5;

    const auto& mesh = this->getMesh();

    H1 sh(order, mesh);
    H1 vh(order, mesh, mesh.getSpaceDimension());

    TrialFunction u(vh);
    TestFunction  v(vh);

    const auto u_exact = polynomialSolution();
    const auto f = polynomialForcing(lambda, mu);

    Problem elasticity(u, v);
    elasticity = Integral(lambda * Div(u), Div(v))
               + Integral(
                   mu * (Jacobian(u) + Jacobian(u).T()),
                   0.5 * (Jacobian(v) + Jacobian(v).T()))
               - Integral(f, v)
               + DirichletBC(u, u_exact);

    CG(elasticity).solve();

    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution() - u_exact), 2);

    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  template <size_t M>
  class Constant_LinearElasticity3D_H1_Test
    : public ::testing::TestWithParam<Polytope::Type>
  {
  protected:
    void SetUp() override
    {
      m_mesh = Mesh().UniformGrid(GetParam(), { M, M, M });
      m_mesh.scale(1.0 / (M - 1));
      m_mesh.getConnectivity().compute(2, 3);
      m_mesh.getConnectivity().compute(3, 2);
      m_mesh.getConnectivity().compute(2, 1);
      m_mesh.getConnectivity().compute(1, 0);
    }

    const auto& getMesh() const { return m_mesh; }

    static auto constantSolution()
    {
      // Any constant vector field: strain = 0, div = 0, so f = 0.
      return VectorFunction{ 1.0, -2.0, 0.5 };
    }

    static auto zeroForcing()
    {
      return VectorFunction{ 0.0, 0.0, 0.0 };
    }

  private:
    Mesh<Context::Local> m_mesh;
  };

  /// @brief Helper used by the manufactured tests to Constant Linear Elasticity 3 D H1 Test 8.
  using Constant_LinearElasticity3D_H1_Test_8  = Constant_LinearElasticity3D_H1_Test<8>;
  /// @brief Helper used by the manufactured tests to Constant Linear Elasticity 3 D H1 Test 16.
  using Constant_LinearElasticity3D_H1_Test_16 = Constant_LinearElasticity3D_H1_Test<16>;

  /// @brief Verifies constant solution H1 1 for constant linear elasticity 3 D H1 test 8 by checking tolerance-based numerical results, solver behavior.
  TEST_P(Constant_LinearElasticity3D_H1_Test_8, ConstantSolution_H1_1)
  {
    constexpr auto order = std::integral_constant<size_t, 1>{};
    const Real lambda = 2.0;
    const Real mu     = 1.0;

    const auto& mesh = this->getMesh();

    H1 sh(order, mesh);
    H1 vh(order, mesh, mesh.getSpaceDimension());

    TrialFunction u(vh);
    TestFunction  v(vh);

    const auto u_exact = constantSolution();
    const auto f       = zeroForcing();

    Problem elasticity(u, v);
    elasticity = Integral(lambda * Div(u), Div(v))
               + Integral(
                   mu * (Jacobian(u) + Jacobian(u).T()),
                   0.5 * (Jacobian(v) + Jacobian(v).T()))
               - Integral(f, v)
               + DirichletBC(u, u_exact);

    CG(elasticity).solve();

    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution() - u_exact), 2);

    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /// @brief Verifies constant solution H1 2 for constant linear elasticity 3 D H1 test 16 by checking tolerance-based numerical results, solver behavior.
  TEST_P(Constant_LinearElasticity3D_H1_Test_16, ConstantSolution_H1_2)
  {
    constexpr auto order = std::integral_constant<size_t, 2>{};
    const Real lambda = 1.5;
    const Real mu     = 0.5;

    const auto& mesh = this->getMesh();

    H1 sh(order, mesh);
    H1 vh(order, mesh, mesh.getSpaceDimension());

    TrialFunction u(vh);
    TestFunction  v(vh);

    const auto u_exact = constantSolution();
    const auto f       = zeroForcing();

    Problem elasticity(u, v);
    elasticity = Integral(lambda * Div(u), Div(v))
               + Integral(
                   mu * (Jacobian(u) + Jacobian(u).T()),
                   0.5 * (Jacobian(v) + Jacobian(v).T()))
               - Integral(f, v)
               + DirichletBC(u, u_exact);

    CG(elasticity).solve();

    GridFunction diff(sh);
    diff = Pow(Frobenius(u.getSolution() - u_exact), 2);

    const Real error = Integral(diff).compute();
    EXPECT_NEAR(error, 0, RODIN_FUZZY_CONSTANT);
  }

  /// @brief Instantiates Manufactured Linear Elasticity 3 D H1 Test 8 over the Polytope Coverage 3 D parameter coverage.
  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage3D,
    Manufactured_LinearElasticity3D_H1_Test_8,
    ::testing::Values(
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron,
      Polytope::Type::Pyramid,
      Polytope::Type::Wedge)
  );

  /// @brief Instantiates Constant Linear Elasticity 3 D H1 Test 8 over the Polytope Coverage 3 D parameter coverage.
  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage3D,
    Constant_LinearElasticity3D_H1_Test_8,
    ::testing::Values(
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron,
      Polytope::Type::Pyramid,
      Polytope::Type::Wedge)
  );

  /// @brief Instantiates Constant Linear Elasticity 3 D H1 Test 16 over the Polytope Coverage 3 D parameter coverage.
  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage3D,
    Constant_LinearElasticity3D_H1_Test_16,
    ::testing::Values(
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron,
      Polytope::Type::Pyramid,
      Polytope::Type::Wedge)
  );
}
