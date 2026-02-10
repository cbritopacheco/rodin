/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <algorithm>
#include <type_traits>

#include "Rodin/Assembly.h"
#include "Rodin/Solver/ForwardDecls.h"
#include "Rodin/Variational.h"
#include "Rodin/Solver/CG.h"
#include "Rodin/Solver/BiCGSTAB.h"
#include "Rodin/Solver/GMRES.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Solver;

namespace Rodin::Tests::Manufactured::P1H1
{
  template <size_t NX, size_t NY, size_t NZ = 1, size_t _K = 1>
  class Manufactured_P1H1_Mixed_Test : public ::testing::TestWithParam<Polytope::Type>
  {
    public:
      static constexpr std::integral_constant<size_t, _K> K{};

      void SetUp() override
      {
        const auto geom = GetParam();
        if (geom == Polytope::Type::Tetrahedron || geom == Polytope::Type::Hexahedron || geom == Polytope::Type::Wedge)
        {
          m_mesh = Mesh().UniformGrid(geom, { NX, NY, NZ });
          m_mesh.scale(1.0 / (NX - 1));
          m_mesh.getConnectivity().compute(2, 3);
          m_mesh.getConnectivity().compute(3, 2);
          m_mesh.getConnectivity().compute(2, 1);
        }
        else
        {
          m_mesh = Mesh().UniformGrid(geom, { NX, NY });
          m_mesh.scale(1.0 / (NX - 1));
          m_mesh.getConnectivity().compute(1, 2);
          m_mesh.getConnectivity().compute(2, 1);
        }
        m_mesh.getConnectivity().compute(1, 0);
      }

      const Mesh<Context::Local>& getMesh() const { return m_mesh; }

      static auto rhs_polynomial_2d()
      {
        return RealFunction([](const Geometry::Point& p) -> double
        {
          return 2.0 * (p.x() * (1.0 - p.x()) + p.y() * p.y());
        });
      }

      static auto rhs_polynomial_3d()
      {
        return RealFunction([](const Geometry::Point& p) -> double
        {
          return 2.0 * (p.x() * (1.0 - p.x()) + p.y() * p.y() + p.z() * p.z());
        });
      }

      static auto rhs_sine_3d()
      {
        auto pi = Rodin::Math::Constants::pi();
        return RealFunction([pi](const Geometry::Point& p) -> double
        {
          return std::sin(pi * p.x()) * std::sin(pi * p.y()) * std::sin(pi * p.z());
        });
      }

    private:
      Mesh<Context::Local> m_mesh;
  };

  using Manufactured_P1H1_Mixed_Test_10x10_K1 =
    Manufactured_P1H1_Mixed_Test<10, 10, 1, 1>;
  using Manufactured_P1H1_Mixed_Test_10x10_K2 =
    Manufactured_P1H1_Mixed_Test<10, 10, 1, 2>;
  using Manufactured_P1H1_Mixed_Test_20_K1 =
    Manufactured_P1H1_Mixed_Test<20, 20, 20, 1>;
  using Manufactured_P1H1_Mixed_Test_6x6x6_K2 =
    Manufactured_P1H1_Mixed_Test<6, 6, 6, 2>;

  TEST_P(Manufactured_P1H1_Mixed_Test_10x10_K1, P1H1_Mixed_P1ExactManufacturedSolution)
  {
    const auto& mesh = getMesh();

    P1 p1h(mesh);
    H1 h1(std::integral_constant<size_t, K>{}, mesh);

    auto exact_solution = 2 * F::x + 3 * F::y + 1; // affine P1 field to exercise constant gradients
    // Mixed saddle-point system weakly enforces u ≈ p and p ≈ exact_solution; picking
    // exact_solution as the forcing makes (u, p) = (exact_solution, exact_solution)
    // the exact discrete solution when combined with the Dirichlet conditions below.

    TrialFunction u(h1);
    TrialFunction p(p1h);
    TestFunction  v(h1);
    TestFunction  q(p1h);

    Problem mixed(u, p, v, q);
    mixed = Integral(u, v)
          - Integral(p, v)
          + Integral(p, q)
          - Integral(exact_solution, q)
          + DirichletBC(u, exact_solution)
          + DirichletBC(p, exact_solution);

    BiCGSTAB(mixed).solve();

    GridFunction u_exact_coeffs(h1); // grid-function coefficients of the exact solution for u
    u_exact_coeffs = exact_solution;
    GridFunction p_exact_coeffs(p1h); // grid-function coefficients of the exact solution for p
    p_exact_coeffs = exact_solution;

    auto& A = mixed.getLinearSystem().getOperator();
    auto& b = mixed.getLinearSystem().getVector();
    auto& x = mixed.getLinearSystem().getSolution();

    auto x_exact = x;
    const auto uSize = u_exact_coeffs.getData().size();
    const auto pSize = p_exact_coeffs.getData().size();
    x_exact.head(uSize) = u_exact_coeffs.getData();
    x_exact.tail(pSize) = p_exact_coeffs.getData();

    auto r = A * x - b;
    auto re = A * x_exact - b;

    const Real scale = std::max<Real>(b.norm(), 1);
    EXPECT_NEAR(r.norm() / scale, 0, 1e-10);
    EXPECT_NEAR(re.norm() / scale, 0, 1e-12);

    GridFunction diff_u(h1);
    diff_u = Pow(u.getSolution() - exact_solution, 2);
    EXPECT_NEAR(Integral(diff_u).compute(), 0, 1e-12);

    GridFunction diff_p(p1h);
    diff_p = Pow(p.getSolution() - exact_solution, 2);
    EXPECT_NEAR(Integral(diff_p).compute(), 0, 1e-12);
  }

  template <class Fixture, class FHandle>
  void runMixedTest(const Fixture& fixture, const FHandle& rhs)
  {
    const auto& mesh = fixture.getMesh();

    P1 p1h_scalar(mesh);
    H1 h1(Fixture::K, mesh);

    TrialFunction u0(h1);
    TrialFunction u(h1);
    TestFunction  v(h1);

    TrialFunction p(p1h_scalar);
    TestFunction  q(p1h_scalar);

    const auto f = rhs();

    Problem p_l2(p, q);
    p_l2 = Integral(p, q) - Integral(f, q);
    BiCGSTAB(p_l2).solve();

    Problem u_l2(u0, v);
    u_l2 = Integral(u0, v) - Integral(p.getSolution(), v);
    BiCGSTAB(u_l2).solve();

    Problem mixed(u, p, v, q);
    mixed = Integral(u, v)
          - Integral(p, v)
          + Integral(p, q)
          - Integral(f, q);

    GMRES(mixed).solve();

    GridFunction diff(h1);
    diff = Pow(u.getSolution() - u0.getSolution(), 2);
    const Real error = Integral(diff).compute();
    EXPECT_NEAR(0.1 * error, 0, RODIN_FUZZY_CONSTANT);
  }

  // K=1
  TEST_P(Manufactured_P1H1_Mixed_Test_10x10_K1, P1H1_Mixed_PolynomialRHS)
  {
    runMixedTest(*this, rhs_polynomial_2d);
  }

  TEST_P(Manufactured_P1H1_Mixed_Test_20_K1, P1H1_Mixed_PolynomialRHS_3D)
  {
    runMixedTest(*this, rhs_polynomial_3d);
  }

  TEST_P(Manufactured_P1H1_Mixed_Test_20_K1, P1H1_Mixed_SineRHS_3D)
  {
    runMixedTest(*this, rhs_sine_3d);
  }

  // K=2
  TEST_P(Manufactured_P1H1_Mixed_Test_10x10_K2, P1H1_Mixed_PolynomialRHS_K2)
  {
    runMixedTest(*this, rhs_polynomial_2d);
  }

  TEST_P(Manufactured_P1H1_Mixed_Test_6x6x6_K2, P1H1_Mixed_PolynomialRHS_3D_K2)
  {
    runMixedTest(*this, rhs_polynomial_2d);
  }

  TEST_P(Manufactured_P1H1_Mixed_Test_6x6x6_K2, P1H1_Mixed_SineRHS_3D_K2)
  {
    runMixedTest(*this, rhs_sine_3d);
  }

  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage2D_K1,
    Manufactured_P1H1_Mixed_Test_10x10_K1,
    ::testing::Values(
      Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral)
  );

  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage3D_K1,
    Manufactured_P1H1_Mixed_Test_20_K1,
    ::testing::Values(
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron,
      Polytope::Type::Wedge
      )
  );

  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage2D_K2,
    Manufactured_P1H1_Mixed_Test_10x10_K2,
    ::testing::Values(
      Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral)
  );

  INSTANTIATE_TEST_SUITE_P(
    PolytopeCoverage3D_K2,
    Manufactured_P1H1_Mixed_Test_6x6x6_K2,
    ::testing::Values(
      Polytope::Type::Tetrahedron,
      Polytope::Type::Hexahedron,
      Polytope::Type::Wedge)
  );
}
