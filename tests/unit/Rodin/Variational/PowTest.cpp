#include <gtest/gtest.h>
#include <cmath>
#include <limits>
#include "Rodin/Test/Random.h"

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

namespace Rodin::Tests::Unit
{
  /// @brief Polynomial degrees handle exponent types and overflow without evaluating huge powers.
  TEST(Rodin_Variational_Pow, PolynomialDegreeArithmeticAndOverflowBoundaries)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {2});
    P1 space(mesh);
    GridFunction field(space);
    field = RealFunction(1);
    const auto cell = mesh.getCell();
    const auto check = [&](auto exponent, Optional<size_t> expected) {
      const auto expression = Pow(field, exponent);
      EXPECT_EQ(expression.getOrder(*cell), expected);
      EXPECT_EQ(expression.getExponent(), exponent);
    };
    check(std::numeric_limits<Integer>::min(), std::nullopt);
    check(Integer(-1), std::nullopt);
    check(Integer(0), size_t(0));
    check(Integer(1), size_t(1));
    check(Integer(3), size_t(3));
    check(std::numeric_limits<Integer>::max(),
      static_cast<size_t>(std::numeric_limits<Integer>::max()));
    check(size_t(0), size_t(0));
    check(size_t(1), size_t(1));
    check(size_t(3), size_t(3));
    check(false, size_t(0));
    check(true, size_t(1));
    check(Real(0), std::nullopt);
    check(Real(2), std::nullopt);
    check(Real(0.5), std::nullopt);

    constexpr size_t limit = std::numeric_limits<size_t>::max();
    const auto maximumDegree = Pow(field, limit);
    EXPECT_EQ(maximumDegree.getOrder(*cell), Optional<size_t>(limit));
    EXPECT_EQ(Pow(maximumDegree, 0).getOrder(*cell), Optional<size_t>(0));
    EXPECT_EQ(Pow(maximumDegree, 1).getOrder(*cell), Optional<size_t>(limit));
    EXPECT_FALSE(Pow(maximumDegree, 2).getOrder(*cell));
    const auto belowOverflow = Pow(field, limit / 3);
    const auto aboveOverflow = Pow(field, limit / 3 + 1);
    EXPECT_EQ(Pow(belowOverflow, 3u).getOrder(*cell), Optional<size_t>((limit / 3) * 3));
    EXPECT_FALSE(Pow(aboveOverflow, 3u).getOrder(*cell));

    const RealFunction unknown([](const Point&) { return Real(1); });
    EXPECT_FALSE(Pow(unknown, 2).getOrder(*cell));
    EXPECT_FALSE(Pow(unknown, 0).getOrder(*cell));
  }

  /// @brief Automatic quadrature integrates integer powers of an affine field exactly.
  TEST(Rodin_Variational_Pow, PolynomialDegreeControlsAutomaticQuadrature)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    mesh.scale(Real(0.5));
    mesh.getConnectivity().compute(0, 1);
    P1 space(mesh);
    GridFunction field(space);
    field = RealFunction([](const Point& point) { return 1 + point.x(); });
    P0 tests(mesh);
    TestFunction v(tests);
    for (const Integer exponent : {0, 1, 2, 3, 6, 9})
    {
      SCOPED_TRACE(exponent);
      const auto expression = Variational::pow(field, exponent);
      for (auto cell = mesh.getCell(); cell; ++cell)
        EXPECT_EQ(expression.getOrder(*cell), Optional<size_t>(exponent));
      LinearForm load(v);
      load = Integral(expression, v);
      load.assemble();
      const Real exact = (std::ldexp(Real(1), exponent + 1) - 1) / (exponent + 1);
      EXPECT_NEAR(load.getVector().sum(), exact, 5e-12 * std::max<Real>(1, exact));
    }
  }

  /// @brief Verifies real function square root for variational pow by checking tolerance-based numerical results.
  TEST(Rodin_Variational_Pow, RealFunction_SquareRoot)
  {
    RealFunction f(4.0);
    auto pow_result = Pow(f, 0.5);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 2.0, 1e-10);
  }

  /// @brief Verifies real function square for variational pow by checking tolerance-based numerical results.
  TEST(Rodin_Variational_Pow, RealFunction_Square)
  {
    RealFunction f(3.0);
    auto pow_result = Pow(f, 2.0);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 9.0, 1e-10);
  }

  /// @brief Verifies real function cube for variational pow by checking tolerance-based numerical results.
  TEST(Rodin_Variational_Pow, RealFunction_Cube)
  {
    RealFunction f(2.0);
    auto pow_result = Pow(f, 3.0);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.2}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 8.0, 1e-10);
  }

  /// @brief Verifies real function zero power for variational pow by checking tolerance-based numerical results.
  TEST(Rodin_Variational_Pow, RealFunction_ZeroPower)
  {
    RealFunction f(5.0);
    auto pow_result = Pow(f, 0.0);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.9}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 1.0, 1e-10);
  }

  /// @brief Verifies real function power of one for variational pow by checking tolerance-based numerical results.
  TEST(Rodin_Variational_Pow, RealFunction_PowerOfOne)
  {
    RealFunction f(7.0);
    auto pow_result = Pow(f, 1.0);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 7.0, 1e-10);
  }

  /// @brief Verifies real function negative power for variational pow by checking tolerance-based numerical results.
  TEST(Rodin_Variational_Pow, RealFunction_NegativePower)
  {
    RealFunction f(2.0);
    auto pow_result = Pow(f, -1.0);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.8}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 0.5, 1e-10);
  }

  /// @brief Verifies real function fractional power for variational pow by checking tolerance-based numerical results.
  TEST(Rodin_Variational_Pow, RealFunction_FractionalPower)
  {
    RealFunction f(8.0);
    auto pow_result = Pow(f, 1.0 / 3.0);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.7, 0.3}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 2.0, 1e-10);
  }

  /// @brief Verifies grid function power operation for variational pow by checking tolerance-based numerical results, grid-function projection.
  TEST(Rodin_Variational_Pow, GridFunction_PowerOperation)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
    P1 fes(mesh);
    GridFunction gf(fes);

    // Project a constant function
    RealFunction f(3.0);
    gf.project(f);

    auto pow_result = Pow(gf, 2.0);

    // Test at a point
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.25, 0.25}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 9.0, 1e-10);
  }

  /// @brief Verifies chained power operations for variational pow by checking tolerance-based numerical results.
  TEST(Rodin_Variational_Pow, ChainedPowerOperations)
  {
    RealFunction f(2.0);
    auto pow1 = Pow(f, 2.0); // f^2 = 4
    auto pow2 = Pow(pow1, 0.5); // (f^2)^0.5 = f = 2

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.4, 0.6}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow2.getValue(p), 2.0, 1e-10);
  }

  /// @brief Verifies power of zero for variational pow by checking tolerance-based numerical results.
  TEST(Rodin_Variational_Pow, PowerOfZero)
  {
    RealFunction f(0.0);
    auto pow_result = Pow(f, 2.0);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.9, 0.1}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 0.0, 1e-10);
  }

  /// @brief Verifies power of one for variational pow by checking tolerance-based numerical results.
  TEST(Rodin_Variational_Pow, PowerOfOne)
  {
    RealFunction f(1.0);
    auto pow_result = Pow(f, 100.0);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 0.0}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 1.0, 1e-10);
  }

  /// @brief Verifies large exponent for variational pow by checking tolerance-based numerical results.
  TEST(Rodin_Variational_Pow, LargeExponent)
  {
    RealFunction f(2.0);
    auto pow_result = Pow(f, 10.0);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 1024.0, 1e-10);
  }

  /// @brief Verifies negative base power for variational pow by checking tolerance-based numerical results.
  TEST(Rodin_Variational_Pow, NegativeBasePower)
  {
    RealFunction f(-2.0);
    auto pow_result = Pow(f, 3.0);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.1}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), -8.0, 1e-10);
  }

  /// @brief Verifies even power negative base for variational pow by checking tolerance-based numerical results.
  TEST(Rodin_Variational_Pow, EvenPowerNegativeBase)
  {
    RealFunction f(-3.0);
    auto pow_result = Pow(f, 2.0);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.8}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 9.0, 1e-10);
  }

  /// @brief Verifies power law validation for variational pow by checking tolerance-based numerical results.
  TEST(Rodin_Variational_Pow, PowerLawValidation)
  {
    RealFunction f(3.0);

    // Test: (a^m)^n = a^(m*n)
    auto pow_mn = Pow(f, 2.0 * 3.0); // f^6
    auto pow_m_pow_n = Pow(Pow(f, 2.0), 3.0); // (f^2)^3 = f^6

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_mn.getValue(p), pow_m_pow_n.getValue(p), 1e-9);
    EXPECT_NEAR(pow_mn.getValue(p), 729.0, 1e-10); // 3^6 = 729
  }

  /// @brief Verifies integer exponent for variational pow by checking tolerance-based numerical results.
  TEST(Rodin_Variational_Pow, IntegerExponent)
  {
    RealFunction f(5.0);
    auto pow_result = Pow(f, 4); // Using integer exponent

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 625.0, 1e-10); // 5^4 = 625
  }

  /// @brief Verifies small fractional exponent for variational pow by checking tolerance-based numerical results.
  TEST(Rodin_Variational_Pow, SmallFractionalExponent)
  {
    RealFunction f(16.0);
    auto pow_result = Pow(f, 0.25); // Fourth root

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{1.0, 0.0}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 2.0, 1e-10); // 16^(1/4) = 2
  }
}
