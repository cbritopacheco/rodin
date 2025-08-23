#include <gtest/gtest.h>
#include "Rodin/Test/Random.h"

#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

namespace Rodin::Tests::Unit
{
  TEST(Rodin_Variational_Pow, RealFunction_SquareRoot)
  {
    RealFunction f(4.0);
    auto pow_result = Pow(f, 0.5);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 2.0, 1e-10);
  }

  TEST(Rodin_Variational_Pow, RealFunction_Square)
  {
    RealFunction f(3.0);
    auto pow_result = Pow(f, 2.0);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 9.0, 1e-10);
  }

  TEST(Rodin_Variational_Pow, RealFunction_Cube)
  {
    RealFunction f(2.0);
    auto pow_result = Pow(f, 3.0);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.2}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 8.0, 1e-10);
  }

  TEST(Rodin_Variational_Pow, RealFunction_ZeroPower)
  {
    RealFunction f(5.0);
    auto pow_result = Pow(f, 0.0);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.9}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 1.0, 1e-10);
  }

  TEST(Rodin_Variational_Pow, RealFunction_PowerOfOne)
  {
    RealFunction f(7.0);
    auto pow_result = Pow(f, 1.0);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 7.0, 1e-10);
  }

  TEST(Rodin_Variational_Pow, RealFunction_NegativePower)
  {
    RealFunction f(2.0);
    auto pow_result = Pow(f, -1.0);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.8}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 0.5, 1e-10);
  }

  TEST(Rodin_Variational_Pow, RealFunction_FractionalPower)
  {
    RealFunction f(8.0);
    auto pow_result = Pow(f, 1.0/3.0);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.7, 0.3}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 2.0, 1e-10);
  }

  TEST(Rodin_Variational_Pow, GridFunction_PowerOperation)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
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

  TEST(Rodin_Variational_Pow, ChainedPowerOperations)
  {
    RealFunction f(2.0);
    auto pow1 = Pow(f, 2.0);  // f^2 = 4
    auto pow2 = Pow(pow1, 0.5); // (f^2)^0.5 = f = 2

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.4, 0.6}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow2.getValue(p), 2.0, 1e-10);
  }

  TEST(Rodin_Variational_Pow, PowerOfZero)
  {
    RealFunction f(0.0);
    auto pow_result = Pow(f, 2.0);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.9, 0.1}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Pow, PowerOfOne)
  {
    RealFunction f(1.0);
    auto pow_result = Pow(f, 100.0);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 0.0}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 1.0, 1e-10);
  }

  TEST(Rodin_Variational_Pow, LargeExponent)
  {
    RealFunction f(2.0);
    auto pow_result = Pow(f, 10.0);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 1024.0, 1e-10);
  }

  TEST(Rodin_Variational_Pow, NegativeBasePower)
  {
    RealFunction f(-2.0);
    auto pow_result = Pow(f, 3.0);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.1}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), -8.0, 1e-10);
  }

  TEST(Rodin_Variational_Pow, EvenPowerNegativeBase)
  {
    RealFunction f(-3.0);
    auto pow_result = Pow(f, 2.0);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.8}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 9.0, 1e-10);
  }

  TEST(Rodin_Variational_Pow, PowerLawValidation)
  {
    RealFunction f(3.0);

    // Test: (a^m)^n = a^(m*n)
    auto pow_mn = Pow(f, 2.0 * 3.0);     // f^6
    auto pow_m_pow_n = Pow(Pow(f, 2.0), 3.0); // (f^2)^3 = f^6

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_mn.getValue(p), pow_m_pow_n.getValue(p), 1e-9);
    EXPECT_NEAR(pow_mn.getValue(p), 729.0, 1e-10); // 3^6 = 729
  }

  TEST(Rodin_Variational_Pow, IntegerExponent)
  {
    RealFunction f(5.0);
    auto pow_result = Pow(f, 4);  // Using integer exponent

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 625.0, 1e-10); // 5^4 = 625
  }

  TEST(Rodin_Variational_Pow, SmallFractionalExponent)
  {
    RealFunction f(16.0);
    auto pow_result = Pow(f, 0.25);  // Fourth root

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{1.0, 0.0}};
    Point p(polytope, rc);

    EXPECT_NEAR(pow_result.getValue(p), 2.0, 1e-10); // 16^(1/4) = 2
  }
}
