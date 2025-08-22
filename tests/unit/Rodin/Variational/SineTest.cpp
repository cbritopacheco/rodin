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
  TEST(Rodin_Variational_Sin, RealFunction_Zero)
  {
    RealFunction f(0.0);
    auto sin_result = Sin(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    EXPECT_NEAR(sin_result.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Sin, RealFunction_PiOverTwo)
  {
    RealFunction f(M_PI / 2.0);
    auto sin_result = Sin(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);

    EXPECT_NEAR(sin_result.getValue(p), 1.0, 1e-10);
  }

  TEST(Rodin_Variational_Sin, RealFunction_Pi)
  {
    RealFunction f(M_PI);
    auto sin_result = Sin(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.2}};
    Point p(polytope, rc);

    EXPECT_NEAR(sin_result.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Sin, RealFunction_ThreePiOverTwo)
  {
    RealFunction f(3.0 * M_PI / 2.0);
    auto sin_result = Sin(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.9}};
    Point p(polytope, rc);

    EXPECT_NEAR(sin_result.getValue(p), -1.0, 1e-10);
  }

  TEST(Rodin_Variational_Sin, RealFunction_TwoPi)
  {
    RealFunction f(2.0 * M_PI);
    auto sin_result = Sin(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);

    EXPECT_NEAR(sin_result.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Sin, RealFunction_PiOverSix)
  {
    RealFunction f(M_PI / 6.0);
    auto sin_result = Sin(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.8}};
    Point p(polytope, rc);

    EXPECT_NEAR(sin_result.getValue(p), 0.5, 1e-10);
  }

  TEST(Rodin_Variational_Sin, RealFunction_PiOverFour)
  {
    RealFunction f(M_PI / 4.0);
    auto sin_result = Sin(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.7, 0.3}};
    Point p(polytope, rc);

    EXPECT_NEAR(sin_result.getValue(p), std::sqrt(2.0) / 2.0, 1e-10);
  }

  TEST(Rodin_Variational_Sin, RealFunction_PiOverThree)
  {
    RealFunction f(M_PI / 3.0);
    auto sin_result = Sin(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.9, 0.1}};
    Point p(polytope, rc);

    EXPECT_NEAR(sin_result.getValue(p), std::sqrt(3.0) / 2.0, 1e-10);
  }

  TEST(Rodin_Variational_Sin, RealFunction_NegativeValue)
  {
    RealFunction f(-M_PI / 6.0);
    auto sin_result = Sin(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 0.0}};
    Point p(polytope, rc);

    EXPECT_NEAR(sin_result.getValue(p), -0.5, 1e-10);
  }

  TEST(Rodin_Variational_Sin, GridFunction_SinWave)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    P1 fes(mesh);
    GridFunction gf(fes);

    // Project a constant function
    RealFunction f(M_PI / 2.0);
    gf.project(f);

    auto sin_result = Sin(gf);

    // Test at a point
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.25, 0.25}};
    Point p(polytope, rc);

    EXPECT_NEAR(sin_result.getValue(p), 1.0, 1e-10);
  }

  TEST(Rodin_Variational_Sin, ChainedOperations)
  {
    RealFunction f(M_PI / 6.0);
    auto sin_f = sin(f);          // sin(π/6) = 0.5
    auto sin_sin_f = sin(sin_f);  // sin(0.5) ≈ 0.479426

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.4, 0.6}};
    Point p(polytope, rc);

    EXPECT_NEAR(sin_sin_f.getValue(p), std::sin(0.5), 1e-10);
  }

  TEST(Rodin_Variational_Sin, Periodicity)
  {
    RealFunction f1(M_PI / 4.0);
    RealFunction f2(M_PI / 4.0 + 2.0 * M_PI);
    auto sin1 = Sin(f1);
    auto sin2 = Sin(f2);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{1.0, 0.0}};
    Point p(polytope, rc);

    // sin(x) = sin(x + 2π)
    EXPECT_NEAR(sin1.getValue(p), sin2.getValue(p), 1e-10);
  }

  TEST(Rodin_Variational_Sin, SinOfSum)
  {
    RealFunction f1(M_PI / 6.0);  // π/6
    RealFunction f2(M_PI / 3.0);  // π/3
    auto sum = f1 + f2;           // π/6 + π/3 = π/2
    auto sin_sum = Sin(sum);     // sin(π/2) = 1

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    EXPECT_NEAR(sin_sum.getValue(p), 1.0, 1e-10);
  }

  TEST(Rodin_Variational_Sin, SmallAngle)
  {
    RealFunction f(0.1);  // Small angle in radians
    auto sin_result = Sin(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.8}};
    Point p(polytope, rc);

    // For small angles, sin(x) ≈ x
    EXPECT_NEAR(sin_result.getValue(p), std::sin(0.1), 1e-10);
    EXPECT_NEAR(sin_result.getValue(p), 0.1, 1e-2);  // Approximate equality
  }

  TEST(Rodin_Variational_Sin, LargeAngle)
  {
    RealFunction f(10.0 * M_PI);
    auto sin_result = Sin(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.1}};
    Point p(polytope, rc);

    // sin(10π) = 0 (multiple of π)
    EXPECT_NEAR(sin_result.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Sin, SinSquaredIdentity)
  {
    RealFunction f(M_PI / 3.0);
    auto sin_f = Sin(f);
    auto sin_squared = sin_f * sin_f;

    // Also create cos²(x) to test sin²(x) + cos²(x) = 1 in a separate test if CoSin is available
    // For now, just test that sin²(π/3) = 3/4

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);

    // sin²(π/3) = (√3/2)² = 3/4
    EXPECT_NEAR(sin_squared.getValue(p), 3.0 / 4.0, 1e-10);
  }

  TEST(Rodin_Variational_Sin, SinWithMultiplication)
  {
    RealFunction f(M_PI / 6.0);
    RealFunction scalar(2.0);
    auto product = scalar * f;  // 2 * π/6 = π/3
    auto sin_product = Sin(product);  // sin(π/3) = √3/2

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);

    EXPECT_NEAR(sin_product.getValue(p), std::sqrt(3.0) / 2.0, 1e-10);
  }

  TEST(Rodin_Variational_Sin, SinOddFunction)
  {
    RealFunction f(M_PI / 4.0);
    RealFunction neg_f(-M_PI / 4.0);
    auto sin_f = Sin(f);
    auto sin_neg_f = Sin(neg_f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 1.0}};
    Point p(polytope, rc);

    // sin(-x) = -sin(x) (odd function)
    EXPECT_NEAR(sin_neg_f.getValue(p), -sin_f.getValue(p), 1e-10);
  }
}
