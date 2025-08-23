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
  TEST(Rodin_Variational_Sqrt, RealFunction_PositiveValue)
  {
    RealFunction f(4.0);
    auto sqrt_result = Sqrt(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    EXPECT_NEAR(sqrt_result.getValue(p), 2.0, 1e-10);
  }

  TEST(Rodin_Variational_Sqrt, RealFunction_Zero)
  {
    RealFunction f(0.0);
    auto sqrt_result = Sqrt(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);

    EXPECT_NEAR(sqrt_result.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Sqrt, RealFunction_One)
  {
    RealFunction f(1.0);
    auto sqrt_result = Sqrt(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.2}};
    Point p(polytope, rc);

    EXPECT_NEAR(sqrt_result.getValue(p), 1.0, 1e-10);
  }

  TEST(Rodin_Variational_Sqrt, RealFunction_PerfectSquare)
  {
    RealFunction f(9.0);
    auto sqrt_result = Sqrt(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.9}};
    Point p(polytope, rc);

    EXPECT_NEAR(sqrt_result.getValue(p), 3.0, 1e-10);
  }

  TEST(Rodin_Variational_Sqrt, RealFunction_LargeValue)
  {
    RealFunction f(100.0);
    auto sqrt_result = Sqrt(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);

    EXPECT_NEAR(sqrt_result.getValue(p), 10.0, 1e-10);
  }

  TEST(Rodin_Variational_Sqrt, RealFunction_FractionalValue)
  {
    RealFunction f(0.25);
    auto sqrt_result = Sqrt(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.8}};
    Point p(polytope, rc);

    EXPECT_NEAR(sqrt_result.getValue(p), 0.5, 1e-10);
  }

  TEST(Rodin_Variational_Sqrt, GridFunction_SquareRoot)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    P1 fes(mesh);
    GridFunction gf(fes);

    // Project a constant function
    RealFunction f(16.0);
    gf.project(f);

    auto sqrt_result = Sqrt(gf);

    // Test at a point
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.25, 0.25}};
    Point p(polytope, rc);

    EXPECT_NEAR(sqrt_result.getValue(p), 4.0, 1e-10);
  }

  TEST(Rodin_Variational_Sqrt, ChainedOperations)
  {
    RealFunction f(2.0);
    auto squared = f * f;  // f^2 = 4
    auto sqrt_result = Sqrt(squared);  // sqrt(f^2) = |f| = f (for positive f)

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.4, 0.6}};
    Point p(polytope, rc);

    EXPECT_NEAR(sqrt_result.getValue(p), 2.0, 1e-10);
  }

  TEST(Rodin_Variational_Sqrt, NestedSquareRoot)
  {
    RealFunction f(256.0);
    auto sqrt1 = sqrt(f);       // sqrt(256) = 16
    auto sqrt2 = sqrt(sqrt1);   // sqrt(16) = 4

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.7, 0.3}};
    Point p(polytope, rc);

    EXPECT_NEAR(f.getValue(p), 256.0, 1e-10);
    EXPECT_NEAR(sqrt1.getValue(p), 16.0, 1e-10);
    EXPECT_NEAR(sqrt2.getValue(p), 4.0, 1e-10);
  }

  TEST(Rodin_Variational_Sqrt, IrrationalSquareRoot)
  {
    RealFunction f(2.0);
    auto sqrt_result = Sqrt(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.9, 0.1}};
    Point p(polytope, rc);

    EXPECT_NEAR(sqrt_result.getValue(p), std::sqrt(2.0), 1e-10);
  }

  TEST(Rodin_Variational_Sqrt, SmallValue)
  {
    RealFunction f(1e-6);
    auto sqrt_result = Sqrt(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 0.0}};
    Point p(polytope, rc);

    EXPECT_NEAR(sqrt_result.getValue(p), 1e-3, 1e-10);
  }

  TEST(Rodin_Variational_Sqrt, VeryLargeValue)
  {
    RealFunction f(1e6);
    auto sqrt_result = Sqrt(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{1.0, 0.0}};
    Point p(polytope, rc);

    EXPECT_NEAR(sqrt_result.getValue(p), 1e3, 1e-6);
  }

  TEST(Rodin_Variational_Sqrt, CombinedWithOtherOperations)
  {
    RealFunction f1(3.0);
    RealFunction f2(4.0);
    auto sum = f1 + f2;  // 3 + 4 = 7
    auto sqrt_sum = Sqrt(sum);  // sqrt(7)

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    EXPECT_NEAR(sqrt_sum.getValue(p), std::sqrt(7.0), 1e-10);
  }

  TEST(Rodin_Variational_Sqrt, SqrtOfProduct)
  {
    RealFunction f1(2.0);
    RealFunction f2(8.0);
    auto product = f1 * f2;  // 2 * 8 = 16
    auto sqrt_product = Sqrt(product);  // sqrt(16) = 4

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.8}};
    Point p(polytope, rc);

    EXPECT_NEAR(sqrt_product.getValue(p), 4.0, 1e-10);
  }

  TEST(Rodin_Variational_Sqrt, MathematicalIdentities)
  {
    RealFunction f(5.0);

    // Test: sqrt(a) * sqrt(a) = a
    auto sqrt_f = Sqrt(f);
    auto sqrt_squared = sqrt_f * sqrt_f;

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.1}};
    Point p(polytope, rc);

    EXPECT_NEAR(sqrt_squared.getValue(p), 5.0, 1e-10);
  }

  TEST(Rodin_Variational_Sqrt, EquivalenceToPow)
  {
    RealFunction f(64.0);
    auto sqrt_result = Sqrt(f);
    auto pow_result = Pow(f, 0.5);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);

    EXPECT_NEAR(sqrt_result.getValue(p), pow_result.getValue(p), 1e-10);
    EXPECT_NEAR(sqrt_result.getValue(p), 8.0, 1e-10);
  }
}
