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
  TEST(Rodin_Variational_Minus, RealFunction_Subtraction)
  {
    RealFunction f1(7.0);
    RealFunction f2(3.0);
    auto minus_result = f1 - f2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(minus_result.getValue(p), 4.0, 1e-10);
  }

  TEST(Rodin_Variational_Minus, RealFunction_ZeroResult)
  {
    RealFunction f1(5.0);
    RealFunction f2(5.0);
    auto minus_result = f1 - f2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(minus_result.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Minus, RealFunction_NegativeResult)
  {
    RealFunction f1(2.0);
    RealFunction f2(8.0);
    auto minus_result = f1 - f2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.2}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(minus_result.getValue(p), -6.0, 1e-10);
  }

  TEST(Rodin_Variational_Minus, RealFunction_SubtractZero)
  {
    RealFunction f1(42.0);
    RealFunction f2(0.0);
    auto minus_result = f1 - f2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.9}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(minus_result.getValue(p), 42.0, 1e-10);
  }

  TEST(Rodin_Variational_Minus, RealFunction_NegativeNumbers)
  {
    RealFunction f1(-3.0);
    RealFunction f2(-7.0);
    auto minus_result = f1 - f2;  // -3 - (-7) = -3 + 7 = 4
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(minus_result.getValue(p), 4.0, 1e-10);
  }

  TEST(Rodin_Variational_Minus, RealFunction_MixedSigns)
  {
    RealFunction f1(10.0);
    RealFunction f2(-4.0);
    auto minus_result = f1 - f2;  // 10 - (-4) = 10 + 4 = 14
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.8}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(minus_result.getValue(p), 14.0, 1e-10);
  }

  TEST(Rodin_Variational_Minus, VectorFunction_Subtraction)
  {
    VectorFunction vf1{5.0, 8.0};
    VectorFunction vf2{2.0, 3.0};
    auto minus_result = vf1 - vf2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.7, 0.3}};
    Point p(polytope, rc);
    
    auto result = minus_result.getValue(p);
    EXPECT_NEAR(result(0), 3.0, 1e-10);
    EXPECT_NEAR(result(1), 5.0, 1e-10);
  }

  TEST(Rodin_Variational_Minus, VectorFunction_ZeroVector)
  {
    VectorFunction vf1{3.0, -2.0, 7.0};
    VectorFunction vf2{0.0, 0.0, 0.0};
    auto minus_result = vf1 - vf2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.9, 0.1}};
    Point p(polytope, rc);
    
    auto result = minus_result.getValue(p);
    EXPECT_NEAR(result(0), 3.0, 1e-10);
    EXPECT_NEAR(result(1), -2.0, 1e-10);
    EXPECT_NEAR(result(2), 7.0, 1e-10);
  }

  TEST(Rodin_Variational_Minus, GridFunction_Subtraction)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    P1 fes(mesh);
    GridFunction gf1(fes);
    GridFunction gf2(fes);
    
    // Project constant functions
    RealFunction f1(9.0);
    RealFunction f2(4.0);
    gf1.project(f1);
    gf2.project(f2);
    
    auto minus_result = gf1 - gf2;
    
    // Test at a point
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.25, 0.25}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(minus_result.getValue(p), 5.0, 1e-10);
  }

  TEST(Rodin_Variational_Minus, ChainedSubtraction)
  {
    RealFunction f1(20.0);
    RealFunction f2(8.0);
    RealFunction f3(5.0);
    auto minus_result = f1 - f2 - f3;  // 20 - 8 - 5 = 7
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.4, 0.6}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(minus_result.getValue(p), 7.0, 1e-10);
  }

  TEST(Rodin_Variational_Minus, SubtractionWithAddition)
  {
    RealFunction f1(12.0);
    RealFunction f2(5.0);
    RealFunction f3(3.0);
    
    auto result1 = f1 - f2 + f3;  // 12 - 5 + 3 = 10
    auto result2 = f1 + f3 - f2;  // 12 + 3 - 5 = 10
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 0.0}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(result1.getValue(p), result2.getValue(p), 1e-10);
    EXPECT_NEAR(result1.getValue(p), 10.0, 1e-10);
  }

  TEST(Rodin_Variational_Minus, SubtractionSelfInverse)
  {
    RealFunction f1(15.0);
    RealFunction f2(7.0);
    
    auto minus_result = f1 - f2;  // 15 - 7 = 8
    auto plus_back = minus_result + f2;  // 8 + 7 = 15
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{1.0, 0.0}};
    Point p(polytope, rc);
    
    // (a - b) + b = a
    EXPECT_NEAR(plus_back.getValue(p), f1.getValue(p), 1e-10);
  }

  TEST(Rodin_Variational_Minus, SubtractionNonCommutative)
  {
    RealFunction f1(10.0);
    RealFunction f2(6.0);
    
    auto minus1 = f1 - f2;  // 10 - 6 = 4
    auto minus2 = f2 - f1;  // 6 - 10 = -4
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);
    
    // a - b ≠ b - a (in general)
    EXPECT_NEAR(minus1.getValue(p), 4.0, 1e-10);
    EXPECT_NEAR(minus2.getValue(p), -4.0, 1e-10);
    EXPECT_NEAR(minus1.getValue(p), -minus2.getValue(p), 1e-10);
  }

  TEST(Rodin_Variational_Minus, SubtractionWithMultiplication)
  {
    RealFunction f1(8.0);
    RealFunction f2(3.0);
    RealFunction f3(2.0);
    
    auto result1 = (f1 - f2) * f3;  // (8-3)*2 = 5*2 = 10
    auto result2 = f1 * f3 - f2 * f3;  // 8*2 - 3*2 = 16 - 6 = 10
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.8}};
    Point p(polytope, rc);
    
    // (a - b) * c = a*c - b*c
    EXPECT_NEAR(result1.getValue(p), result2.getValue(p), 1e-10);
    EXPECT_NEAR(result1.getValue(p), 10.0, 1e-10);
  }

  TEST(Rodin_Variational_Minus, SubtractionWithDivision)
  {
    RealFunction f1(15.0);
    RealFunction f2(9.0);
    RealFunction f3(3.0);
    
    auto result1 = (f1 - f2) / f3;  // (15-9)/3 = 6/3 = 2
    auto result2 = f1 / f3 - f2 / f3;  // 15/3 - 9/3 = 5 - 3 = 2
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.1}};
    Point p(polytope, rc);
    
    // (a - b) / c = a/c - b/c
    EXPECT_NEAR(result1.getValue(p), result2.getValue(p), 1e-10);
    EXPECT_NEAR(result1.getValue(p), 2.0, 1e-10);
  }

  TEST(Rodin_Variational_Minus, LargeValues)
  {
    RealFunction f1(1e6);
    RealFunction f2(3e5);
    auto minus_result = f1 - f2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(minus_result.getValue(p), 7e5, 1e-3);
  }

  TEST(Rodin_Variational_Minus, SmallValues)
  {
    RealFunction f1(1e-3);
    RealFunction f2(5e-4);
    auto minus_result = f1 - f2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(minus_result.getValue(p), 5e-4, 1e-15);
  }

  TEST(Rodin_Variational_Minus, SubtractFromSelf)
  {
    RealFunction f(25.0);
    auto minus_result = f - f;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 1.0}};
    Point p(polytope, rc);
    
    // a - a = 0
    EXPECT_NEAR(minus_result.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Minus, SubtractionWithSqrt)
  {
    RealFunction f1(25.0);
    RealFunction f2(9.0);
    auto sqrt_f1 = Sqrt(f1);   // sqrt(25) = 5
    auto sqrt_f2 = Sqrt(f2);   // sqrt(9) = 3
    auto minus_result = sqrt_f1 - sqrt_f2;  // 5 - 3 = 2
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.7, 0.7}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(minus_result.getValue(p), 2.0, 1e-10);
  }

  TEST(Rodin_Variational_Minus, SubtractionWithAbs)
  {
    RealFunction f1(3.0);
    RealFunction f2(-7.0);
    auto abs_f2 = Abs(f2);     // |-7| = 7
    auto minus_result = f1 - abs_f2;  // 3 - 7 = -4
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.2}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(minus_result.getValue(p), -4.0, 1e-10);
  }

  TEST(Rodin_Variational_Minus, VectorComponentSubtraction)
  {
    VectorFunction vf1{10.0, 15.0, 20.0};
    VectorFunction vf2{4.0, 6.0, 8.0};
    
    auto x_diff = vf1.x() - vf2.x();  // 10 - 4 = 6
    auto y_diff = vf1.y() - vf2.y();  // 15 - 6 = 9
    auto z_diff = vf1.z() - vf2.z();  // 20 - 8 = 12
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.9, 0.9}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(x_diff.getValue(p), 6.0, 1e-10);
    EXPECT_NEAR(y_diff.getValue(p), 9.0, 1e-10);
    EXPECT_NEAR(z_diff.getValue(p), 12.0, 1e-10);
  }

  TEST(Rodin_Variational_Minus, SubtractionFractionalValues)
  {
    RealFunction f1(2.75);
    RealFunction f2(1.25);
    auto minus_result = f1 - f2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.4, 0.4}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(minus_result.getValue(p), 1.5, 1e-10);
  }
}