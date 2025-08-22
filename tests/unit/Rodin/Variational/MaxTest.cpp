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
  TEST(Rodin_Variational_Max, RealFunction_TwoValues)
  {
    RealFunction f1(3.0);
    RealFunction f2(7.0);
    auto max_result = Max(f1, f2);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(max_result.getValue(p), 7.0, 1e-10);
  }

  TEST(Rodin_Variational_Max, RealFunction_EqualValues)
  {
    RealFunction f1(5.0);
    RealFunction f2(5.0);
    auto max_result = Max(f1, f2);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(max_result.getValue(p), 5.0, 1e-10);
  }

  TEST(Rodin_Variational_Max, RealFunction_NegativeValues)
  {
    RealFunction f1(-2.0);
    RealFunction f2(-8.0);
    auto max_result = Max(f1, f2);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.2}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(max_result.getValue(p), -2.0, 1e-10);
  }

  TEST(Rodin_Variational_Max, RealFunction_MixedSigns)
  {
    RealFunction f1(-3.0);
    RealFunction f2(4.0);
    auto max_result = Max(f1, f2);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.9}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(max_result.getValue(p), 4.0, 1e-10);
  }

  TEST(Rodin_Variational_Max, RealFunction_WithZero)
  {
    RealFunction f1(0.0);
    RealFunction f2(-10.0);
    auto max_result = Max(f1, f2);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(max_result.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Max, RealFunction_SmallValues)
  {
    RealFunction f1(1e-6);
    RealFunction f2(1e-5);
    auto max_result = Max(f1, f2);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.8}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(max_result.getValue(p), 1e-5, 1e-16);
  }

  TEST(Rodin_Variational_Max, GridFunction_Maximum)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    P1 fes(mesh);
    GridFunction gf1(fes);
    GridFunction gf2(fes);
    
    // Project constant functions
    RealFunction f1(2.0);
    RealFunction f2(5.0);
    gf1.project(f1);
    gf2.project(f2);
    
    auto max_result = Max(gf1, gf2);
    
    // Test at a point
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.25, 0.25}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(max_result.getValue(p), 5.0, 1e-10);
  }

  TEST(Rodin_Variational_Max, ChainedMaximum)
  {
    RealFunction f1(10.0);
    RealFunction f2(5.0);
    RealFunction f3(7.0);
    auto max_result = Max(Max(f1, f2), f3);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.4, 0.6}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(max_result.getValue(p), 10.0, 1e-10);
  }

  TEST(Rodin_Variational_Max, Commutativity)
  {
    RealFunction f1(8.0);
    RealFunction f2(3.0);
    
    auto max1 = Max(f1, f2);
    auto max2 = Max(f2, f1);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.7, 0.3}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(max1.getValue(p), max2.getValue(p), 1e-10);
    EXPECT_NEAR(max1.getValue(p), 8.0, 1e-10);
  }

  TEST(Rodin_Variational_Max, Associativity)
  {
    RealFunction f1(6.0);
    RealFunction f2(2.0);
    RealFunction f3(9.0);
    
    auto max1 = Max(Max(f1, f2), f3);  // max(max(6,2), 9) = max(6, 9) = 9
    auto max2 = Max(f1, Max(f2, f3));  // max(6, max(2,9)) = max(6, 9) = 9
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.9, 0.1}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(max1.getValue(p), max2.getValue(p), 1e-10);
    EXPECT_NEAR(max1.getValue(p), 9.0, 1e-10);
  }

  TEST(Rodin_Variational_Max, Idempotence)
  {
    RealFunction f(12.0);
    auto max_result = Max(f, f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 0.0}};
    Point p(polytope, rc);
    
    // max(a, a) = a
    EXPECT_NEAR(max_result.getValue(p), 12.0, 1e-10);
  }

  TEST(Rodin_Variational_Max, LargeValues)
  {
    RealFunction f1(1e6);
    RealFunction f2(2e6);
    auto max_result = Max(f1, f2);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{1.0, 0.0}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(max_result.getValue(p), 2e6, 1e-3);
  }

  TEST(Rodin_Variational_Max, MaxWithOperations)
  {
    RealFunction f1(4.0);
    RealFunction f2(3.0);
    RealFunction f3(2.0);
    
    auto sum = f2 + f3;  // 3 + 2 = 5
    auto max_result = Max(f1, sum);  // max(4, 5) = 5
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(max_result.getValue(p), 5.0, 1e-10);
  }

  TEST(Rodin_Variational_Max, MaxMonotonicity)
  {
    RealFunction f1(3.0);
    RealFunction f2(7.0);
    RealFunction f3(5.0);  // f1 < f3
    
    auto max1 = Max(f1, f2);  // max(3, 7) = 7
    auto max2 = Max(f3, f2);  // max(5, 7) = 7
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.8}};
    Point p(polytope, rc);
    
    // Since f1 < f3, we should have max(f1, f2) <= max(f3, f2)
    EXPECT_LE(max1.getValue(p), max2.getValue(p) + 1e-10);
    EXPECT_NEAR(max1.getValue(p), 7.0, 1e-10);
    EXPECT_NEAR(max2.getValue(p), 7.0, 1e-10);
  }

  TEST(Rodin_Variational_Max, MaxWithAbsoluteValue)
  {
    RealFunction f1(-5.0);
    RealFunction f2(3.0);
    
    auto abs_f1 = Abs(f1);  // |-5| = 5
    auto max_result = Max(abs_f1, f2);  // max(5, 3) = 5
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.1}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(max_result.getValue(p), 5.0, 1e-10);
  }

  TEST(Rodin_Variational_Max, FractionalValues)
  {
    RealFunction f1(0.3);
    RealFunction f2(0.7);
    auto max_result = Max(f1, f2);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(max_result.getValue(p), 0.7, 1e-10);
  }

  TEST(Rodin_Variational_Max, MaxMinRelationship)
  {
    RealFunction f1(4.0);
    RealFunction f2(9.0);
    
    auto max_result = Max(f1, f2);  // max(4, 9) = 9
    auto min_result = Min(f1, f2);  // min(4, 9) = 4
    auto sum = max_result + min_result;  // 9 + 4 = 13
    auto direct_sum = f1 + f2;  // 4 + 9 = 13
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);
    
    // max(a,b) + min(a,b) = a + b
    EXPECT_NEAR(sum.getValue(p), direct_sum.getValue(p), 1e-10);
    EXPECT_NEAR(sum.getValue(p), 13.0, 1e-10);
  }

  TEST(Rodin_Variational_Max, MaxWithPow)
  {
    RealFunction f1(2.0);
    RealFunction f2(3.0);
    
    auto pow_f1 = Pow(f1, 2.0);  // 2^2 = 4
    auto pow_f2 = Pow(f2, 2.0);  // 3^2 = 9
    auto max_result = Max(pow_f1, pow_f2);  // max(4, 9) = 9
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 1.0}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(max_result.getValue(p), 9.0, 1e-10);
  }

  TEST(Rodin_Variational_Max, NegativeZeroComparison)
  {
    RealFunction f1(-0.0);
    RealFunction f2(0.0);
    auto max_result = Max(f1, f2);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{1.0, 1.0}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(max_result.getValue(p), 0.0, 1e-10);
  }
}