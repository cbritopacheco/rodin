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
  TEST(Rodin_Variational_Min, RealFunction_TwoValues)
  {
    RealFunction f1(3.0);
    RealFunction f2(7.0);
    auto min_result = Min(f1, f2);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(min_result.getValue(p), 3.0, 1e-10);
  }

  TEST(Rodin_Variational_Min, RealFunction_EqualValues)
  {
    RealFunction f1(5.0);
    RealFunction f2(5.0);
    auto min_result = Min(f1, f2);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(min_result.getValue(p), 5.0, 1e-10);
  }

  TEST(Rodin_Variational_Min, RealFunction_NegativeValues)
  {
    RealFunction f1(-2.0);
    RealFunction f2(-8.0);
    auto min_result = Min(f1, f2);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.2}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(min_result.getValue(p), -8.0, 1e-10);
  }

  TEST(Rodin_Variational_Min, RealFunction_MixedSigns)
  {
    RealFunction f1(-3.0);
    RealFunction f2(4.0);
    auto min_result = Min(f1, f2);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.9}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(min_result.getValue(p), -3.0, 1e-10);
  }

  TEST(Rodin_Variational_Min, RealFunction_WithZero)
  {
    RealFunction f1(0.0);
    RealFunction f2(10.0);
    auto min_result = Min(f1, f2);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(min_result.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Min, RealFunction_SmallValues)
  {
    RealFunction f1(1e-6);
    RealFunction f2(1e-5);
    auto min_result = Min(f1, f2);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.8}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(min_result.getValue(p), 1e-6, 1e-16);
  }

  TEST(Rodin_Variational_Min, GridFunction_Minimum)
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
    
    auto min_result = Min(gf1, gf2);
    
    // Test at a point
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.25, 0.25}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(min_result.getValue(p), 2.0, 1e-10);
  }

  TEST(Rodin_Variational_Min, ChainedMinimum)
  {
    RealFunction f1(10.0);
    RealFunction f2(5.0);
    RealFunction f3(7.0);
    auto min_result = Min(Min(f1, f2), f3);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.4, 0.6}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(min_result.getValue(p), 5.0, 1e-10);
  }

  TEST(Rodin_Variational_Min, Commutativity)
  {
    RealFunction f1(8.0);
    RealFunction f2(3.0);
    
    auto min1 = Min(f1, f2);
    auto min2 = Min(f2, f1);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.7, 0.3}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(min1.getValue(p), min2.getValue(p), 1e-10);
    EXPECT_NEAR(min1.getValue(p), 3.0, 1e-10);
  }

  TEST(Rodin_Variational_Min, Associativity)
  {
    RealFunction f1(6.0);
    RealFunction f2(2.0);
    RealFunction f3(9.0);
    
    auto min1 = Min(Min(f1, f2), f3);  // min(min(6,2), 9) = min(2, 9) = 2
    auto min2 = Min(f1, Min(f2, f3));  // min(6, min(2,9)) = min(6, 2) = 2
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.9, 0.1}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(min1.getValue(p), min2.getValue(p), 1e-10);
    EXPECT_NEAR(min1.getValue(p), 2.0, 1e-10);
  }

  TEST(Rodin_Variational_Min, Idempotence)
  {
    RealFunction f(12.0);
    auto min_result = Min(f, f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 0.0}};
    Point p(polytope, rc);
    
    // min(a, a) = a
    EXPECT_NEAR(min_result.getValue(p), 12.0, 1e-10);
  }

  TEST(Rodin_Variational_Min, LargeValues)
  {
    RealFunction f1(1e6);
    RealFunction f2(2e6);
    auto min_result = Min(f1, f2);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{1.0, 0.0}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(min_result.getValue(p), 1e6, 1e-3);
  }

  TEST(Rodin_Variational_Min, MinWithOperations)
  {
    RealFunction f1(4.0);
    RealFunction f2(6.0);
    RealFunction f3(2.0);
    
    auto sum = f1 + f3;  // 4 + 2 = 6
    auto min_result = Min(sum, f2);  // min(6, 6) = 6
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(min_result.getValue(p), 6.0, 1e-10);
  }

  TEST(Rodin_Variational_Min, MinMonotonicity)
  {
    RealFunction f1(3.0);
    RealFunction f2(7.0);
    RealFunction f3(1.0);  // f3 < f1
    
    auto min1 = Min(f1, f2);  // min(3, 7) = 3
    auto min2 = Min(f3, f2);  // min(1, 7) = 1
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.8}};
    Point p(polytope, rc);
    
    // Since f3 < f1, we should have min(f3, f2) <= min(f1, f2)
    EXPECT_LE(min2.getValue(p), min1.getValue(p) + 1e-10);
    EXPECT_NEAR(min1.getValue(p), 3.0, 1e-10);
    EXPECT_NEAR(min2.getValue(p), 1.0, 1e-10);
  }

  TEST(Rodin_Variational_Min, MinWithAbsoluteValue)
  {
    RealFunction f1(-5.0);
    RealFunction f2(3.0);
    
    auto abs_f1 = Abs(f1);  // |-5| = 5
    auto min_result = Min(abs_f1, f2);  // min(5, 3) = 3
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.1}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(min_result.getValue(p), 3.0, 1e-10);
  }

  TEST(Rodin_Variational_Min, FractionalValues)
  {
    RealFunction f1(0.3);
    RealFunction f2(0.7);
    auto min_result = Min(f1, f2);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(min_result.getValue(p), 0.3, 1e-10);
  }
}