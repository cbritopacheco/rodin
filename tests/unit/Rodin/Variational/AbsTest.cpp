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
  TEST(Rodin_Variational_Abs, RealFunction_PositiveValue)
  {
    RealFunction f(3.14);
    auto abs_result = Abs(f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(abs_result.getValue(p), 3.14, 1e-10);
  }

  TEST(Rodin_Variational_Abs, RealFunction_NegativeValue)
  {
    RealFunction f(-2.71);
    auto abs_result = Abs(f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(abs_result.getValue(p), 2.71, 1e-10);
  }

  TEST(Rodin_Variational_Abs, RealFunction_ZeroValue)
  {
    RealFunction f(0.0);
    auto abs_result = Abs(f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.2}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(abs_result.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Abs, RealFunction_LargePositiveValue)
  {
    RealFunction f(1000.0);
    auto abs_result = Abs(f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.9}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(abs_result.getValue(p), 1000.0, 1e-10);
  }

  TEST(Rodin_Variational_Abs, RealFunction_LargeNegativeValue)
  {
    RealFunction f(-999.0);
    auto abs_result = Abs(f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(abs_result.getValue(p), 999.0, 1e-10);
  }

  TEST(Rodin_Variational_Abs, RealFunction_SmallPositiveValue)
  {
    RealFunction f(1e-6);
    auto abs_result = Abs(f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.8}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(abs_result.getValue(p), 1e-6, 1e-16);
  }

  TEST(Rodin_Variational_Abs, RealFunction_SmallNegativeValue)
  {
    RealFunction f(-1e-6);
    auto abs_result = Abs(f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.7, 0.3}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(abs_result.getValue(p), 1e-6, 1e-16);
  }

  TEST(Rodin_Variational_Abs, GridFunction_AbsoluteValue)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    P1 fes(mesh);
    GridFunction gf(fes);
    
    // Project a negative constant function
    RealFunction f(-5.0);
    gf.project(f);
    
    auto abs_result = Abs(gf);
    
    // Test at a point
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.25, 0.25}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(abs_result.getValue(p), 5.0, 1e-10);
  }

  TEST(Rodin_Variational_Abs, ChainedOperations)
  {
    RealFunction f1(3.0);
    RealFunction f2(7.0);
    auto diff = f1 - f2;  // 3 - 7 = -4
    auto abs_diff = Abs(diff);  // |-4| = 4
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.4, 0.6}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(abs_diff.getValue(p), 4.0, 1e-10);
  }

  TEST(Rodin_Variational_Abs, DoubleAbsoluteValue)
  {
    RealFunction f(-42.0);
    auto abs1 = Abs(f);      // |-42| = 42
    auto abs2 = Abs(abs1);   // |42| = 42
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.9, 0.1}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(abs2.getValue(p), 42.0, 1e-10);
    // Double absolute value should equal single absolute value for any input
    EXPECT_NEAR(abs1.getValue(p), abs2.getValue(p), 1e-10);
  }

  TEST(Rodin_Variational_Abs, AbsoluteValueIdentity)
  {
    RealFunction pos_f(15.0);
    RealFunction neg_f(-15.0);
    
    auto abs_pos = Abs(pos_f);
    auto abs_neg = Abs(neg_f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 0.0}};
    Point p(polytope, rc);
    
    // |a| = |-a| for any a
    EXPECT_NEAR(abs_pos.getValue(p), abs_neg.getValue(p), 1e-10);
    EXPECT_NEAR(abs_pos.getValue(p), 15.0, 1e-10);
  }

  TEST(Rodin_Variational_Abs, TriangleInequality)
  {
    RealFunction f1(3.0);
    RealFunction f2(-5.0);
    
    auto sum = f1 + f2;           // 3 + (-5) = -2
    auto abs_sum = Abs(sum);      // |-2| = 2
    
    auto abs_f1 = Abs(f1);        // |3| = 3
    auto abs_f2 = Abs(f2);        // |-5| = 5
    auto sum_abs = abs_f1 + abs_f2; // 3 + 5 = 8
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{1.0, 0.0}};
    Point p(polytope, rc);
    
    // Triangle inequality: |a + b| <= |a| + |b|
    EXPECT_LE(abs_sum.getValue(p), sum_abs.getValue(p) + 1e-10);
    EXPECT_NEAR(abs_sum.getValue(p), 2.0, 1e-10);
    EXPECT_NEAR(sum_abs.getValue(p), 8.0, 1e-10);
  }

  TEST(Rodin_Variational_Abs, AbsoluteValueOfProduct)
  {
    RealFunction f1(-2.0);
    RealFunction f2(3.0);
    
    auto product = f1 * f2;        // (-2) * 3 = -6
    auto abs_product = Abs(product); // |-6| = 6
    
    auto abs_f1 = Abs(f1);         // |-2| = 2
    auto abs_f2 = Abs(f2);         // |3| = 3
    auto product_abs = abs_f1 * abs_f2; // 2 * 3 = 6
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);
    
    // |a * b| = |a| * |b|
    EXPECT_NEAR(abs_product.getValue(p), product_abs.getValue(p), 1e-10);
    EXPECT_NEAR(abs_product.getValue(p), 6.0, 1e-10);
  }

  TEST(Rodin_Variational_Abs, FractionalValues)
  {
    RealFunction f(-0.5);
    auto abs_result = Abs(f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.8}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(abs_result.getValue(p), 0.5, 1e-10);
  }

  TEST(Rodin_Variational_Abs, AbsoluteValueNonNegativity)
  {
    RealFunction f(-100.0);
    auto abs_result = Abs(f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.1}};
    Point p(polytope, rc);
    
    // Absolute value is always non-negative
    EXPECT_GE(abs_result.getValue(p), 0.0);
    EXPECT_NEAR(abs_result.getValue(p), 100.0, 1e-10);
  }

  TEST(Rodin_Variational_Abs, AbsoluteValueOfDifference)
  {
    RealFunction f1(10.0);
    RealFunction f2(7.0);
    
    auto diff1 = f1 - f2;         // 10 - 7 = 3
    auto diff2 = f2 - f1;         // 7 - 10 = -3
    auto abs_diff1 = Abs(diff1);  // |3| = 3
    auto abs_diff2 = Abs(diff2);  // |-3| = 3
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);
    
    // |a - b| = |b - a|
    EXPECT_NEAR(abs_diff1.getValue(p), abs_diff2.getValue(p), 1e-10);
    EXPECT_NEAR(abs_diff1.getValue(p), 3.0, 1e-10);
  }
}