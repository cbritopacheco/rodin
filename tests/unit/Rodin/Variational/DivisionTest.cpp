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
  TEST(Rodin_Variational_Division, RealFunction_SimpleDivision)
  {
    RealFunction f1(8.0);
    RealFunction f2(2.0);
    auto div_result = f1 / f2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(div_result.getValue(p), 4.0, 1e-10);
  }

  TEST(Rodin_Variational_Division, RealFunction_DivisionByOne)
  {
    RealFunction f1(42.0);
    RealFunction f2(1.0);
    auto div_result = f1 / f2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(div_result.getValue(p), 42.0, 1e-10);
  }

  TEST(Rodin_Variational_Division, RealFunction_FractionalResult)
  {
    RealFunction f1(3.0);
    RealFunction f2(4.0);
    auto div_result = f1 / f2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.2}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(div_result.getValue(p), 0.75, 1e-10);
  }

  TEST(Rodin_Variational_Division, RealFunction_LargeNumerator)
  {
    RealFunction f1(1000.0);
    RealFunction f2(10.0);
    auto div_result = f1 / f2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.9}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(div_result.getValue(p), 100.0, 1e-10);
  }

  TEST(Rodin_Variational_Division, RealFunction_SmallNumerator)
  {
    RealFunction f1(0.5);
    RealFunction f2(2.0);
    auto div_result = f1 / f2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(div_result.getValue(p), 0.25, 1e-10);
  }

  TEST(Rodin_Variational_Division, RealFunction_NegativeValues)
  {
    RealFunction f1(-6.0);
    RealFunction f2(3.0);
    auto div_result = f1 / f2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.8}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(div_result.getValue(p), -2.0, 1e-10);
  }

  TEST(Rodin_Variational_Division, RealFunction_BothNegative)
  {
    RealFunction f1(-12.0);
    RealFunction f2(-4.0);
    auto div_result = f1 / f2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.7, 0.3}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(div_result.getValue(p), 3.0, 1e-10);
  }

  TEST(Rodin_Variational_Division, GridFunction_Division)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    P1 fes(mesh);
    GridFunction gf1(fes);
    GridFunction gf2(fes);
    
    // Project constant functions
    RealFunction f1(15.0);
    RealFunction f2(3.0);
    gf1.project(f1);
    gf2.project(f2);
    
    auto div_result = gf1 / gf2;
    
    // Test at a point
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.25, 0.25}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(div_result.getValue(p), 5.0, 1e-10);
  }

  TEST(Rodin_Variational_Division, ChainedDivision)
  {
    RealFunction f1(24.0);
    RealFunction f2(3.0);
    RealFunction f3(2.0);
    auto div_result = (f1 / f2) / f3;  // (24/3)/2 = 8/2 = 4
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.4, 0.6}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(div_result.getValue(p), 4.0, 1e-10);
  }

  TEST(Rodin_Variational_Division, DivisionWithMultiplication)
  {
    RealFunction f1(12.0);
    RealFunction f2(3.0);
    RealFunction f3(2.0);
    
    auto div_first = (f1 / f2) * f3;  // (12/3)*2 = 4*2 = 8
    auto mult_first = f1 * (f3 / f2); // 12*(2/3) = 12*0.667 = 8
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.9, 0.1}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(div_first.getValue(p), 8.0, 1e-10);
    EXPECT_NEAR(mult_first.getValue(p), 8.0, 1e-10);
  }

  TEST(Rodin_Variational_Division, DivisionOfSum)
  {
    RealFunction f1(5.0);
    RealFunction f2(3.0);
    RealFunction f3(2.0);
    auto sum = f1 + f2;        // 5 + 3 = 8
    auto div_result = sum / f3; // 8 / 2 = 4
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 0.0}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(div_result.getValue(p), 4.0, 1e-10);
  }

  TEST(Rodin_Variational_Division, ReciprocalProperty)
  {
    RealFunction f1(7.0);
    RealFunction f2(3.0);
    
    auto div1 = f1 / f2;       // 7/3
    auto div2 = f2 / f1;       // 3/7
    auto product = div1 * div2; // (7/3) * (3/7) = 1
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{1.0, 0.0}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(product.getValue(p), 1.0, 1e-10);
  }

  TEST(Rodin_Variational_Division, DivisionByItself)
  {
    RealFunction f(15.0);
    auto div_result = f / f;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);
    
    // a/a = 1 (for a != 0)
    EXPECT_NEAR(div_result.getValue(p), 1.0, 1e-10);
  }

  TEST(Rodin_Variational_Division, LargeDivision)
  {
    RealFunction f1(1e6);
    RealFunction f2(1e3);
    auto div_result = f1 / f2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.8}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(div_result.getValue(p), 1e3, 1e-6);
  }

  TEST(Rodin_Variational_Division, SmallDivision)
  {
    RealFunction f1(1e-6);
    RealFunction f2(1e-3);
    auto div_result = f1 / f2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.1}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(div_result.getValue(p), 1e-3, 1e-15);
  }

  TEST(Rodin_Variational_Division, DivisionWithSqrt)
  {
    RealFunction f1(16.0);
    RealFunction f2(4.0);
    auto sqrt_f1 = Sqrt(f1);   // sqrt(16) = 4
    auto div_result = sqrt_f1 / f2; // 4 / 4 = 1
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(div_result.getValue(p), 1.0, 1e-10);
  }

  TEST(Rodin_Variational_Division, DivisionWithAbs)
  {
    RealFunction f1(-8.0);
    RealFunction f2(2.0);
    auto abs_f1 = Abs(f1);     // |-8| = 8
    auto div_result = abs_f1 / f2; // 8 / 2 = 4
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(div_result.getValue(p), 4.0, 1e-10);
  }

  TEST(Rodin_Variational_Division, DistributiveProperty)
  {
    RealFunction f1(6.0);
    RealFunction f2(3.0);
    RealFunction f3(2.0);
    
    auto sum = f1 + f2;           // 6 + 3 = 9
    auto div_sum = sum / f3;      // (6+3)/2 = 9/2 = 4.5
    
    auto div1 = f1 / f3;          // 6/2 = 3
    auto div2 = f2 / f3;          // 3/2 = 1.5
    auto sum_div = div1 + div2;   // 3 + 1.5 = 4.5
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 1.0}};
    Point p(polytope, rc);
    
    // (a+b)/c = a/c + b/c
    EXPECT_NEAR(div_sum.getValue(p), sum_div.getValue(p), 1e-10);
    EXPECT_NEAR(div_sum.getValue(p), 4.5, 1e-10);
  }

  TEST(Rodin_Variational_Division, DivisionOrder)
  {
    RealFunction f1(12.0);
    RealFunction f2(4.0);
    RealFunction f3(3.0);
    
    auto div1 = (f1 / f2) / f3;  // (12/4)/3 = 3/3 = 1
    auto div2 = f1 / (f2 * f3);  // 12/(4*3) = 12/12 = 1
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.7, 0.7}};
    Point p(polytope, rc);
    
    // (a/b)/c = a/(b*c)
    EXPECT_NEAR(div1.getValue(p), div2.getValue(p), 1e-10);
    EXPECT_NEAR(div1.getValue(p), 1.0, 1e-10);
  }

  TEST(Rodin_Variational_Division, FractionalDivision)
  {
    RealFunction f1(0.75);
    RealFunction f2(0.25);
    auto div_result = f1 / f2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.2}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(div_result.getValue(p), 3.0, 1e-10);
  }
}