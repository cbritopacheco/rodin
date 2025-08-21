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
  TEST(Rodin_Variational_Mult, RealFunction_Multiplication)
  {
    RealFunction f1(3.0);
    RealFunction f2(7.0);
    auto product = f1 * f2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(product.getValue(p), 21.0, 1e-10);
  }

  TEST(Rodin_Variational_Mult, RealFunction_ZeroMultiplication)
  {
    RealFunction f1(42.0);
    RealFunction f2(0.0);
    auto product = f1 * f2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(product.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Mult, RealFunction_OneMultiplication)
  {
    RealFunction f1(15.0);
    RealFunction f2(1.0);
    auto product = f1 * f2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.2}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(product.getValue(p), 15.0, 1e-10);
  }

  TEST(Rodin_Variational_Mult, RealFunction_NegativeMultiplication)
  {
    RealFunction f1(4.0);
    RealFunction f2(-3.0);
    auto product = f1 * f2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.9}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(product.getValue(p), -12.0, 1e-10);
  }

  TEST(Rodin_Variational_Mult, ScalarVectorMultiplication)
  {
    RealFunction scalar(2.0);
    VectorFunction vector{3.0, 4.0, 5.0};
    auto product = scalar * vector;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);
    
    auto result = product.getValue(p);
    EXPECT_NEAR(result(0), 6.0, 1e-10);
    EXPECT_NEAR(result(1), 8.0, 1e-10);
    EXPECT_NEAR(result(2), 10.0, 1e-10);
  }

  TEST(Rodin_Variational_Mult, VectorScalarMultiplication)
  {
    VectorFunction vector{1.5, 2.5};
    RealFunction scalar(4.0);
    auto product = vector * scalar;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.8}};
    Point p(polytope, rc);
    
    auto result = product.getValue(p);
    EXPECT_NEAR(result(0), 6.0, 1e-10);
    EXPECT_NEAR(result(1), 10.0, 1e-10);
  }

  TEST(Rodin_Variational_Mult, GridFunction_Multiplication)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    P1 fes(mesh);
    GridFunction gf1(fes);
    GridFunction gf2(fes);
    
    // Project constant functions
    RealFunction f1(2.0);
    RealFunction f2(3.0);
    gf1.project(f1);
    gf2.project(f2);
    
    auto product = gf1 * gf2;
    
    // Test at a point
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.25, 0.25}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(product.getValue(p), 6.0, 1e-10);
  }

  TEST(Rodin_Variational_Mult, TrialFunction_ScalarMultiplication)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TrialFunction u(fes);
    RealFunction scalar(5.0);
    
    auto product = scalar * u;
    
    // Test that the product is constructed properly
    // Detailed mathematical testing would require assembly context
  }

  TEST(Rodin_Variational_Mult, TestFunction_ScalarMultiplication)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    P1 fes(mesh);
    TestFunction v(fes);
    RealFunction scalar(0.5);
    
    auto product = scalar * v;
    
    // Test that the product is constructed properly
    // Detailed mathematical testing would require assembly context
  }

  TEST(Rodin_Variational_Mult, ChainedMultiplication)
  {
    RealFunction f1(2.0);
    RealFunction f2(3.0);
    RealFunction f3(5.0);
    auto product = f1 * f2 * f3;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.4, 0.6}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(product.getValue(p), 30.0, 1e-10);
  }

  TEST(Rodin_Variational_Mult, AssociativityTest)
  {
    RealFunction f1(2.0);
    RealFunction f2(3.0);
    RealFunction f3(7.0);
    
    auto product1 = (f1 * f2) * f3;
    auto product2 = f1 * (f2 * f3);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.7, 0.3}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(product1.getValue(p), product2.getValue(p), 1e-10);
    EXPECT_NEAR(product1.getValue(p), 42.0, 1e-10);
  }

  TEST(Rodin_Variational_Mult, CommutativityTest)
  {
    RealFunction f1(6.0);
    RealFunction f2(9.0);
    
    auto product1 = f1 * f2;
    auto product2 = f2 * f1;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.9, 0.1}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(product1.getValue(p), product2.getValue(p), 1e-10);
    EXPECT_NEAR(product1.getValue(p), 54.0, 1e-10);
  }

  TEST(Rodin_Variational_Mult, FractionalMultiplication)
  {
    RealFunction f1(0.5);
    RealFunction f2(0.25);
    auto product = f1 * f2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(product.getValue(p), 0.125, 1e-10);
  }

  TEST(Rodin_Variational_Mult, LargeValueMultiplication)
  {
    RealFunction f1(1e3);
    RealFunction f2(2e3);
    auto product = f1 * f2;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(product.getValue(p), 2e6, 1e-3);
  }

  TEST(Rodin_Variational_Mult, VectorComponentMultiplication)
  {
    VectorFunction vector{2.0, 3.0, 4.0};
    RealFunction scalar(5.0);
    
    auto x_mult = vector.x() * scalar;
    auto y_mult = vector.y() * scalar;
    auto z_mult = vector.z() * scalar;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 1.0}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(x_mult.getValue(p), 10.0, 1e-10);
    EXPECT_NEAR(y_mult.getValue(p), 15.0, 1e-10);
    EXPECT_NEAR(z_mult.getValue(p), 20.0, 1e-10);
  }
}