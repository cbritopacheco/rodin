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
  TEST(Rodin_Variational_Sum, RealFunction_Addition)
  {
    RealFunction f1(3.0);
    RealFunction f2(7.0);
    auto sum = f1 + f2;

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    EXPECT_NEAR(sum.getValue(p), 10.0, 1e-10);
  }

  TEST(Rodin_Variational_Sum, RealFunction_ZeroAddition)
  {
    RealFunction f1(42.0);
    RealFunction f2(0.0);
    auto sum = f1 + f2;

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);

    EXPECT_NEAR(sum.getValue(p), 42.0, 1e-10);
  }

  TEST(Rodin_Variational_Sum, RealFunction_NegativeAddition)
  {
    RealFunction f1(10.0);
    RealFunction f2(-3.0);
    auto sum = f1 + f2;

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.2}};
    Point p(polytope, rc);

    EXPECT_NEAR(sum.getValue(p), 7.0, 1e-10);
  }

  TEST(Rodin_Variational_Sum, VectorFunction_Addition)
  {
    VectorFunction vf1{1.0, 2.0};
    VectorFunction vf2{3.0, 4.0};
    auto sum = vf1 + vf2;

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.9}};
    Point p(polytope, rc);

    auto result = sum.getValue(p);
    EXPECT_NEAR(result(0), 4.0, 1e-10);
    EXPECT_NEAR(result(1), 6.0, 1e-10);
  }

  TEST(Rodin_Variational_Sum, VectorFunction_ZeroVector_Addition)
  {
    VectorFunction vf1{5.0, -2.0, 8.0};
    VectorFunction vf2{0.0, 0.0, 0.0};
    auto sum = vf1 + vf2;

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);

    auto result = sum.getValue(p);
    EXPECT_NEAR(result(0), 5.0, 1e-10);
    EXPECT_NEAR(result(1), -2.0, 1e-10);
    EXPECT_NEAR(result(2), 8.0, 1e-10);
  }

  TEST(Rodin_Variational_Sum, GridFunction_Addition)
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

    auto sum = gf1 + gf2;

    // Test at a point
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.25, 0.25}};
    Point p(polytope, rc);

    EXPECT_NEAR(sum.getValue(p), 5.0, 1e-10);
  }

  TEST(Rodin_Variational_Sum, ChainedAddition)
  {
    RealFunction f1(1.0);
    RealFunction f2(2.0);
    RealFunction f3(3.0);
    auto sum = f1 + f2 + f3;

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.4, 0.6}};
    Point p(polytope, rc);

    EXPECT_NEAR(sum.getValue(p), 6.0, 1e-10);
  }

  TEST(Rodin_Variational_Sum, MixedScalarVector_Components)
  {
    RealFunction scalar1(5.0);
    RealFunction scalar2(3.0);
    VectorFunction vector1{1.0, 2.0};

    // Test scalar + scalar
    auto scalar_sum = scalar1 + scalar2;

    // Test that we can access components
    auto x_comp = vector1.x() + scalar1;
    auto y_comp = vector1.y() + scalar2;

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 1.0}};
    Point p(polytope, rc);

    EXPECT_NEAR(scalar_sum.getValue(p), 8.0, 1e-10);
    EXPECT_NEAR(x_comp.getValue(p), 6.0, 1e-10);
    EXPECT_NEAR(y_comp.getValue(p), 5.0, 1e-10);
  }

  TEST(Rodin_Variational_Sum, AssociativityTest)
  {
    RealFunction f1(10.0);
    RealFunction f2(20.0);
    RealFunction f3(30.0);

    auto sum1 = (f1 + f2) + f3;
    auto sum2 = f1 + (f2 + f3);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.7, 0.3}};
    Point p(polytope, rc);

    EXPECT_NEAR(sum1.getValue(p), sum2.getValue(p), 1e-10);
    EXPECT_NEAR(sum1.getValue(p), 60.0, 1e-10);
  }

  TEST(Rodin_Variational_Sum, CommutativityTest)
  {
    RealFunction f1(15.0);
    RealFunction f2(25.0);

    auto sum1 = f1 + f2;
    auto sum2 = f2 + f1;

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.9, 0.1}};
    Point p(polytope, rc);

    EXPECT_NEAR(sum1.getValue(p), sum2.getValue(p), 1e-10);
    EXPECT_NEAR(sum1.getValue(p), 40.0, 1e-10);
  }

  TEST(Rodin_Variational_Sum, LargeValueAddition)
  {
    RealFunction f1(1e6);
    RealFunction f2(2e6);
    auto sum = f1 + f2;

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    EXPECT_NEAR(sum.getValue(p), 3e6, 1e-4);
  }
}
