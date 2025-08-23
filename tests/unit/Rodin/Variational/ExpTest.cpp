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
  TEST(Rodin_Variational_Exp, RealFunction_Zero)
  {
    RealFunction f(0.0);
    auto exp_result = Exp(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    EXPECT_NEAR(exp_result.getValue(p), 1.0, 1e-10);
  }

  TEST(Rodin_Variational_Exp, RealFunction_One)
  {
    RealFunction f(1.0);
    auto exp_result = Exp(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);

    EXPECT_NEAR(exp_result.getValue(p), std::exp(1.0), 1e-10);
    EXPECT_NEAR(exp_result.getValue(p), M_E, 1e-10);
  }

  TEST(Rodin_Variational_Exp, RealFunction_NegativeValue)
  {
    RealFunction f(-1.0);
    auto exp_result = Exp(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.2}};
    Point p(polytope, rc);

    EXPECT_NEAR(exp_result.getValue(p), 1.0 / M_E, 1e-10);
  }

  TEST(Rodin_Variational_Exp, RealFunction_LargeValue)
  {
    RealFunction f(2.0);
    auto exp_result = Exp(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.9}};
    Point p(polytope, rc);

    EXPECT_NEAR(exp_result.getValue(p), std::exp(2.0), 1e-10);
  }

  TEST(Rodin_Variational_Exp, RealFunction_SmallValue)
  {
    RealFunction f(0.1);
    auto exp_result = Exp(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);

    EXPECT_NEAR(exp_result.getValue(p), std::exp(0.1), 1e-10);
    // For small x, exp(x) ≈ 1 + x
    EXPECT_NEAR(exp_result.getValue(p), 1.1, 1e-2);
  }

  TEST(Rodin_Variational_Exp, RealFunction_NaturalLog)
  {
    RealFunction f(std::log(2.0));
    auto exp_result = Exp(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.8}};
    Point p(polytope, rc);

    // exp(ln(2)) = 2
    EXPECT_NEAR(exp_result.getValue(p), 2.0, 1e-10);
  }

  TEST(Rodin_Variational_Exp, GridFunction_Exponential)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    P1 fes(mesh);
    GridFunction gf(fes);

    // Project a constant function
    RealFunction f(std::log(3.0));
    gf.project(f);

    auto exp_result = Exp(gf);

    // Test at a point
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.25, 0.25}};
    Point p(polytope, rc);

    EXPECT_NEAR(exp_result.getValue(p), 3.0, 1e-10);
  }

  TEST(Rodin_Variational_Exp, ChainedOperations)
  {
    RealFunction f(0.5);
    auto exp_f = exp(f);              // exp(0.5)
    auto exp_exp_f = exp(exp_f);      // exp(exp(0.5))

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.4, 0.6}};
    Point p(polytope, rc);

    EXPECT_NEAR(exp_exp_f.getValue(p), std::exp(std::exp(0.5)), 1e-10);
  }

  TEST(Rodin_Variational_Exp, ExponentialSum)
  {
    RealFunction f1(1.0);
    RealFunction f2(2.0);
    auto sum = f1 + f2;       // 1 + 2 = 3
    auto exp_sum = Exp(sum);  // exp(3)

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.7, 0.3}};
    Point p(polytope, rc);

    EXPECT_NEAR(exp_sum.getValue(p), std::exp(3.0), 1e-10);
  }

  TEST(Rodin_Variational_Exp, ExponentialProduct)
  {
    RealFunction f1(2.0);
    RealFunction f2(0.5);
    auto product = f1 * f2;        // 2 * 0.5 = 1
    auto exp_product = Exp(product); // exp(1) = e

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.9, 0.1}};
    Point p(polytope, rc);

    EXPECT_NEAR(exp_product.getValue(p), M_E, 1e-10);
  }

  TEST(Rodin_Variational_Exp, ExponentialLaws)
  {
    RealFunction f1(2.0);
    RealFunction f2(3.0);

    auto exp_f1 = Exp(f1);         // exp(2)
    auto exp_f2 = Exp(f2);         // exp(3)
    auto product_exp = exp_f1 * exp_f2; // exp(2) * exp(3) = exp(5)

    auto sum = f1 + f2;            // 2 + 3 = 5
    auto exp_sum = Exp(sum);       // exp(5)

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 0.0}};
    Point p(polytope, rc);

    // exp(a) * exp(b) = exp(a + b)
    EXPECT_NEAR(product_exp.getValue(p), exp_sum.getValue(p), 1e-10);
  }

  TEST(Rodin_Variational_Exp, ExponentialNegativeSum)
  {
    RealFunction f1(3.0);
    RealFunction f2(-3.0);
    auto sum = f1 + f2;       // 3 + (-3) = 0
    auto exp_sum = Exp(sum);  // exp(0) = 1

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{1.0, 0.0}};
    Point p(polytope, rc);

    EXPECT_NEAR(exp_sum.getValue(p), 1.0, 1e-10);
  }

  TEST(Rodin_Variational_Exp, VerySmallValue)
  {
    RealFunction f(-10.0);
    auto exp_result = Exp(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    EXPECT_NEAR(exp_result.getValue(p), std::exp(-10.0), 1e-15);
    EXPECT_LT(exp_result.getValue(p), 1e-4);  // Should be very small
  }

  TEST(Rodin_Variational_Exp, ExponentialGrowth)
  {
    RealFunction f(5.0);
    auto exp_result = Exp(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.8}};
    Point p(polytope, rc);

    EXPECT_NEAR(exp_result.getValue(p), std::exp(5.0), 1e-9);
    EXPECT_GT(exp_result.getValue(p), 100.0);  // Should be large
  }

  TEST(Rodin_Variational_Exp, ExponentialWithSqrt)
  {
    RealFunction f(4.0);
    auto sqrt_f = Sqrt(f);      // sqrt(4) = 2
    auto exp_sqrt = Exp(sqrt_f); // exp(2)

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.1}};
    Point p(polytope, rc);

    EXPECT_NEAR(exp_sqrt.getValue(p), std::exp(2.0), 1e-10);
  }

  TEST(Rodin_Variational_Exp, ExponentialPowers)
  {
    RealFunction f(2.0);
    auto exp_f = Exp(f);         // exp(2)
    auto exp_squared = Pow(exp_f, 2.0); // (exp(2))^2 = exp(4)

    RealFunction f_double(4.0);
    auto exp_double = Exp(f_double); // exp(4)

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);

    // (exp(a))^n = exp(n*a)
    EXPECT_NEAR(exp_squared.getValue(p), exp_double.getValue(p), 1e-10);
  }

  TEST(Rodin_Variational_Exp, ExponentialMonotonicity)
  {
    RealFunction f1(1.0);
    RealFunction f2(2.0);
    auto exp1 = Exp(f1);
    auto exp2 = Exp(f2);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);

    // Since f1 < f2, we should have exp(f1) < exp(f2)
    EXPECT_LT(exp1.getValue(p), exp2.getValue(p));
  }

  TEST(Rodin_Variational_Exp, ExponentialPositivity)
  {
    RealFunction f(-100.0);  // Very negative value
    auto exp_result = Exp(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 1.0}};
    Point p(polytope, rc);

    // Exponential function is always positive
    EXPECT_GT(exp_result.getValue(p), 0.0);
  }
}
