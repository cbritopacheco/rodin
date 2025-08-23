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
  TEST(Rodin_Variational_Tan, RealFunction_Zero)
  {
    RealFunction f(0.0);
    auto tan_result = Tan(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    EXPECT_NEAR(tan_result.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Tan, RealFunction_PiOverFour)
  {
    RealFunction f(M_PI / 4.0);
    auto tan_result = Tan(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);

    EXPECT_NEAR(tan_result.getValue(p), 1.0, 1e-10);
  }

  TEST(Rodin_Variational_Tan, RealFunction_Pi)
  {
    RealFunction f(M_PI);
    auto tan_result = Tan(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.2}};
    Point p(polytope, rc);

    EXPECT_NEAR(tan_result.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Tan, RealFunction_MinusPiOverFour)
  {
    RealFunction f(-M_PI / 4.0);
    auto tan_result = Tan(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.9}};
    Point p(polytope, rc);

    EXPECT_NEAR(tan_result.getValue(p), -1.0, 1e-10);
  }

  TEST(Rodin_Variational_Tan, RealFunction_PiOverSix)
  {
    RealFunction f(M_PI / 6.0);
    auto tan_result = Tan(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);

    EXPECT_NEAR(tan_result.getValue(p), 1.0 / std::sqrt(3.0), 1e-10);
  }

  TEST(Rodin_Variational_Tan, RealFunction_PiOverThree)
  {
    RealFunction f(M_PI / 3.0);
    auto tan_result = Tan(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.8}};
    Point p(polytope, rc);

    EXPECT_NEAR(tan_result.getValue(p), std::sqrt(3.0), 1e-10);
  }

  TEST(Rodin_Variational_Tan, GridFunction_Tan)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    P1 fes(mesh);
    GridFunction gf(fes);

    // Project a constant function
    RealFunction f(M_PI / 4.0);
    gf.project(f);

    auto tan_result = Tan(gf);

    // Test at a point
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.25, 0.25}};
    Point p(polytope, rc);

    EXPECT_NEAR(tan_result.getValue(p), 1.0, 1e-10);
  }

  TEST(Rodin_Variational_Tan, ChainedOperations)
  {
    RealFunction f(M_PI / 6.0);
    auto tan_f = tan(f);          // tan(π/6) = 1/√3
    auto tan_tan_f = tan(tan_f);  // tan(1/√3)

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.4, 0.6}};
    Point p(polytope, rc);

    EXPECT_NEAR(tan_tan_f.getValue(p), std::tan(1.0 / std::sqrt(3.0)), 1e-10);
  }

  TEST(Rodin_Variational_Tan, Periodicity)
  {
    RealFunction f1(M_PI / 6.0);
    RealFunction f2(M_PI / 6.0 + M_PI);
    auto tan1 = Tan(f1);
    auto tan2 = Tan(f2);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{1.0, 0.0}};
    Point p(polytope, rc);

    // tan(x) = tan(x + π)
    EXPECT_NEAR(tan1.getValue(p), tan2.getValue(p), 1e-10);
  }

  TEST(Rodin_Variational_Tan, SmallAngle)
  {
    RealFunction f(0.1);  // Small angle in radians
    auto tan_result = Tan(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.8}};
    Point p(polytope, rc);

    // For small angles, tan(x) ≈ x
    EXPECT_NEAR(tan_result.getValue(p), std::tan(0.1), 1e-10);
    EXPECT_NEAR(tan_result.getValue(p), 0.1, 1e-2);  // Approximate equality
  }

  TEST(Rodin_Variational_Tan, TanOfSum)
  {
    RealFunction f1(M_PI / 6.0);  // π/6
    RealFunction f2(M_PI / 12.0); // π/12
    auto sum = f1 + f2;           // π/6 + π/12 = π/4
    auto tan_sum = Tan(sum);  // tan(π/4) = 1

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    EXPECT_NEAR(tan_sum.getValue(p), 1.0, 1e-10);
  }

  TEST(Rodin_Variational_Tan, TanOddFunction)
  {
    RealFunction f(M_PI / 6.0);
    RealFunction neg_f(-M_PI / 6.0);
    auto tan_f = Tan(f);
    auto tan_neg_f = Tan(neg_f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 1.0}};
    Point p(polytope, rc);

    // tan(-x) = -tan(x) (odd function)
    EXPECT_NEAR(tan_neg_f.getValue(p), -tan_f.getValue(p), 1e-10);
  }

  TEST(Rodin_Variational_Tan, TanSinCosIdentity)
  {
    RealFunction f(M_PI / 6.0);
    auto sin_f = Sin(f);
    auto cos_f = Cos(f);
    auto sin_over_cos = sin_f / cos_f;  // sin/cos = tan
    auto tan_f = Tan(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.7, 0.3}};
    Point p(polytope, rc);

    // tan(x) = sin(x) / cos(x)
    EXPECT_NEAR(tan_f.getValue(p), sin_over_cos.getValue(p), 1e-10);
  }

  TEST(Rodin_Variational_Tan, TanWithMultiplication)
  {
    RealFunction f(M_PI / 8.0);
    RealFunction scalar(2.0);
    auto product = scalar * f;  // 2 * π/8 = π/4
    auto tan_product = Tan(product);  // tan(π/4) = 1

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);

    EXPECT_NEAR(tan_product.getValue(p), 1.0, 1e-10);
  }

  TEST(Rodin_Variational_Tan, TanLargeAngle)
  {
    RealFunction f(3.0 * M_PI);
    auto tan_result = Tan(f);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.1}};
    Point p(polytope, rc);

    // tan(3π) = 0 (multiple of π)
    EXPECT_NEAR(tan_result.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Tan, TanPositiveNegative)
  {
    RealFunction f1(M_PI / 6.0);   // First quadrant: positive
    RealFunction f2(5.0 * M_PI / 6.0); // Second quadrant: negative
    auto tan1 = Tan(f1);
    auto tan2 = Tan(f2);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);

    EXPECT_GT(tan1.getValue(p), 0.0);  // tan(π/6) > 0
    EXPECT_LT(tan2.getValue(p), 0.0);  // tan(5π/6) < 0
  }

  TEST(Rodin_Variational_Tan, TanSecantIdentity)
  {
    RealFunction f(M_PI / 6.0);
    auto tan_f = Tan(f);
    auto tan_squared = tan_f * tan_f;

    auto cos_f = Cos(f);
    auto sec_squared = RealFunction(1.0) / (cos_f * cos_f);  // sec²(x) = 1/cos²(x)
    auto identity = sec_squared - tan_squared;  // sec²(x) - tan²(x) = 1

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.9, 0.1}};
    Point p(polytope, rc);

    // sec²(x) - tan²(x) = 1
    EXPECT_NEAR(identity.getValue(p), 1.0, 1e-10);
  }

  TEST(Rodin_Variational_Tan, TanWithAbs)
  {
    RealFunction f(-M_PI / 4.0);
    auto abs_f = Abs(f);        // |−π/4| = π/4
    auto tan_abs = Tan(abs_f);   // tan(π/4) = 1

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);

    EXPECT_NEAR(tan_abs.getValue(p), 1.0, 1e-10);
  }

  TEST(Rodin_Variational_Tan, TanSpecialAngles)
  {
    // Test multiple special angles
    std::vector<std::pair<Real, Real>> angles_and_values = {
      {0.0, 0.0},
      {M_PI / 6.0, 1.0 / std::sqrt(3.0)},
      {M_PI / 4.0, 1.0},
      {M_PI / 3.0, std::sqrt(3.0)},
      {M_PI, 0.0}
    };

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.8}};
    Point p(polytope, rc);

    for (const auto& [angle, expected] : angles_and_values) {
      RealFunction f(angle);
      auto tan_result = Tan(f);
      EXPECT_NEAR(tan_result.getValue(p), expected, 1e-10) 
        << "Failed for angle " << angle;
    }
  }

  TEST(Rodin_Variational_Tan, TanContinuity)
  {
    // Test values approaching but not at π/2 where Tan has asymptotes
    RealFunction f1(M_PI / 2.0 - 0.1);
    RealFunction f2(M_PI / 2.0 - 0.01);
    auto tan1 = Tan(f1);
    auto tan2 = Tan(f2);

    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.2}};
    Point p(polytope, rc);

    // tan should be large and positive as we approach π/2 from the left
    EXPECT_GT(tan1.getValue(p), 9.0);    // tan(π/2 - 0.1) ≈ 9.97
    EXPECT_GT(tan2.getValue(p), 99.0);   // tan(π/2 - 0.01) ≈ 99.997
    EXPECT_GT(tan2.getValue(p), tan1.getValue(p));  // Should be increasing
  }
}
