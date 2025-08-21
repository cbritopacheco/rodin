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
  TEST(Rodin_Variational_Cosine, RealFunction_Zero)
  {
    RealFunction f(0.0);
    auto cos_result = Cosine(f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(cos_result.getValue(p), 1.0, 1e-10);
  }

  TEST(Rodin_Variational_Cosine, RealFunction_PiOverTwo)
  {
    RealFunction f(M_PI / 2.0);
    auto cos_result = Cosine(f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(cos_result.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Cosine, RealFunction_Pi)
  {
    RealFunction f(M_PI);
    auto cos_result = Cosine(f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.2}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(cos_result.getValue(p), -1.0, 1e-10);
  }

  TEST(Rodin_Variational_Cosine, RealFunction_ThreePiOverTwo)
  {
    RealFunction f(3.0 * M_PI / 2.0);
    auto cos_result = Cosine(f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.9}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(cos_result.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Cosine, RealFunction_TwoPi)
  {
    RealFunction f(2.0 * M_PI);
    auto cos_result = Cosine(f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(cos_result.getValue(p), 1.0, 1e-10);
  }

  TEST(Rodin_Variational_Cosine, RealFunction_PiOverSix)
  {
    RealFunction f(M_PI / 6.0);
    auto cos_result = Cosine(f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.8}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(cos_result.getValue(p), std::sqrt(3.0) / 2.0, 1e-10);
  }

  TEST(Rodin_Variational_Cosine, RealFunction_PiOverFour)
  {
    RealFunction f(M_PI / 4.0);
    auto cos_result = Cosine(f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.7, 0.3}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(cos_result.getValue(p), std::sqrt(2.0) / 2.0, 1e-10);
  }

  TEST(Rodin_Variational_Cosine, RealFunction_PiOverThree)
  {
    RealFunction f(M_PI / 3.0);
    auto cos_result = Cosine(f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.9, 0.1}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(cos_result.getValue(p), 0.5, 1e-10);
  }

  TEST(Rodin_Variational_Cosine, RealFunction_NegativeValue)
  {
    RealFunction f(-M_PI / 3.0);
    auto cos_result = Cosine(f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 0.0}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(cos_result.getValue(p), 0.5, 1e-10);
  }

  TEST(Rodin_Variational_Cosine, GridFunction_CosineWave)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    P1 fes(mesh);
    GridFunction gf(fes);
    
    // Project a constant function
    RealFunction f(M_PI / 3.0);
    gf.project(f);
    
    auto cos_result = Cosine(gf);
    
    // Test at a point
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.25, 0.25}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(cos_result.getValue(p), 0.5, 1e-10);
  }

  TEST(Rodin_Variational_Cosine, ChainedOperations)
  {
    RealFunction f(M_PI / 3.0);
    auto cos_f = Cosine(f);          // cos(π/3) = 0.5
    auto cos_cos_f = Cosine(cos_f);  // cos(0.5) ≈ 0.877583
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.4, 0.6}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(cos_cos_f.getValue(p), std::cos(0.5), 1e-10);
  }

  TEST(Rodin_Variational_Cosine, Periodicity)
  {
    RealFunction f1(M_PI / 4.0);
    RealFunction f2(M_PI / 4.0 + 2.0 * M_PI);
    auto cos1 = Cosine(f1);
    auto cos2 = Cosine(f2);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{1.0, 0.0}};
    Point p(polytope, rc);
    
    // cos(x) = cos(x + 2π)
    EXPECT_NEAR(cos1.getValue(p), cos2.getValue(p), 1e-10);
  }

  TEST(Rodin_Variational_Cosine, CosineOfSum)
  {
    RealFunction f1(M_PI / 6.0);  // π/6
    RealFunction f2(M_PI / 3.0);  // π/3
    auto sum = f1 + f2;           // π/6 + π/3 = π/2
    auto cos_sum = Cosine(sum);   // cos(π/2) = 0
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(cos_sum.getValue(p), 0.0, 1e-10);
  }

  TEST(Rodin_Variational_Cosine, SmallAngle)
  {
    RealFunction f(0.1);  // Small angle in radians
    auto cos_result = Cosine(f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.8, 0.8}};
    Point p(polytope, rc);
    
    // For small angles, cos(x) ≈ 1 - x²/2
    EXPECT_NEAR(cos_result.getValue(p), std::cos(0.1), 1e-10);
    EXPECT_NEAR(cos_result.getValue(p), 1.0 - 0.01/2.0, 1e-4);  // Approximate equality
  }

  TEST(Rodin_Variational_Cosine, LargeAngle)
  {
    RealFunction f(10.0 * M_PI);
    auto cos_result = Cosine(f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.1, 0.1}};
    Point p(polytope, rc);
    
    // cos(10π) = 1 (even multiple of π)
    EXPECT_NEAR(cos_result.getValue(p), 1.0, 1e-10);
  }

  TEST(Rodin_Variational_Cosine, CosineSquaredIdentity)
  {
    RealFunction f(M_PI / 6.0);
    auto cos_f = Cosine(f);
    auto cos_squared = cos_f * cos_f;
    
    // cos²(π/6) = (√3/2)² = 3/4
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.3, 0.7}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(cos_squared.getValue(p), 3.0 / 4.0, 1e-10);
  }

  TEST(Rodin_Variational_Cosine, CosineWithMultiplication)
  {
    RealFunction f(M_PI / 6.0);
    RealFunction scalar(2.0);
    auto product = scalar * f;  // 2 * π/6 = π/3
    auto cos_product = Cosine(product);  // cos(π/3) = 1/2
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.6, 0.4}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(cos_product.getValue(p), 0.5, 1e-10);
  }

  TEST(Rodin_Variational_Cosine, CosineEvenFunction)
  {
    RealFunction f(M_PI / 4.0);
    RealFunction neg_f(-M_PI / 4.0);
    auto cos_f = Cosine(f);
    auto cos_neg_f = Cosine(neg_f);
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 1.0}};
    Point p(polytope, rc);
    
    // cos(-x) = cos(x) (even function)
    EXPECT_NEAR(cos_neg_f.getValue(p), cos_f.getValue(p), 1e-10);
  }

  TEST(Rodin_Variational_Cosine, SineCosineIdentity)
  {
    RealFunction f(M_PI / 4.0);
    auto sin_f = Sine(f);
    auto cos_f = Cosine(f);
    auto sin_squared = sin_f * sin_f;
    auto cos_squared = cos_f * cos_f;
    auto identity = sin_squared + cos_squared;  // Should equal 1
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.7, 0.7}};
    Point p(polytope, rc);
    
    // sin²(x) + cos²(x) = 1
    EXPECT_NEAR(identity.getValue(p), 1.0, 1e-10);
  }

  TEST(Rodin_Variational_Cosine, CosinePhaseShift)
  {
    RealFunction f(M_PI / 6.0);
    RealFunction phase_shift(M_PI / 2.0);
    auto shifted = f + phase_shift;  // π/6 + π/2 = 2π/3
    auto cos_shifted = Cosine(shifted);  // cos(2π/3) = -1/2
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.2}};
    Point p(polytope, rc);
    
    EXPECT_NEAR(cos_shifted.getValue(p), -0.5, 1e-10);
  }

  TEST(Rodin_Variational_Cosine, CosinePythagoreanIdentity)
  {
    RealFunction f(M_PI / 3.0);
    auto sin_f = Sine(f);     // sin(π/3) = √3/2
    auto cos_f = Cosine(f);   // cos(π/3) = 1/2
    
    auto sin_sq_plus_cos_sq = sin_f * sin_f + cos_f * cos_f;
    
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.9, 0.9}};
    Point p(polytope, rc);
    
    // Fundamental trigonometric identity
    EXPECT_NEAR(sin_sq_plus_cos_sq.getValue(p), 1.0, 1e-10);
  }
}