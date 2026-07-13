/**
 * @file ComplexFunctionTest.cpp
 * @brief Tests for the ComplexFunction class.
 */
#include <gtest/gtest.h>
#include <cmath>
#include <complex>

#include "Rodin/QF/Centroid.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  struct IntegrationPointComplexCallable
  {
    Complex operator()(const Geometry::Point&) const
    {
      return {1.0, 2.0};
    }

    Complex operator()(const IntegrationPoint& ip) const
    {
      const auto idx = static_cast<Real>(ip.getIndex());
      return {10.0 + idx, 20.0 + idx};
    }
  };

  struct IntegrationPointRealPartCallable
  {
    Real operator()(const Geometry::Point&) const
    {
      return 1.0;
    }

    Real operator()(const IntegrationPoint& ip) const
    {
      return 10.0 + static_cast<Real>(ip.getIndex());
    }
  };

  struct IntegrationPointImagPartCallable
  {
    Real operator()(const Geometry::Point&) const
    {
      return 2.0;
    }

    Real operator()(const IntegrationPoint& ip) const
    {
      return 20.0 + static_cast<Real>(ip.getIndex());
    }
  };
}

/// @brief Verifies integer constant for variational complex function by checking tolerance-based numerical results.
TEST(Rodin_Variational_ComplexFunction, IntegerConstant)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  ComplexFunction f(5);

  Complex val = f.getValue(p);
  EXPECT_NEAR(val.real(), 5.0, 1e-10);
  EXPECT_NEAR(val.imag(), 0.0, 1e-10);
}

/// @brief Verifies real constant for variational complex function by checking tolerance-based numerical results.
TEST(Rodin_Variational_ComplexFunction, RealConstant)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  ComplexFunction f(3.14);

  Complex val = f.getValue(p);
  EXPECT_NEAR(val.real(), 3.14, 1e-10);
  EXPECT_NEAR(val.imag(), 0.0, 1e-10);
}

/// @brief Verifies complex constant for variational complex function by checking tolerance-based numerical results.
TEST(Rodin_Variational_ComplexFunction, ComplexConstant)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  ComplexFunction f(Complex(2.0, 3.0));

  Complex val = f.getValue(p);
  EXPECT_NEAR(val.real(), 2.0, 1e-10);
  EXPECT_NEAR(val.imag(), 3.0, 1e-10);
}

/// @brief Verifies from real and imag parts for variational complex function by checking tolerance-based numerical results.
TEST(Rodin_Variational_ComplexFunction, FromRealAndImagParts)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  RealFunction re(4.0);
  RealFunction im(5.0);
  auto f = ComplexFunction(re, im);

  Complex val = f.getValue(p);
  EXPECT_NEAR(val.real(), 4.0, 1e-10);
  EXPECT_NEAR(val.imag(), 5.0, 1e-10);
}

/// @brief Verifies callable for variational complex function by checking tolerance-based numerical results.
TEST(Rodin_Variational_ComplexFunction, Callable)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  auto lambda = [](const Geometry::Point&) -> Complex { return Complex(1.0, 2.0); };
  ComplexFunction f(lambda);

  Complex val = f.getValue(p);
  EXPECT_NEAR(val.real(), 1.0, 1e-10);
  EXPECT_NEAR(val.imag(), 2.0, 1e-10);
}

/// @brief Verifies callable pair for variational complex function by checking tolerance-based numerical results.
TEST(Rodin_Variational_ComplexFunction, CallablePair)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  auto reLambda = [](const Geometry::Point&) -> Real { return 10.0; };
  auto imLambda = [](const Geometry::Point&) -> Real { return 20.0; };
  ComplexFunction f(reLambda, imLambda);

  Complex val = f.getValue(p);
  EXPECT_NEAR(val.real(), 10.0, 1e-10);
  EXPECT_NEAR(val.imag(), 20.0, 1e-10);
}

/// @brief Verifies copy constant for variational complex function by checking tolerance-based numerical results, copy semantics.
TEST(Rodin_Variational_ComplexFunction, CopyConstant)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  ComplexFunction f(Complex(7.0, 8.0));
  auto copy = f;

  Complex val = copy.getValue(p);
  EXPECT_NEAR(val.real(), 7.0, 1e-10);
  EXPECT_NEAR(val.imag(), 8.0, 1e-10);
}

/// @brief Verifies copy composite for variational complex function by checking tolerance-based numerical results, copy semantics.
TEST(Rodin_Variational_ComplexFunction, CopyComposite)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);

  RealFunction re(1.0);
  RealFunction im(2.0);
  auto f = ComplexFunction(re, im);
  auto copy = f;

  Complex val = copy.getValue(p);
  EXPECT_NEAR(val.real(), 1.0, 1e-10);
  EXPECT_NEAR(val.imag(), 2.0, 1e-10);
}

/// @brief Verifies get order constant for variational complex function by checking exact expected values, true predicates.
TEST(Rodin_Variational_ComplexFunction, GetOrderConstant)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);

  ComplexFunction f(Complex(1.0, 0.0));
  auto order = f.getOrder(*it);
  EXPECT_TRUE(order.has_value());
  EXPECT_EQ(*order, 0u);
}

/// @brief Verifies get order callable for variational complex function by checking false predicates.
TEST(Rodin_Variational_ComplexFunction, GetOrderCallable)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);

  auto lambda = [](const Geometry::Point&) -> Complex { return Complex(1.0, 0.0); };
  ComplexFunction f(lambda);
  auto order = f.getOrder(*it);
  EXPECT_FALSE(order.has_value());
}

/// @brief Verifies callable uses integration point fast path for variational complex function by checking tolerance-based numerical results.
TEST(Rodin_Variational_ComplexFunction, Callable_UsesIntegrationPointFastPath)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);
  QF::Centroid qf(Polytope::Type::Triangle);
  IntegrationPoint ip(p, &qf, 3);

  ComplexFunction f(IntegrationPointComplexCallable{});

  EXPECT_NEAR(f.getValue(p).real(), 1.0, 1e-10);
  EXPECT_NEAR(f.getValue(p).imag(), 2.0, 1e-10);
  EXPECT_NEAR(f.getValue(ip).real(), 13.0, 1e-10);
  EXPECT_NEAR(f.getValue(ip).imag(), 23.0, 1e-10);
}

/// @brief Verifies real and imag parts use integration point fast path for variational complex function by checking tolerance-based numerical results.
TEST(Rodin_Variational_ComplexFunction, RealAndImagParts_UseIntegrationPointFastPath)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 2, 2 });
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{ 0.25, 0.25 }};
  Point p(*it, rc);
  QF::Centroid qf(Polytope::Type::Triangle);
  IntegrationPoint ip(p, &qf, 4);

  RealFunction re(IntegrationPointRealPartCallable{});
  RealFunction im(IntegrationPointImagPartCallable{});
  auto f = ComplexFunction(re, im);

  EXPECT_NEAR(f.getValue(p).real(), 1.0, 1e-10);
  EXPECT_NEAR(f.getValue(p).imag(), 2.0, 1e-10);
  EXPECT_NEAR(f.getValue(ip).real(), 14.0, 1e-10);
  EXPECT_NEAR(f.getValue(ip).imag(), 24.0, 1e-10);
}
