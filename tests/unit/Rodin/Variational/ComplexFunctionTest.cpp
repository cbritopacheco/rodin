/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ComplexFunctionTest.cpp
 * @brief Tests for the ComplexFunction class.
 */
#include <gtest/gtest.h>
#include <cmath>
#include <complex>
#include <memory>

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
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{0.25, 0.25}};
  Point p(*it, rc);

  ComplexFunction f(5);

  Complex val = f.getValue(p);
  EXPECT_NEAR(val.real(), 5.0, 1e-10);
  EXPECT_NEAR(val.imag(), 0.0, 1e-10);
}

/// @brief Verifies real constant for variational complex function by checking tolerance-based numerical results.
TEST(Rodin_Variational_ComplexFunction, RealConstant)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{0.25, 0.25}};
  Point p(*it, rc);

  ComplexFunction f(3.14);

  Complex val = f.getValue(p);
  EXPECT_NEAR(val.real(), 3.14, 1e-10);
  EXPECT_NEAR(val.imag(), 0.0, 1e-10);
}

/// @brief Verifies complex constant for variational complex function by checking tolerance-based numerical results.
TEST(Rodin_Variational_ComplexFunction, ComplexConstant)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{0.25, 0.25}};
  Point p(*it, rc);

  ComplexFunction f(Complex(2.0, 3.0));

  Complex val = f.getValue(p);
  EXPECT_NEAR(val.real(), 2.0, 1e-10);
  EXPECT_NEAR(val.imag(), 3.0, 1e-10);
}

/// @brief Verifies from real and imag parts for variational complex function by checking tolerance-based numerical results.
TEST(Rodin_Variational_ComplexFunction, FromRealAndImagParts)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{0.25, 0.25}};
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
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{0.25, 0.25}};
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
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{0.25, 0.25}};
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
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{0.25, 0.25}};
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
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{0.25, 0.25}};
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
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  auto it = mesh.getPolytope(mesh.getDimension(), 0);

  ComplexFunction f(Complex(1.0, 0.0));
  auto order = f.getOrder(*it);
  EXPECT_TRUE(order.has_value());
  EXPECT_EQ(*order, 0u);
}

/// @brief Verifies get order callable for variational complex function by checking false predicates.
TEST(Rodin_Variational_ComplexFunction, GetOrderCallable)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  auto it = mesh.getPolytope(mesh.getDimension(), 0);

  auto lambda = [](const Geometry::Point&) -> Complex { return Complex(1.0, 0.0); };
  ComplexFunction f(lambda);
  auto order = f.getOrder(*it);
  EXPECT_FALSE(order.has_value());
}

/// @brief Verifies callable uses integration point fast path for variational complex function by checking tolerance-based numerical results.
TEST(Rodin_Variational_ComplexFunction, Callable_UsesIntegrationPointFastPath)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{0.25, 0.25}};
  Point p(*it, rc);
  const auto& qf = QF::QuadratureFormula::get<QF::Centroid>(Polytope::Type::Triangle);
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
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  auto it = mesh.getPolytope(mesh.getDimension(), 0);
  const Math::Vector<Real> rc{{0.25, 0.25}};
  Point p(*it, rc);
  const auto& qf = QF::QuadratureFormula::get<QF::Centroid>(Polytope::Type::Triangle);
  IntegrationPoint ip(p, &qf, 4);

  RealFunction re(IntegrationPointRealPartCallable{});
  RealFunction im(IntegrationPointImagPartCallable{});
  auto f = ComplexFunction(re, im);

  EXPECT_NEAR(f.getValue(p).real(), 1.0, 1e-10);
  EXPECT_NEAR(f.getValue(p).imag(), 2.0, 1e-10);
  EXPECT_NEAR(f.getValue(ip).real(), 14.0, 1e-10);
  EXPECT_NEAR(f.getValue(ip).imag(), 24.0, 1e-10);
}

/// @brief Polymorphic constant copies retain values and independent trace selections.
TEST(Rodin_Variational_ComplexFunction, PolymorphicConstantCopiesOwnTraceMetadata)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  const auto cell = *mesh.getCell();
  const auto& qf = QF::QuadratureFormula::get<QF::Centroid>(Polytope::Type::Triangle);
  const Point point(cell, qf.getPoint(0));
  const IntegrationPoint ip(point, &qf, 0);
  const FlatSet<Attribute> selected{Attribute(7), Attribute(11)};
  const auto check = [&](auto source, Complex expected) {
    using Function = decltype(source);
    source.traceOf(selected);
    ComplexFunctionBase<Function>& base = source;
    std::unique_ptr<ComplexFunctionBase<Function>> copied(base.copy());
    source.traceOf(Attribute(19));
    EXPECT_EQ(copied->getValue(point), expected);
    EXPECT_EQ(copied->getValue(ip), expected);
    EXPECT_EQ(copied->getOrder(cell), Optional<size_t>(0));
    EXPECT_EQ(copied->getTraceDomain(), selected);
    copied->traceOf(Attribute(23));
    EXPECT_TRUE(source.getTraceDomain().contains(Attribute(19)));
    EXPECT_FALSE(source.getTraceDomain().contains(Attribute(23)));
  };
  check(ComplexFunction(Integer(-3)), Complex(-3, 0));
  check(ComplexFunction(Real(2.5)), Complex(2.5, 0));
  check(ComplexFunction(Complex(-2, 3)), Complex(-2, 3));
}

/// @brief A copied pair of callables outlives its source and keeps integration-point dispatch.
TEST(Rodin_Variational_ComplexFunction, PolymorphicCallablePairOutlivesSource)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  const auto cell = *mesh.getCell();
  const auto& qf = QF::QuadratureFormula::get<QF::Centroid>(Polytope::Type::Triangle);
  const Point point(cell, qf.getPoint(0));
  const IntegrationPoint ip(point, &qf, 0);
  struct Part
  {
    Real m_offset;
    Real operator()(const Point&) const
    {
      return m_offset;
    }
    Real operator()(const IntegrationPoint& sample) const
    {
      return m_offset + 10 + static_cast<Real>(sample.getIndex());
    }
  };
  using Function = decltype(ComplexFunction(Part{3}, Part{-2}));
  std::unique_ptr<ComplexFunctionBase<Function>> copied;
  {
    Function source(Part{3}, Part{-2});
    source.traceOf(Attribute(31));
    copied.reset(source.copy());
  }
  ASSERT_TRUE(copied);
  EXPECT_EQ(copied->getValue(point), Complex(3, -2));
  EXPECT_EQ(copied->getValue(ip), Complex(13, 8));
  EXPECT_FALSE(copied->getOrder(cell));
  EXPECT_TRUE(copied->getTraceDomain().contains(Attribute(31)));
}

/// @brief Explicit nested wrappers own copies and preserve scalar and quadrature dispatch.
TEST(Rodin_Variational_ComplexFunction, NestedWrapperOwnsAndDispatchesThroughItsBase)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  const auto cell = *mesh.getCell();
  const auto& rule = QF::QuadratureFormula::get<QF::Centroid>(Polytope::Type::Triangle);
  const Point point(cell, rule.getPoint(0));
  const IntegrationPoint ip(point, &rule, 0);
  const auto check = [&]<class Derived>(const FunctionBase<Derived>& nested,
                       Complex expectedPoint, Complex expectedQuadrature,
                       Optional<size_t> expectedOrder) {
    using Wrapper = ComplexFunction<FunctionBase<Derived>>;
    static_assert(std::is_base_of_v<ComplexFunctionBase<Wrapper>, Wrapper>);
    std::unique_ptr<ComplexFunctionBase<Wrapper>> owned;
    {
      Wrapper original(nested);
      Wrapper copied(original);
      Wrapper moved(std::move(copied));
      owned.reset(moved.copy());
    }
    EXPECT_EQ(owned->getValue(point), expectedPoint);
    EXPECT_EQ(owned->getValue(ip), expectedQuadrature);
    EXPECT_EQ((*owned)(point), expectedPoint);
    EXPECT_EQ(owned->getOrder(cell), expectedOrder);
  };
  check(RealFunction(Real(2.5)), Complex(2.5, 0), Complex(2.5, 0), 0);
  check(ComplexFunction(Complex(-2, 3)), Complex(-2, 3), Complex(-2, 3), 0);
  check(ComplexFunction(IntegrationPointComplexCallable{}), Complex(1, 2),
    Complex(10, 20), {});
  EXPECT_EQ(ComplexFunction(Integer(-3)).getValue(), -3);
  EXPECT_EQ(ComplexFunction(Real(2.5)).getValue(), 2.5);
  EXPECT_EQ(ComplexFunction(Complex(-2, 3)).getValue(), Complex(-2, 3));
}

/// @brief Nested wrapper copies keep independent one-sided derivatives across a continuous kink.
TEST(Rodin_Variational_ComplexFunction, NestedWrapperCopiesOwnInterfaceTraceSelection)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
  mesh.getConnectivity().compute(0, 1);
  mesh.setAttribute({1, 0}, Attribute(7));
  mesh.setAttribute({1, 1}, Attribute(11));
  P1<Complex> space(mesh);
  GridFunction field(space);
  field = [](const Point& point) {
    return point.x() <= 1 ? Complex(2, 3) * point.x()
                          : Complex(2, 3) + Complex(-4, 5) * (point.x() - 1);
  };
  const Point point(*mesh.getVertex(1), Math::SpatialPoint{});
  const IntegrationPoint ip(point);
  const auto check = [&]<class Derived>(const FunctionBase<Derived>& operand) {
    using Wrapper = ComplexFunction<FunctionBase<Derived>>;
    Wrapper original(operand);
    EXPECT_EQ(&original.traceOf(Attribute(7)), &original);
    EXPECT_EQ(original.getValue(point), Complex(2, 3));
    std::unique_ptr<Wrapper> copied(original.copy());
    original.traceOf(FlatSet<Attribute>{Attribute(11)});
    EXPECT_EQ(original.getValue(ip), Complex(-4, 5));
    EXPECT_EQ(copied->getValue(ip), Complex(2, 3));
    copied->traceOf(Attribute(11));
    EXPECT_EQ(copied->getValue(point), Complex(-4, 5));
  };
  check(Derivative(size_t(0), field));
  EXPECT_TRUE(field.getTraceDomain().empty());
}

/// @brief Explicit function pairs forward traces, preserve polynomial order, and own copies.
TEST(Rodin_Variational_ComplexFunction, FunctionPairOwnsTraceSelectionAndReportsOrder)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
  mesh.getConnectivity().compute(0, 1);
  mesh.setAttribute({1, 0}, Attribute(7));
  mesh.setAttribute({1, 1}, Attribute(11));
  P1<Real> space(mesh);
  GridFunction field(space);
  field = [](const Point& point) {
    return point.x() <= 1 ? 2 * point.x() : 2 - 4 * (point.x() - 1);
  };
  const auto pair = []<class Re, class Im>(
                      const FunctionBase<Re>& re, const FunctionBase<Im>& im) {
    return ComplexFunction<FunctionBase<Re>, FunctionBase<Im>>(re, im);
  };
  const Point interface(*mesh.getVertex(1), Math::SpatialPoint{});
  const IntegrationPoint ip(interface);
  const auto constant = RealFunction(Real(3));
  auto original = pair(Derivative(size_t(0), field), constant);
  EXPECT_EQ(&original.traceOf(Attribute(7)), &original);
  EXPECT_EQ(original.getValue(interface), Complex(2, 3));
  using Pair = decltype(original);
  std::unique_ptr<ComplexFunctionBase<Pair>> owned;
  {
    Pair copied(original);
    Pair moved(std::move(copied));
    owned.reset(moved.copy());
  }
  original.traceOf(FlatSet<Attribute>{Attribute(11)});
  EXPECT_EQ(original.getValue(ip), Complex(-4, 3));
  EXPECT_EQ(owned->getValue(ip), Complex(2, 3));
  EXPECT_EQ(owned->getOrder(*mesh.getCell()), Optional<size_t>(0));
  EXPECT_EQ(pair(field, constant).getOrder(*mesh.getCell()), Optional<size_t>(1));
  const auto unknown = RealFunction([](const Point& p) { return p.x(); });
  EXPECT_FALSE(pair(unknown, constant).getOrder(*mesh.getCell()));
  EXPECT_FALSE(pair(constant, unknown).getOrder(*mesh.getCell()));
}
