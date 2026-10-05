/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Per-point coefficient evaluation and reassembly regressions. */
#include <gtest/gtest.h>
#include "Rodin/Variational.h"
#ifdef RODIN_USE_OPENMP
#include "Rodin/Assembly/Sequential.h"
#include "Rodin/Assembly/OpenMP.h"
#include "Rodin/Assembly/Default.h"
#endif

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  class LazyVectorFunction final : public FunctionBase<LazyVectorFunction>
  {
    public:
      LazyVectorFunction(size_t& calls, Real& scale)
        : m_calls(calls),
          m_scale(scale),
          m_value{1, 2}
      {}
      auto getValue(const Point&) const
      {
        ++m_calls.get();
        return m_value * m_scale.get();
      }
      Optional<size_t> getOrder(const Polytope&) const
      {
        return 0;
      }
      LazyVectorFunction* copy() const noexcept override
      {
        return new LazyVectorFunction(*this);
      }

    private:
      std::reference_wrapper<size_t> m_calls;
      std::reference_wrapper<Real> m_scale;
      Eigen::Vector2d m_value;
  };

  template <class Expr>
  class CountingShape final : public ShapeFunctionBase<CountingShape<Expr>,
                                typename Expr::FESType, Expr::SpaceType>
  {
    public:
      using FESType = typename Expr::FESType;
      static constexpr ShapeFunctionSpaceType SpaceType = Expr::SpaceType;
      using Parent = ShapeFunctionBase<CountingShape, FESType, SpaceType>;
      CountingShape(Expr expr, size_t& calls)
        : Parent(expr.getFiniteElementSpace()),
          m_expr(std::move(expr)),
          m_calls(calls)
      {}
      CountingShape(const CountingShape& other)
        : Parent(other),
          m_expr(other.m_expr),
          m_calls(other.m_calls)
      {}
      const auto& getLeaf() const
      {
        return m_expr.getLeaf();
      }
      const auto& getFiniteElementSpace() const
      {
        return m_expr.getFiniteElementSpace();
      }
      auto getDOFs(const Polytope& cell) const
      {
        return m_expr.getDOFs(cell);
      }
      const IntegrationPoint& getIntegrationPoint() const
      {
        return m_expr.getIntegrationPoint();
      }
      CountingShape& setIntegrationPoint(const IntegrationPoint& ip)
      {
        m_expr.setIntegrationPoint(ip);
        return *this;
      }
      auto getBasis(size_t i) const
      {
        ++m_calls.get();
        return Variational::Internal::materializeProduct(m_expr.getBasis(i));
      }
      auto getOrder(const Polytope& cell) const
      {
        return m_expr.getOrder(cell);
      }
      CountingShape* copy() const noexcept override
      {
        return new CountingShape(*this);
      }

    private:
      Expr m_expr;
      std::reference_wrapper<size_t> m_calls;
  };

  template <class Expr>
  void checkExpression(Expr expr, const Polytope& cell, size_t& calls)
  {
    const auto& qf = QF::PolytopeQuadratureFormula::get(4, cell.getGeometry());
    const auto& q = cell.getQuadrature(qf);
    const auto count = expr.getDOFs(cell);
    calls = 0;
    for (size_t qp = 0; qp < q.getSize(); ++qp)
    {
      const IntegrationPoint ip(q.getPoint(qp), &qf, qp);
      expr.setIntegrationPoint(ip);
      for (size_t i = 0; i < count; ++i)
      {
        const auto a = Variational::Internal::materializeProduct(expr.getBasis(i));
        const auto b = Variational::Internal::materializeProduct(expr.getBasis(i));
        EXPECT_LE(std::abs(Math::dot(a, a) - Math::dot(b, b)), 1e-12);
      }
      EXPECT_EQ(calls, qp + 1);
    }
  }
}

TEST(CoefficientEvaluation, AllCoefficientShapeOperators)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);
  const auto cell = *mesh.getCell();
  H1<3, Real> fes(std::integral_constant<size_t, 3>{}, mesh);
  TestFunction v(fes);
  size_t calls = 0;
  auto f = RealFunction([&](const Point& p) {
    ++calls;
    return 2.0 + p.x();
  });
  checkExpression(f * v, cell, calls);
  checkExpression(v * f, cell, calls);
  checkExpression(Dot(f, v), cell, calls);
  checkExpression(Dot(v, f), cell, calls);
  checkExpression(v / f, cell, calls);
}

TEST(CoefficientEvaluation, RebindingCopiesAndDirectEvaluation)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);
  P1<Real> fes(mesh);
  TestFunction v(fes);
  size_t calls = 0;
  Real scale = 2;
  auto f = RealFunction([&](const Point&) {
    ++calls;
    return scale;
  });
  auto expr = f * v;
  const auto cell = *mesh.getCell();
  const auto& qf = QF::PolytopeQuadratureFormula::get(4, cell.getGeometry());
  const Point p(cell, qf.getPoint(0));
  const IntegrationPoint ip(p, &qf, 0);
  expr.setIntegrationPoint(ip);
  const auto original = expr.getBasis(0);
  EXPECT_EQ(calls, 1);
  scale = 6;
  expr.setIntegrationPoint(ip);
  EXPECT_NEAR(expr.getBasis(0), 3 * original, 1e-12);
  EXPECT_EQ(calls, 2);
  auto copy = expr;
  copy.setIntegrationPoint(ip);
  EXPECT_NEAR(copy.getBasis(0), 3 * original, 1e-12);
  EXPECT_EQ(calls, 3);
  auto moved = std::move(copy);
  moved.setIntegrationPoint(ip);
  EXPECT_NEAR(moved.getBasis(0), 3 * original, 1e-12);
  const IntegrationPoint direct(p);
  moved.setIntegrationPoint(direct);
  const auto a = moved.getBasis(0);
  scale = 12;
  EXPECT_NEAR(moved.getBasis(0), 2 * a, 1e-12);
}

TEST(CoefficientEvaluation, CacheLifecycleAndIndependentBindings)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  const auto cell = *mesh.getCell();
  const auto& qf = QF::PolytopeQuadratureFormula::get(4, cell.getGeometry());
  const Point point(cell, qf.getPoint(0));
  const IntegrationPoint ip(point, &qf, 0);
  Real scale = 2;
  size_t calls = 0;
  auto f = RealFunction([&](const Point&) {
    ++calls;
    return scale;
  });
  using Cache = decltype(f)::Cache;
  static_assert(Cache::Enabled);
  Cache first;
  Cache second;
  EXPECT_EQ(first.get(), nullptr);
  first.refresh(f, ip);
  ASSERT_NE(first.get(), nullptr);
  EXPECT_EQ(*first.get(), 2);
  EXPECT_EQ(calls, 1);
  scale = 5;
  second.refresh(f, ip);
  EXPECT_EQ(*first.get(), 2);
  ASSERT_NE(second.get(), nullptr);
  EXPECT_EQ(*second.get(), 5);
  first.refresh(f, ip);
  EXPECT_EQ(*first.get(), 5);
  EXPECT_EQ(calls, 3);
  Cache copied(first);
  EXPECT_EQ(copied.get(), nullptr);
  copied.refresh(f, ip);
  copied = first;
  EXPECT_EQ(copied.get(), nullptr);
  Cache moved(std::move(second));
  ASSERT_NE(moved.get(), nullptr);
  EXPECT_EQ(*moved.get(), 5);
  copied = std::move(moved);
  ASSERT_NE(copied.get(), nullptr);
  EXPECT_EQ(*copied.get(), 5);
  const size_t before = calls;
  first.refresh(f, IntegrationPoint(point));
  EXPECT_EQ(first.get(), nullptr);
  EXPECT_EQ(calls, before);
}

TEST(CoefficientEvaluation, LazyFunctionValuesRetainDirectEvaluation)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  P1<Math::SpatialVector<Real>> fes(mesh, 2);
  TestFunction v(fes);
  size_t calls = 0;
  Real scale = 2;
  LazyVectorFunction f(calls, scale);
  static_assert(!LazyVectorFunction::Cache::Enabled);
  const auto cell = *mesh.getCell();
  const auto& qf = QF::PolytopeQuadratureFormula::get(4, cell.getGeometry());
  const Point point(cell, qf.getPoint(0));
  const IntegrationPoint ip(point, &qf, 0);
  auto expr = Dot(f, v);
  expr.setIntegrationPoint(ip);
  EXPECT_EQ(calls, 0);
  const auto initial = expr.getBasis(0);
  EXPECT_EQ(calls, 1);
  scale = 5;
  EXPECT_NEAR(expr.getBasis(0), 2.5 * initial, 1e-12);
  EXPECT_EQ(calls, 2);
}

TEST(CoefficientEvaluation, VectorMatrixAndTensorCallables)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);
  const auto cell = *mesh.getCell();
  H1<2, Math::SpatialVector<Real>> vectors(std::integral_constant<size_t, 2>{}, mesh, 2);
  H1<2, Math::SpatialMatrix<Complex>> matrices(
    std::integral_constant<size_t, 2>{}, mesh, 2, 3);
  TestFunction v(vectors);
  TestFunction t(matrices);
  size_t calls = 0;
  auto vector = VectorFunction(size_t(2), [&](const Point&) {
    ++calls;
    return Math::SpatialVector<Real>{1, 2};
  });
  auto matrix = MatrixFunction(2, 2, [&](const Point&) {
    ++calls;
    return Math::SpatialMatrix<Real>::Identity(2, 2);
  });
  auto tensor = TensorFunction([&](const Point&) {
    ++calls;
    Math::SpatialTensor<Complex, 3> value(2, 3, 2);
    for (size_t i = 0; i < value.size(); ++i)
      value[i] = Complex(1 + i, 2);
    return value;
  });
  checkExpression(Dot(vector, v), cell, calls);
  checkExpression(Dot(v, vector), cell, calls);
  checkExpression(matrix * v, cell, calls);
  checkExpression(Dot(tensor, Jacobian(t)), cell, calls);
  checkExpression(Dot(Jacobian(t), tensor), cell, calls);
}

namespace
{
  template <class TrialFES, class TestFES>
  void checkGenericRule(const TrialFES& fes, const TestFES& testFes)
  {
    const auto& mesh = fes.getMesh();
    TrialFunction u(fes);
    TestFunction v(testFes);
    size_t trialCalls = 0, testCalls = 0, coefficientCalls = 0;
    Real scale = 2;
    auto f = RealFunction([&](const Point&) {
      ++coefficientCalls;
      return scale;
    });
    CountingShape trial(1.0 * u, trialCalls);
    CountingShape test(f * v, testCalls);
    auto rule = Integral(trial, test);
    static_assert(!decltype(rule)::Specialized);
    rule.setOrder(4);
    const auto cell = *mesh.getCell();
    const auto& qf = QF::PolytopeQuadratureFormula::get(4, cell.getGeometry());
    const auto n = fes.getFiniteElement(cell.getDimension(), cell.getIndex()).getCount();
    const auto nte =
      testFes.getFiniteElement(cell.getDimension(), cell.getIndex()).getCount();
    rule.setPolytope(cell);
    EXPECT_EQ(trialCalls, qf.getSize() * n);
    EXPECT_EQ(testCalls, qf.getSize() * nte);
    EXPECT_EQ(coefficientCalls, qf.getSize());
    auto reference = Integral(u, v);
    reference.setOrder(4);
    reference.setPolytope(cell);
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < nte; ++j)
        EXPECT_NEAR(rule.integrate(i, j), scale * reference.integrate(i, j), 1e-12);
    scale = 5;
    rule.setPolytope(cell);
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < nte; ++j)
        EXPECT_NEAR(rule.integrate(i, j), scale * reference.integrate(i, j), 1e-12);
    EXPECT_EQ(coefficientCalls, 2 * qf.getSize());
  }
}

TEST(CoefficientEvaluation, GenericRuleEvaluatesEachBasisOnceAndReassembles)
{
  for (auto geometry :
    {Polytope::Type::Point, Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Wedge, Polytope::Type::Hexahedron})
  {
    SCOPED_TRACE(static_cast<int>(geometry));
    const auto d = Polytope::Traits(geometry).getDimension();
    auto mesh = d == 0
      ? LocalMesh::Builder().initialize(1).nodes(1).vertex({0}).finalize()
      : d == 1 ? LocalMesh::UniformGrid(geometry, {2})
      : d == 2 ? LocalMesh::UniformGrid(geometry, {2, 2})
               : LocalMesh::UniformGrid(geometry, {2, 2, 2});
    for (size_t dim = 1; dim <= d; ++dim)
      for (size_t lower = 0; lower < dim; ++lower)
        mesh.getConnectivity().compute(dim, lower);
    H1<1, Real> p1(std::integral_constant<size_t, 1>{}, mesh);
    H1<2, Real> p2(std::integral_constant<size_t, 2>{}, mesh);
    H1<3, Real> p3(std::integral_constant<size_t, 3>{}, mesh);
    checkGenericRule(p1, p2);
    checkGenericRule(p2, p3);
  }
}

namespace
{
  template <class ScalarSpace, class VectorSpace>
  void checkSpecialized(const ScalarSpace& scalar, const VectorSpace& vector)
  {
    TrialFunction u(scalar);
    TestFunction v(scalar);
    TrialFunction w(vector);
    TestFunction z(vector);
    size_t calls = 0;
    auto f = RealFunction([&](const Point& p) {
      ++calls;
      return 2 + p.x();
    });
    auto b = VectorFunction(size_t(2), [&](const Point&) {
      ++calls;
      return Math::SpatialVector<Real>{1, 2};
    });
    auto A = MatrixFunction(2, 2, [&](const Point&) {
      ++calls;
      return Math::SpatialMatrix<Real>::Identity(2, 2);
    });
    const auto cell = *scalar.getMesh().getCell();
    const auto& qf = QF::PolytopeQuadratureFormula::get(4, cell.getGeometry());
    const auto check = [&](auto rule, size_t coefficients = 1) {
      rule.setOrder(4);
      calls = 0;
      rule.setPolytope(cell);
      EXPECT_EQ(calls, coefficients * qf.getSize());
    };
    check(Integral(f, v));
    check(Integral(f * u, v));
    check(Integral(f * Grad(u), Grad(v)));
    check(Integral(f * Jacobian(w), Jacobian(z)));
    check(Integral(Jacobian(w) * b, z));
    check(Integral(Dot(A, Jacobian(w)), Dot(A, Jacobian(z))));
  }
}

TEST(CoefficientEvaluation, SpecializedP1AndH1HoistCoefficients)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);
  checkSpecialized(P1<Real>(mesh), P1<Math::SpatialVector<Real>>(mesh, 2));
  checkSpecialized(H1<2, Real>(std::integral_constant<size_t, 2>{}, mesh),
    H1<2, Math::SpatialVector<Real>>(std::integral_constant<size_t, 2>{}, mesh, 2));
}

TEST(CoefficientEvaluation, P0P0gAndP1UseTheSameCoefficientBinding)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);
  const auto cell = *mesh.getCell();
  size_t calls = 0;
  auto f = RealFunction([&](const Point& p) {
    ++calls;
    return 2 + p.x();
  });
  const auto check = [&](const auto& fes) {
    TestFunction v(fes);
    checkExpression(f * v, cell, calls);
    checkExpression(v * f, cell, calls);
    checkExpression(Dot(f, v), cell, calls);
    checkExpression(Dot(v, f), cell, calls);
    checkExpression(v / f, cell, calls);
  };
  check(P0<Real>(mesh));
  check(P0g<Real, LocalMesh>(mesh));
  check(P1<Real>(mesh));
}

#ifdef RODIN_USE_OPENMP
TEST(CoefficientEvaluation, GenericParallelAssemblyMatchesSequential)
{
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);
  const auto check = [&](const auto& fes) {
    TrialFunction u(fes);
    TestFunction v(fes);
    auto f = RealFunction([](const Point& p) { return 2 + p.x(); });
    BilinearForm form(u, v);
    form = Integral(u, f * v);
    using Form = decltype(form);
    using Matrix = Math::SparseMatrix<Real>;
    const auto input = Assembly::BilinearFormAssemblyInput{
      fes, fes, form.getLocalIntegrators(), form.getGlobalIntegrators()};
    Matrix sequential;
    Assembly::Sequential<Matrix, Form>{}.execute(sequential, input);
    for (size_t pass = 0; pass < 3; ++pass)
    {
      Matrix parallel;
      Assembly::OpenMP<Matrix, Form>{}.execute(parallel, input);
      EXPECT_LE((parallel - sequential).norm(), 1e-12 * (1 + sequential.norm()));
    }
  };
  check(H1<2, Real>(std::integral_constant<size_t, 2>{}, mesh));
  check(H1<2, Math::SpatialVector<Real>>(std::integral_constant<size_t, 2>{}, mesh, 2));
  check(
    H1<2, Math::SpatialMatrix<Real>>(std::integral_constant<size_t, 2>{}, mesh, 2, 3));
}
#endif
