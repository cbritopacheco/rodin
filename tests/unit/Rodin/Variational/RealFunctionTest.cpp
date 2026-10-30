#include <gtest/gtest.h>
#include <memory>
#include <type_traits>
#include "Rodin/Test/Random.h"

#include "Rodin/QF/Centroid.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::IO;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::Test::Random;

namespace Rodin::Tests::Unit
{
  namespace
  {
    struct IntegrationPointRealCallable
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
  }

  /// @brief Verifies constant real construction for variational real function by checking tolerance-based numerical results.
  TEST(Rodin_Variational_RealFunction, ConstantReal_Construction)
  {
    RealFunction f(3.14);
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);
    EXPECT_NEAR(f.getValue(p), 3.14, 1e-10);
  }

  /// @brief Verifies constant integer construction for variational real function by checking tolerance-based numerical results.
  TEST(Rodin_Variational_RealFunction, ConstantInteger_Construction)
  {
    RealFunction f(42);
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5, 0.5}};
    Point p(polytope, rc);
    EXPECT_NEAR(f.getValue(p), 42.0, 1e-10);
  }

  /// @brief Verifies copy constructor real for variational real function by checking tolerance-based numerical results, copy semantics.
  TEST(Rodin_Variational_RealFunction, CopyConstructor_Real)
  {
    RealFunction f(2.718);
    RealFunction f_copy(f);
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 0.0}};
    Point p(polytope, rc);
    EXPECT_NEAR(f_copy.getValue(p), f.getValue(p), 1e-10);
  }

  /// @brief Verifies move constructor real for variational real function by checking tolerance-based numerical results, move semantics.
  TEST(Rodin_Variational_RealFunction, MoveConstructor_Real)
  {
    RealFunction f(1.414);
    Real original_value = f.getValue(
      Point(*LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2}).getPolytope(2, 0),
        Math::Vector<Real>{{0.0, 0.0}}));
    RealFunction f_moved(std::move(f));
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 0.0}};
    Point p(polytope, rc);
    EXPECT_NEAR(f_moved.getValue(p), original_value, 1e-10);
  }

  /// @brief Verifies copy real for variational real function by checking tolerance-based numerical results, copy semantics.
  TEST(Rodin_Variational_RealFunction, Copy_Real)
  {
    RealFunction f(123.456);
    auto copied = f.copy();
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 0.0}};
    Point p(polytope, rc);
    EXPECT_NE(copied, nullptr);
    EXPECT_NEAR(copied->getValue(p), f.getValue(p), 1e-10);
    delete copied;
  }

  /// @brief Verifies copy integer for variational real function by checking tolerance-based numerical results, copy semantics.
  TEST(Rodin_Variational_RealFunction, Copy_Integer)
  {
    RealFunction f(789);
    auto copied = f.copy();
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.0, 0.0}};
    Point p(polytope, rc);
    EXPECT_NE(copied, nullptr);
    EXPECT_NEAR(copied->getValue(p), 789.0, 1e-10);
    delete copied;
  }

  /// @brief Verifies zero value for variational real function by checking tolerance-based numerical results.
  TEST(Rodin_Variational_RealFunction, ZeroValue)
  {
    RealFunction f(0.0);
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{1.0, 1.0}};
    Point p(polytope, rc);
    EXPECT_NEAR(f.getValue(p), 0.0, 1e-10);
  }

  /// @brief Verifies negative value for variational real function by checking tolerance-based numerical results.
  TEST(Rodin_Variational_RealFunction, NegativeValue)
  {
    RealFunction f(-99.99);
    // Create a simple point for testing
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.2, 0.8}};
    Point p(polytope, rc);
    EXPECT_NEAR(f.getValue(p), -99.99, 1e-10);
  }

  /// @brief Verifies constant function traced to boundary for variational real function by checking tolerance-based numerical results.
  TEST(Rodin_Variational_RealFunction, ConstantFunction_TracedToBoundary)
  {
    const Attribute interior_attr = 1;
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D - 1, D);
    for (auto it = mesh.getCell(); !it.end(); ++it)
      mesh.setAttribute(it.key(), interior_attr);
    RealFunction f(5.0);
    auto traced_f = f.traceOf(interior_attr);
    // Create a boundary point for testing
    auto it = mesh.getPolytope(D - 1, 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{0.5}};
    Point p(polytope, rc);
    EXPECT_NEAR(traced_f.getValue(p), 5.0, 1e-10);
  }

  /// @brief Verifies callable uses integration point fast path for variational real function by checking tolerance-based numerical results.
  TEST(Rodin_Variational_RealFunction, Callable_UsesIntegrationPointFastPath)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const Math::Vector<Real> rc{{0.25, 0.25}};
    Point p(*it, rc);
    const auto& qf = QF::QuadratureFormula::get<QF::Centroid>(Polytope::Type::Triangle);
    IntegrationPoint ip(p, &qf, 7);

    RealFunction f(IntegrationPointRealCallable{});

    EXPECT_NEAR(f.getValue(p), 1.0, 1e-10);
    EXPECT_NEAR(f.getValue(ip), 17.0, 1e-10);
  }

  /// @brief Verifies expression propagates integration point fast path for variational real function by checking tolerance-based numerical results.
  TEST(Rodin_Variational_RealFunction, Expression_PropagatesIntegrationPointFastPath)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const Math::Vector<Real> rc{{0.25, 0.25}};
    Point p(*it, rc);
    const auto& qf = QF::QuadratureFormula::get<QF::Centroid>(Polytope::Type::Triangle);
    IntegrationPoint ip(p, &qf, 8);

    auto f = abs(-RealFunction(IntegrationPointRealCallable{})) + 1.0;

    EXPECT_NEAR(f.getValue(p), 2.0, 1e-10);
    EXPECT_NEAR(f.getValue(ip), 19.0, 1e-10);
  }
}

namespace Rodin::Tests::Unit
{
  /// @brief Nested real wrappers own their operands and retain dispatch through the CRTP base.
  TEST(Rodin_Variational_RealFunction, NestedWrapperOwnsOperandAndPreservesDispatch)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const auto cell = *mesh.getCell();
    const auto& rule = QF::QuadratureFormula::get<QF::Centroid>(Polytope::Type::Triangle);
    const Point point(cell, rule.getPoint(0));
    const IntegrationPoint ip(point, &rule, 0);
    const auto clone = []<class Derived>(const FunctionBase<Derived>& operand) {
      using Wrapper = RealFunction<FunctionBase<Derived>>;
      static_assert(std::is_base_of_v<RealFunctionBase<Wrapper>, Wrapper>);
      Wrapper original(operand);
      Wrapper copied(original);
      Wrapper moved(std::move(copied));
      return std::unique_ptr<RealFunctionBase<Wrapper>>(moved.copy());
    };
    // Each temporary operand is destroyed before evaluating its independent copy.
    const auto constant = clone(RealFunction(Real(2.5)));
    EXPECT_EQ(constant->getValue(point), 2.5);
    EXPECT_EQ(constant->getValue(ip), 2.5);
    EXPECT_EQ((*constant)(point), 2.5);
    EXPECT_EQ(constant->getOrder(cell), Optional<size_t>(0));
    const auto callable = clone(RealFunction(IntegrationPointRealCallable{}));
    EXPECT_EQ(callable->getValue(point), 1);
    EXPECT_EQ(callable->getValue(ip), 10);
    EXPECT_FALSE(callable->getOrder(cell));

    RealFunction original(Real(3.5));
    const RealFunctionBase<decltype(original)>& base = original;
    const auto deduced = RealFunction(base);
    EXPECT_EQ(deduced.getValue(point), 3.5);
    EXPECT_EQ(deduced.getOrder(cell), Optional<size_t>(0));
  }

  /// @brief A wrapped P1 derivative preserves independent one-sided trace selections.
  TEST(Rodin_Variational_RealFunction, NestedWrapperCopiesOwnInterfaceTraceSelection)
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
    const Point point(*mesh.getVertex(1), Math::SpatialPoint{});
    const IntegrationPoint ip(point);
    const auto check = [&]<class Derived>(const FunctionBase<Derived>& operand) {
      using Wrapper = RealFunction<FunctionBase<Derived>>;
      Wrapper original(operand);
      EXPECT_EQ(&original.traceOf(Attribute(7)), &original);
      EXPECT_EQ(original.getValue(point), 2);
      std::unique_ptr<Wrapper> copied(original.copy());
      original.traceOf(FlatSet<Attribute>{Attribute(11)});
      EXPECT_EQ(original.getValue(ip), -4);
      EXPECT_EQ(copied->getValue(ip), 2);
      copied->traceOf(Attribute(11));
      EXPECT_EQ(copied->getValue(point), -4);
      EXPECT_EQ(copied->getOrder(*mesh.getCell()), Optional<size_t>(0));
    };
    check(Derivative(size_t(0), field));
    EXPECT_TRUE(field.getTraceDomain().empty());
  }
}
