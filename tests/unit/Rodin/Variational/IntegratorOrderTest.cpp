#include <gtest/gtest.h>

#include <set>

#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"
#include "Rodin/Assembly/Default.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  namespace
  {
    /// @brief Unit square meshed with @p n x @p n cells of the given geometry.
    LocalMesh unitSquare(Polytope::Type g, size_t n)
    {
      LocalMesh mesh = LocalMesh::UniformGrid(g, {n, n});
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      return mesh;
    }

    /**
     * @brief Unit square as one triangle pair plus one quadrilateral, so that a
     * single mesh carries two polytope geometries.
     *
     *   (0,1) 2---3 (1,1)      (1,1) 3---5 (2,1)
     *         | \ |                  |   |
     *   (0,0) 0---1 (1,0)      (1,0) 1---4 (2,0)
     */
    LocalMesh mixedGeometryMesh()
    {
      LocalMesh mesh = Mesh<Rodin::Context::Local>::Builder()
                         .initialize(2)
                         .nodes(6)
                         .vertex({0, 0})
                         .vertex({1, 0})
                         .vertex({0, 1})
                         .vertex({1, 1})
                         .vertex({2, 0})
                         .vertex({2, 1})
                         .polytope(Polytope::Type::Triangle, {0, 1, 2})
                         .polytope(Polytope::Type::Triangle, {1, 3, 2})
                         .polytope(Polytope::Type::Quadrilateral, {1, 4, 5, 3})
                         .finalize();
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      return mesh;
    }

    /**
     * @brief Integral of @p f over the mesh, assembled through a load form.
     *
     * The nodal basis is a partition of unity, so summing the assembled vector
     * gives @f$\sum_i\int f\varphi_i=\int f@f$. That makes the exactness of the
     * quadrature rule directly observable.
     */
    template <class Integrand>
    Real assembledIntegral(Integrand&& integrand, const auto& v)
    {
      LinearForm lf(v);
      lf = std::forward<Integrand>(integrand);
      lf.assemble();
      return lf.getVector().sum();
    }
  }

  // ==========================================================================
  // Order resolution: nullopt / constant / rule.
  // ==========================================================================

  /// @brief By default an integrator infers its order, reporting none of its own.
  TEST(Rodin_Variational_IntegratorOrder, Default_InfersOrder)
  {
    auto mesh = unitSquare(Polytope::Type::Triangle, 2);
    H1 fes(std::integral_constant<size_t, 1>{}, mesh);
    TestFunction v(fes);

    auto integ = Integral(RealFunction(1.0), v);
    const auto it = mesh.getPolytope(2, 0);
    EXPECT_FALSE(integ.getOrder(*it).has_value());
  }

  /// @brief A constant order is reported for every polytope.
  TEST(Rodin_Variational_IntegratorOrder, Constant_AppliesEverywhere)
  {
    auto mesh = unitSquare(Polytope::Type::Triangle, 2);
    H1 fes(std::integral_constant<size_t, 1>{}, mesh);
    TestFunction v(fes);

    auto integ = Integral(RealFunction(1.0), v);
    integ.setOrder(5);

    for (auto it = mesh.getCell(); it; ++it)
      EXPECT_EQ(integ.getOrder(*it).value_or(0), 5u);
  }

  /// @brief A rule is evaluated per polytope, and sees each one.
  TEST(Rodin_Variational_IntegratorOrder, Rule_EvaluatedPerPolytope)
  {
    auto mesh = unitSquare(Polytope::Type::Triangle, 2);
    H1 fes(std::integral_constant<size_t, 1>{}, mesh);
    TestFunction v(fes);

    std::set<Index> seen;
    auto integ = Integral(RealFunction(1.0), v);
    integ.setOrder([&seen](const Polytope& p) -> size_t {
      seen.insert(p.getIndex());
      return 2 + p.getIndex() % 3;
    });

    for (auto it = mesh.getCell(); it; ++it)
      EXPECT_EQ(integ.getOrder(*it).value_or(0), 2 + it->getIndex() % 3);
    EXPECT_EQ(seen.size(), mesh.getCellCount());
  }

  /// @brief Passing nullopt restores inference after an explicit order was set.
  TEST(Rodin_Variational_IntegratorOrder, Nullopt_RestoresInference)
  {
    auto mesh = unitSquare(Polytope::Type::Triangle, 2);
    H1 fes(std::integral_constant<size_t, 1>{}, mesh);
    TestFunction v(fes);

    auto integ = Integral(RealFunction(1.0), v);
    const auto it = mesh.getPolytope(2, 0);

    integ.setOrder(4);
    EXPECT_TRUE(integ.getOrder(*it).has_value());
    integ.setOrder(std::nullopt);
    EXPECT_FALSE(integ.getOrder(*it).has_value());
  }

  /// @brief The order survives the copy every form assembly makes of its integrator.
  TEST(Rodin_Variational_IntegratorOrder, Order_SurvivesCopy)
  {
    auto mesh = unitSquare(Polytope::Type::Triangle, 2);
    H1 fes(std::integral_constant<size_t, 1>{}, mesh);
    TestFunction v(fes);

    auto integ = Integral(RealFunction(1.0), v);
    integ.setOrder(6);

    auto copied = integ;
    const auto it = mesh.getPolytope(2, 0);
    EXPECT_EQ(copied.getOrder(*it).value_or(0), 6u);

    std::unique_ptr<Integrator> cloned(integ.copy());
    EXPECT_EQ(cloned->getOrder(*it).value_or(0), 6u);
  }

  // ==========================================================================
  // The order must reach the quadrature rule. A rule of degree >= deg(f)
  // integrates f exactly; a deficient one does not. This is the property whose
  // silent loss degrades a solve without any diagnostic.
  // ==========================================================================

  /// @brief A sufficient order integrates a quartic exactly; a deficient one does not.
  TEST(Rodin_Variational_IntegratorOrder, Order_ReachesQuadrature)
  {
    auto mesh = unitSquare(Polytope::Type::Triangle, 2);
    H1 fes(std::integral_constant<size_t, 1>{}, mesh);
    TestFunction v(fes);

      // f(x, y) = x^4, whose integral over the unit square is 1/5.
    auto f = RealFunction([](const Geometry::Point& p) -> Real {
      const Real x = p.getPhysicalCoordinates()(0);
      return x * x * x * x;
    });
    constexpr Real exact = Real(1) / Real(5);

    auto accurate = Integral(f, v);
    accurate.setOrder(8);
    EXPECT_NEAR(assembledIntegral(accurate, v), exact, 1e-12);

    auto deficient = Integral(f, v);
    deficient.setOrder(1);
    const Real coarse = assembledIntegral(deficient, v);
    EXPECT_GT(std::abs(coarse - exact), 1e-6)
      << "an order-1 rule must not integrate a quartic exactly; if it does, "
         "setOrder is not reaching the quadrature formula";
  }

  /// @brief A per-polytope rule reaches the quadrature just as a constant does.
  TEST(Rodin_Variational_IntegratorOrder, Rule_ReachesQuadrature)
  {
    auto mesh = unitSquare(Polytope::Type::Triangle, 2);
    H1 fes(std::integral_constant<size_t, 1>{}, mesh);
    TestFunction v(fes);

    auto f = RealFunction([](const Geometry::Point& p) -> Real {
      const Real x = p.getPhysicalCoordinates()(0);
      return x * x * x * x;
    });
    constexpr Real exact = Real(1) / Real(5);

    auto integ = Integral(f, v);
    integ.setOrder([](const Polytope&) -> size_t { return 8; });
    EXPECT_NEAR(assembledIntegral(integ, v), exact, 1e-12);
  }

  // ==========================================================================
  // Uniformity across integrator kinds and mesh configurations.
  // ==========================================================================

  /// @brief setOrder is honoured by bilinear integrators too.
  TEST(Rodin_Variational_IntegratorOrder, BilinearIntegrator_HonoursOrder)
  {
    auto mesh = unitSquare(Polytope::Type::Triangle, 2);
    H1 fes(std::integral_constant<size_t, 1>{}, mesh);
    TrialFunction u(fes);
    TestFunction v(fes);

    auto integ = Integral(Dot(u, v));
    const auto it = mesh.getPolytope(2, 0);
    EXPECT_FALSE(integ.getOrder(*it).has_value());

    integ.setOrder(7);
    EXPECT_EQ(integ.getOrder(*it).value_or(0), 7u);

      // The mass matrix of a partition of unity sums to the domain measure,
      // and an over-integrated rule must not change that.
    BilinearForm bf(u, v);
    bf = integ;
    bf.assemble();
    EXPECT_NEAR(bf.getOperator().sum(), 1.0, 1e-12);
  }

  /// @brief setOrder is honoured by boundary integrators too.
  TEST(Rodin_Variational_IntegratorOrder, BoundaryIntegrator_HonoursOrder)
  {
    auto mesh = unitSquare(Polytope::Type::Triangle, 2);
    H1 fes(std::integral_constant<size_t, 1>{}, mesh);
    TestFunction v(fes);

    auto integ = BoundaryIntegral(RealFunction(1.0), v);
    const auto it = mesh.getPolytope(1, 0);
    EXPECT_FALSE(integ.getOrder(*it).has_value());

    integ.setOrder(3);
    EXPECT_EQ(integ.getOrder(*it).value_or(0), 3u);
  }

  /// @brief Inference differs with the polynomial degree of the space.
  TEST(Rodin_Variational_IntegratorOrder, MixedOrder_InferenceFollowsDegree)
  {
    auto mesh = unitSquare(Polytope::Type::Triangle, 2);
    H1 linear(std::integral_constant<size_t, 1>{}, mesh);
    H1 quadratic(std::integral_constant<size_t, 2>{}, mesh);

    TrialFunction u1(linear);
    TestFunction v1(linear);
    TrialFunction u2(quadratic);
    TestFunction v2(quadratic);

      // Both spaces integrate their own mass matrix exactly, whatever order
      // each one infers.
    BilinearForm bf1(u1, v1);
    bf1 = Integral(Dot(u1, v1));
    bf1.assemble();
    EXPECT_NEAR(bf1.getOperator().sum(), 1.0, 1e-12);

    BilinearForm bf2(u2, v2);
    bf2 = Integral(Dot(u2, v2));
    bf2.assemble();
    EXPECT_NEAR(bf2.getOperator().sum(), 1.0, 1e-12);
  }

  /// @brief An explicit order overrides inference on a higher-degree space.
  TEST(Rodin_Variational_IntegratorOrder, MixedOrder_ExplicitOverridesInference)
  {
    auto mesh = unitSquare(Polytope::Type::Triangle, 2);
    H1 quadratic(std::integral_constant<size_t, 2>{}, mesh);
    TestFunction v(quadratic);

    auto integ = Integral(RealFunction(1.0), v);
    const auto it = mesh.getPolytope(2, 0);
    EXPECT_FALSE(integ.getOrder(*it).has_value());

    integ.setOrder(2);
    EXPECT_EQ(integ.getOrder(*it).value_or(0), 2u);
  }

  /// @brief On a quadrilateral mesh the order machinery behaves identically.
  TEST(Rodin_Variational_IntegratorOrder, Quadrilateral_HonoursOrder)
  {
    auto mesh = unitSquare(Polytope::Type::Quadrilateral, 2);
    H1 fes(std::integral_constant<size_t, 1>{}, mesh);
    TestFunction v(fes);

    auto f = RealFunction([](const Geometry::Point& p) -> Real {
      const Real x = p.getPhysicalCoordinates()(0);
      return x * x * x * x;
    });

    auto integ = Integral(f, v);
    integ.setOrder(8);
    EXPECT_NEAR(assembledIntegral(integ, v), Real(1) / Real(5), 1e-12);
  }

  /**
   * @brief On a mesh carrying two geometries a rule sees both, and may give
   * each its own order.
   *
   * This is the case a constant order cannot express: the inferred order is a
   * function of the polytope, so only a rule can follow it.
   */
  TEST(Rodin_Variational_IntegratorOrder, MixedPolytopes_RuleSeesEachGeometry)
  {
    auto mesh = mixedGeometryMesh();
    H1 fes(std::integral_constant<size_t, 1>{}, mesh);
    TestFunction v(fes);

    std::set<Polytope::Type> geometries;
    auto integ = Integral(RealFunction(1.0), v);
    integ.setOrder([&geometries](const Polytope& p) -> size_t {
      geometries.insert(p.getGeometry());
      return p.getGeometry() == Polytope::Type::Quadrilateral ? 4 : 2;
    });

    for (auto it = mesh.getCell(); it; ++it)
    {
      const size_t expected =
        it->getGeometry() == Polytope::Type::Quadrilateral ? 4u : 2u;
      EXPECT_EQ(integ.getOrder(*it).value_or(0), expected);
    }

    EXPECT_EQ(geometries.size(), 2u)
      << "the rule must observe both polytope geometries in the mesh";
  }

  /// @brief A per-geometry rule integrates correctly across a mixed mesh.
  TEST(Rodin_Variational_IntegratorOrder, MixedPolytopes_IntegratesExactly)
  {
    auto mesh = mixedGeometryMesh();
    H1 fes(std::integral_constant<size_t, 1>{}, mesh);
    TestFunction v(fes);

      // The mesh covers [0,2] x [0,1]; the integral of x^4 over it is 32/5.
    auto f = RealFunction([](const Geometry::Point& p) -> Real {
      const Real x = p.getPhysicalCoordinates()(0);
      return x * x * x * x;
    });

    auto integ = Integral(f, v);
    integ.setOrder([](const Polytope& p) -> size_t {
      return p.getGeometry() == Polytope::Type::Quadrilateral ? 8 : 6;
    });
    EXPECT_NEAR(assembledIntegral(integ, v), Real(32) / Real(5), 1e-10);
  }
}
