/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Physical evaluation reuse is confined to explicit cell bindings.
 */
#include <gtest/gtest.h>
#include "Rodin/Variational.h"
#include "../../../convergence/Convergence.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  /// @brief Mutable affine map without any external invalidation token.
  class BindingMap final : public PolytopeTransformation
  {
    public:
      /// @param dimension Reference and physical dimension.
      explicit BindingMap(size_t dimension)
        : PolytopeTransformation(dimension, dimension),
          m_scale(1),
          m_jacobians(0)
      {}

      /// @param scale New geometric dilation, published through fresh bindings.
      void setScale(Real scale)
      {
        m_scale = scale;
      }

      /// @returns Number of Jacobian evaluations performed by this map.
      size_t getJacobianCount() const
      {
        return m_jacobians;
      }

      /// @returns The affine order.
      Optional<size_t> getDegree() const override
      {
        return 1;
      }

      /**
       * @brief Applies the geometric dilation.
       * @param pc Physical coordinates.
       * @param rc Reference coordinates.
       */
      void transform(Math::SpatialPoint& pc, const Math::SpatialPoint& rc) const override
      {
        pc = m_scale * rc;
      }

      /**
       * @brief Evaluates the constant Jacobian.
       * @param matrix Physical Jacobian.
       * @param rc Reference coordinates, unused for an affine map.
       */
      void jacobian(Math::SpatialMatrix<Real>& matrix,
        [[maybe_unused]] const Math::SpatialPoint& rc) const override
      {
        ++m_jacobians;
        matrix = Math::SpatialMatrix<Real>::Identity(
          getPhysicalDimension(), getReferenceDimension());
        matrix *= m_scale;
      }

      /// @returns An independently owned copy of the map.
      BindingMap* copy() const noexcept override
      {
        return new BindingMap(*this);
      }

    private:
      Real m_scale;
      mutable size_t m_jacobians;
  };

  template <class Scalar>
  class BindingTest : public ::testing::Test
  {};
  using BindingScalars = ::testing::Types<Real, Complex>;
  TYPED_TEST_SUITE(BindingTest, BindingScalars);

  /// @brief Nested shape nodes prepare affine derivatives only once per cell binding.
  TYPED_TEST(BindingTest, NestedExpressionsReuseAffineDerivatives)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const auto cell = *mesh.getCell();
    auto map = std::make_unique<BindingMap>(2);
    auto& geometryMap = *map;
    mesh.setPolytopeTransformation({2, cell.getIndex()}, map.release());
    P1<TypeParam, LocalMesh> space(mesh);
    TestFunction v(space);
    const auto check = [&](auto expression) {
      expression.setPolytope(cell);
      const auto initial = geometryMap.getJacobianCount();
      for (const Real coordinate : {Real(0.125), Real(0.25), Real(0.5)})
      {
        const Point point(cell, Math::SpatialPoint{coordinate, Real(0.125)});
        const IntegrationPoint ip(point);
        expression.setIntegrationPoint(ip);
        EXPECT_GT(expression.getBasis(0).norm(), 0);
        EXPECT_EQ(geometryMap.getJacobianCount(), initial + 1);
      }
      expression.setPolytope(cell);
      const Point point(cell, Math::SpatialPoint{Real(0.25), Real(0.125)});
      const IntegrationPoint ip(point);
      expression.setIntegrationPoint(ip);
      EXPECT_EQ(geometryMap.getJacobianCount(), initial + 2);
    };
    check(TypeParam(2) * (Grad(v) + Grad(v)));
  }

  /// @brief Repeated cell bindings and unbound point calls observe in-place map changes.
  TYPED_TEST(BindingTest, P1DerivativesRefreshAllGeometries)
  {
    using Scalar = TypeParam;
    using G = Polytope::Type;
    for (const G geometry : {G::Segment, G::Triangle, G::Quadrilateral, G::Tetrahedron,
           G::Pyramid, G::Hexahedron, G::Wedge})
    {
      SCOPED_TRACE(static_cast<int>(geometry));
      auto mesh = Convergence::UniformGrid(geometry).makeMesh(2);
      const auto cell = *mesh.getCell();
      auto map = std::make_unique<BindingMap>(cell.getDimension());
      auto& geometryMap = *map;
      mesh.setPolytopeTransformation(
        {cell.getDimension(), cell.getIndex()}, map.release());
      P1<Scalar, LocalMesh> scalar(mesh);
      P1<Math::SpatialVector<Scalar>, LocalMesh> vector(mesh, cell.getDimension());
      TestFunction v(scalar);
      TestFunction w(vector);
      const auto& formula = QF::QuadratureFormula::get(4, geometry);
      const auto check = [&](auto make) {
        auto expression = make();
        using Value = std::remove_cvref_t<decltype(expression.getBasis(size_t(0)))>;
        const auto entries = [](const Value& value) {
          std::vector<Scalar> result;
          if constexpr (requires {
                          value.rows();
                          value.cols();
                        })
          {
            for (Eigen::Index i = 0; i < value.rows(); ++i)
            {
              for (Eigen::Index j = 0; j < value.cols(); ++j)
                result.push_back(value(i, j));
            }
          }
          else if constexpr (requires { value.size(); })
          {
            for (Eigen::Index i = 0; i < value.size(); ++i)
              result.push_back(value(i));
          }
          else
            result.push_back(value);
          return result;
        };
        geometryMap.setScale(1);
        const PolytopeQuadrature first(cell, formula);
        expression.setPolytope(cell);
        const IntegrationPoint firstPoint(first, 0);
        expression.setIntegrationPoint(firstPoint);
        std::vector<std::vector<Scalar>> original;
        for (size_t local = 0; local < expression.getDOFs(cell); ++local)
          original.push_back(entries(expression.getBasis(local)));
        geometryMap.setScale(2);
        const PolytopeQuadrature second(cell, formula);
        expression.setPolytope(cell);
        const IntegrationPoint bound(second, 0);
        expression.setIntegrationPoint(bound);
        for (size_t local = 0; local < original.size(); ++local)
        {
          const auto actual = entries(expression.getBasis(local));
          ASSERT_EQ(actual.size(), original[local].size());
          for (size_t k = 0; k < actual.size(); ++k)
            EXPECT_EQ(Scalar(2) * actual[k], original[local][k]);
        }
        auto direct = make();
        direct.setIntegrationPoint(bound);
        geometryMap.setScale(4);
        const PolytopeQuadrature third(cell, formula);
        const IntegrationPoint rebound(third, 0);
        direct.setIntegrationPoint(rebound);
        auto copied = expression;
        copied.setIntegrationPoint(rebound);
        for (size_t local = 0; local < original.size(); ++local)
        {
          const auto actual = entries(direct.getBasis(local));
          const auto cloned = entries(copied.getBasis(local));
          for (size_t k = 0; k < actual.size(); ++k)
          {
            EXPECT_EQ(Scalar(4) * actual[k], original[local][k]);
            EXPECT_EQ(actual[k], cloned[k]);
          }
        }
      };
      check([&] { return Grad(v); });
      check([&] { return Div(w); });
      check([&] { return Jacobian(w); });
    }
  }

  /// @brief Cell binding reaches nested weighted derivatives before replacing quadrature.
  TYPED_TEST(BindingTest, CompositeQuadratureRebindsSameCell)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const auto cell = *mesh.getCell();
    auto map = std::make_unique<BindingMap>(2);
    auto& geometryMap = *map;
    mesh.setPolytopeTransformation({2, cell.getIndex()}, map.release());
    P1<TypeParam, LocalMesh> space(mesh);
    TrialFunction u(space);
    TestFunction v(space);
    auto expression = TypeParam(2) * (Grad(u) + Grad(u));
    auto rule = Integral(Dot(expression, Grad(v)));
    rule.setPolytope(cell);
    const auto before = rule.integrate(0, 0);
    ASSERT_GT(std::abs(before), 0);
    geometryMap.setScale(2);
    rule.setPolytope(cell);
    auto fresh = Integral(Dot(TypeParam(2) * (Grad(u) + Grad(u)), Grad(v)));
    fresh.setPolytope(cell);
    EXPECT_EQ(rule.integrate(0, 0), fresh.integrate(0, 0));
    EXPECT_EQ(rule.integrate(0, 0), before);
  }

  /// @brief Unary, binary, matrix and component nodes forward fresh cell bindings.
  TYPED_TEST(BindingTest, NestedExpressionsPropagateBindings)
  {
    using Scalar = TypeParam;
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const auto cell = *mesh.getCell();
    auto map = std::make_unique<BindingMap>(2);
    auto& geometryMap = *map;
    mesh.setPolytopeTransformation({2, cell.getIndex()}, map.release());
    P1<Scalar, LocalMesh> scalar(mesh);
    P1<Math::SpatialVector<Scalar>, LocalMesh> vector(mesh, 2);
    TestFunction v(scalar);
    TestFunction w(vector);
    const auto& formula = QF::QuadratureFormula::get(4, cell.getGeometry());
    const auto check = [&](auto expression) {
      using Value = std::remove_cvref_t<decltype(expression.getBasis(size_t(0)))>;
      geometryMap.setScale(1);
      const PolytopeQuadrature first(cell, formula);
      const IntegrationPoint firstPoint(first, 0);
      expression.setPolytope(cell).setIntegrationPoint(firstPoint);
      std::vector<Value> original;
      for (size_t local = 0; local < expression.getDOFs(cell); ++local)
        original.push_back(expression.getBasis(local));
      geometryMap.setScale(2);
      const PolytopeQuadrature second(cell, formula);
      const IntegrationPoint secondPoint(second, 0);
      expression.setPolytope(cell).setIntegrationPoint(secondPoint);
      for (size_t local = 0; local < original.size(); ++local)
      {
        const Value actual = expression.getBasis(local);
        if constexpr (requires { actual.norm(); })
          EXPECT_EQ((Scalar(2) * actual - original[local]).norm(), 0);
        else
          EXPECT_EQ(Scalar(2) * actual, original[local]);
      }
    };
    check(-((Grad(v) + Grad(v)) / RealFunction(2)));
    check(Transpose(Jacobian(w)));
    check(Conjugate(Component(Transpose(Jacobian(w)), 0, 0)));
    check(Trace(Jacobian(w)));
  }

  /// @brief A copied trial reference owns its binding without mutating the source trial.
  TYPED_TEST(BindingTest, TrialReferencesOwnTheirEvaluationBinding)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const auto cell = *mesh.getCell();
    P1<TypeParam, LocalMesh> space(mesh);
    TrialFunction u(space);
    const auto& formula = QF::QuadratureFormula::get(4, cell.getGeometry());
    const PolytopeQuadrature quadrature(cell, formula);
    const IntegrationPoint first(quadrature, 0);
    const IntegrationPoint second(quadrature, 1);
    u.setPolytope(cell).setIntegrationPoint(first);
    const auto original = u.getBasis(0);
    std::unique_ptr<std::remove_pointer_t<decltype(u.copy())>> copy(u.copy());
    copy->setPolytope(cell).setIntegrationPoint(second);
    EXPECT_EQ(&u.getIntegrationPoint(), &first);
    EXPECT_EQ(&copy->getIntegrationPoint(), &second);
    EXPECT_EQ(u.getBasis(0), original);
    EXPECT_EQ(copy->getBasis(0),
      space.getFiniteElement(2, cell.getIndex())
        .getBasis(0)(second.getPoint().getReferenceCoordinates()));
  }

  /// @brief Failed point binding never publishes the preceding coefficient snapshot.
  TEST(BindingCache, RefreshAndFailureInvalidateOldValue)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const auto cell = *mesh.getCell();
    const auto& formula = QF::QuadratureFormula::get(4, cell.getGeometry());
    const PolytopeQuadrature quadrature(cell, formula);
    const IntegrationPoint point(quadrature, 0);
    Real value = 2;
    bool reject = false;
    auto function = RealFunction([&](const Point&) {
      if (reject)
        Alert::Exception() << "Rejected coefficient binding" << Alert::Raise;
      return value;
    });
    typename decltype(function)::Cache snapshot(function);
    snapshot.setPolytope(cell);
    EXPECT_EQ(snapshot.get(), nullptr);
    snapshot.setIntegrationPoint(point);
    ASSERT_NE(snapshot.get(), nullptr);
    EXPECT_EQ(*snapshot.get(), 2);
    value = 4;
    snapshot.setPolytope(cell);
    EXPECT_EQ(snapshot.get(), nullptr);
    snapshot.setIntegrationPoint(point);
    EXPECT_EQ(*snapshot.get(), 4);
    reject = true;
    EXPECT_THROW(snapshot.setIntegrationPoint(point), Alert::Exception);
    EXPECT_EQ(snapshot.get(), nullptr);
    reject = false;
    value = 8;
    snapshot.setIntegrationPoint(point);
    ASSERT_NE(snapshot.get(), nullptr);
    EXPECT_EQ(*snapshot.get(), 8);
  }

  /**
   * @brief Curved maps must not identify varying P1 gradients as constants.
   */
  TEST(Variational_Binding, CurvedPhysicalGradientIsNotConstant)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    const auto cell = *mesh.getCell();
    const RealH1Element<2> element(Polytope::Type::Triangle);
    PointCloud coefficients(2, element.getCount());
    for (size_t a = 0; a < element.getCount(); ++a)
    {
      const auto reference = element.getNode(a);
      Math::SpatialPoint physical;
      cell.getTransformation().transform(physical, reference);
      coefficients(0, a) = physical[0] + Real(0.25) * reference[0] * reference[0];
      coefficients(1, a) = physical[1];
    }
    mesh.setPolytopeTransformation({2, cell.getIndex()},
      new ParametricTransformation(std::move(coefficients), element));
    P1<Real> space(mesh);
    GridFunction field(space);
    field = RealFunction([](const Point& p) { return p.x(); });
    const auto gradient = Grad(field);
    const auto first = gradient(Point(cell, Math::SpatialPoint{0.2, 0.2}));
    const auto second = gradient(Point(cell, Math::SpatialPoint{0.6, 0.2}));
    EXPECT_GT((first - second).norm(), Real(1e-3));
    ASSERT_TRUE(gradient.getOrder(cell));
    EXPECT_GT(*gradient.getOrder(cell), 0);
    TestFunction test(space);
    const auto testOrder = Grad(test).getOrder(cell);
    ASSERT_TRUE(testOrder);
    EXPECT_GT(*testOrder, 0);
  }

}
