/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Reference-coordinate and sampling-provenance regressions.
 */
#include <gtest/gtest.h>
#include "Rodin/Variational.h"
#include "Rodin/Variational/P0g.h"
#include "Rodin/QF/GaussLegendre.h"
#include "../../../convergence/Convergence.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  namespace
  {
    /**
     * @brief Compares every scalar entry without a coordinate or value tolerance.
     * @param actual Evaluated scalar, vector or matrix.
     * @param expected Independently bound direct-point result.
     * @returns Whether dimensions and all entries agree exactly.
     */
    bool equal(const auto& actual, const auto& expected)
    {
      if constexpr (requires {
                      actual.rows();
                      actual.cols();
                    })
      {
        if (actual.rows() != expected.rows() || actual.cols() != expected.cols())
          return false;
        for (Eigen::Index row = 0; row < actual.rows(); ++row)
        {
          for (Eigen::Index col = 0; col < actual.cols(); ++col)
          {
            if (actual(row, col) != expected(row, col))
              return false;
          }
        }
        return true;
      }
      else if constexpr (requires {
                           actual.size();
                           actual(size_t(0));
                         })
      {
        if (actual.size() != expected.size())
          return false;
        for (size_t i = 0; i < actual.size(); ++i)
        {
          if (actual(i) != expected(i))
            return false;
        }
        return true;
      }
      else
        return actual == expected;
    }

    /**
     * @brief Checks two coordinates carrying the same sampling metadata.
     *
     * Sampling provenance remains formula/index zero while the supplied point
     * traverses distinct reference samples. Shape values and derivatives must
     * use that point, not substitute the origin formula's first coordinates.
     * Grid-function comparisons additionally retain coefficients across the
     * same transition. The direct oracle has no formula association.
     *
     * @tparam FieldDerivatives Whether nonconstant differential expressions are checked.
     * @tparam Directional Whether scalar H1 partial derivatives are checked.
     * @param space Supported scalar, vector or matrix finite-element space.
     */
    template <bool FieldDerivatives = true, bool Directional = false, class Space>
    void checkUnprovenSamples(Space& space)
    {
      using Scalar = typename FormLanguage::Traits<Space>::ScalarType;
      using Range = typename FormLanguage::Traits<Space>::RangeType;
      const auto cell = *space.getMesh().getCell();
      const auto& formula = QF::QuadratureFormula::get(4, cell.getGeometry());
      ASSERT_GT(formula.getSize(), 1);
      TestFunction u(space);
      const auto exercise = [&](auto make) {
        auto reused = make(u);
        for (const size_t sample : {size_t(0), size_t(1), size_t(0)})
        {
          const Point point(cell, formula.getPoint(sample));
          const IntegrationPoint associated(point, &formula, 0);
          const IntegrationPoint direct(point);
          reused.setIntegrationPoint(associated);
          auto fresh = make(u);
          fresh.setIntegrationPoint(direct);
          for (size_t a = 0; a < reused.getDOFs(cell); ++a)
          {
            const auto actual = reused.getBasis(a);
            const auto expected = fresh.getBasis(a);
            ASSERT_TRUE(equal(actual, expected)) << "sample=" << sample << " basis=" << a;
          }
        }
      };
      exercise([](auto& unknown) { return unknown + unknown; });
      GridFunction field(space);
      for (size_t a = 0; a < field.getSize(); ++a)
      {
        if constexpr (std::is_same_v<Scalar, Complex>)
          field.getData()(a) = Scalar(a + 1, 2 * a + 1);
        else
          field.getData()(a) = Scalar(a + 1);
      }
      const auto checkField = [&](const auto& expression) {
        for (const size_t sample : {size_t(0), size_t(1), size_t(0)})
        {
          const Point point(cell, formula.getPoint(sample));
          const IntegrationPoint associated(point, &formula, 0);
          const auto actual = expression.getValue(associated);
          const auto expected = expression.getValue(point);
          ASSERT_TRUE(equal(actual, expected)) << "field sample=" << sample;
        }
      };
      checkField(field);
      if constexpr (FieldDerivatives && FormLanguage::IsVectorRange<Range>::Value)
      {
        exercise([](auto& unknown) { return Jacobian(unknown); });
        exercise([](auto& unknown) { return Div(unknown); });
        if constexpr (FieldDerivatives)
        {
          checkField(Jacobian(field));
          checkField(Div(field));
        }
      }
      else if constexpr (FieldDerivatives && !FormLanguage::IsMatrixRange<Range>::Value)
      {
        exercise([](auto& unknown) { return Grad(unknown); });
        if constexpr (FieldDerivatives)
          checkField(Grad(field));
        if constexpr (Directional)
        {
          for (size_t axis = 0; axis < space.getMesh().getSpaceDimension(); ++axis)
            checkField(Derivative(axis, field));
        }
      }
    }

    class IntegrationPointTest : public ::testing::TestWithParam<Polytope::Type>
    {};
  }

  /// @brief Origin formula/index metadata cannot authorize reference-coordinate reuse.
  TEST_P(IntegrationPointTest, SamplingMetadataDoesNotReplaceReferenceCoordinates)
  {
    auto mesh = Convergence::UniformGrid(GetParam()).makeMesh(2);
    const auto fields = [&]<class Scalar>() {
      const auto ranges = [&]<class Range>() {
        const auto check = [&]<class FES>(FES space) { checkUnprovenSamples(space); };
        if constexpr (FormLanguage::IsMatrixRange<Range>::Value)
        {
          auto constant = P0<Range>(mesh, 2, 3);
          checkUnprovenSamples<false>(constant);
          auto global = P0g<Range, LocalMesh>(mesh, 2, 3);
          checkUnprovenSamples<false>(global);
          check(P1<Range>(mesh, 2, 3));
        }
        else if constexpr (FormLanguage::IsVectorRange<Range>::Value)
        {
          auto constant = P0<Range>(mesh, mesh.getDimension());
          checkUnprovenSamples<false>(constant);
          auto global = P0g<Range, LocalMesh>(mesh, mesh.getDimension());
          checkUnprovenSamples<false>(global);
          check(P1<Range>(mesh, mesh.getDimension()));
        }
        else
        {
          auto constant = P0<Range>(mesh);
          checkUnprovenSamples<false>(constant);
          auto global = P0g<Range, LocalMesh>(mesh);
          checkUnprovenSamples<false>(global);
          check(P1<Range>(mesh));
        }
        const auto degree = [&]<size_t K>(std::integral_constant<size_t, K> order) {
          SCOPED_TRACE(::testing::Message() << "degree=" << K);
          auto highOrder = [&]() {
            if constexpr (FormLanguage::IsMatrixRange<Range>::Value)
              return H1<K, Range>(order, mesh, 2, 3);
            else if constexpr (FormLanguage::IsVectorRange<Range>::Value)
              return H1<K, Range>(order, mesh, mesh.getDimension());
            else
              return H1<K, Range>(order, mesh);
          }();
          checkUnprovenSamples<true, true>(highOrder);
        };
        degree(std::integral_constant<size_t, 1>{});
        degree(std::integral_constant<size_t, 2>{});
        degree(std::integral_constant<size_t, 3>{});
        degree(std::integral_constant<size_t, 4>{});
        degree(std::integral_constant<size_t, 5>{});
        degree(std::integral_constant<size_t, 6>{});
      };
      ranges.template operator()<Scalar>();
      ranges.template operator()<Math::SpatialVector<Scalar>>();
      ranges.template operator()<Math::SpatialMatrix<Scalar>>();
    };
    fields.template operator()<Real>();
    fields.template operator()<Complex>();
  }

  /**
   * @brief Only owned samples retain reference association through value lifetimes.
   *
   * Origin metadata and coefficient-snapshot eligibility are preserved for
   * legacy bindings. Formula assignment invalidates reference association in
   * existing bindings and in new bindings to an old mapped quadrature. Moving
   * the quadrature preserves its owned sample storage and identity snapshot.
   */
  TEST_P(IntegrationPointTest, OwnedSamplesAndFormulaVersions)
  {
    auto mesh = Convergence::UniformGrid(GetParam()).makeMesh(2);
    const auto cell = *mesh.getCell();
    const auto& formula =
      QF::QuadratureFormula::get<QF::GaussLegendre>(cell.getGeometry(), size_t{2});
    PolytopeQuadrature points(cell, formula);
    const IntegrationPoint original(points, 0);
    IntegrationPoint copied(original);
    const IntegrationPoint moved(std::move(copied));
    EXPECT_EQ(original.getReferenceQuadratureFormula(), &formula);
    EXPECT_EQ(moved.getReferenceQuadratureFormula(), &formula);
    const IntegrationPoint legacy(original.getPoint(), &formula, 0);
    EXPECT_EQ(legacy.getQuadratureFormula(), &formula);
    EXPECT_EQ(legacy.getReferenceQuadratureFormula(), nullptr);
    EXPECT_EQ(legacy.getIndex(), 0);
    EXPECT_EQ(
      IntegrationPoint(original.getPoint()).getReferenceQuadratureFormula(), nullptr);
    P1<Real> space(mesh);
    TestFunction u(space);
    size_t calls = 0;
    RealFunction coefficient([&](const Point&) {
      ++calls;
      return Real(2);
    });
    auto expression = coefficient * (u + u);
    expression.setIntegrationPoint(legacy);
    for (size_t a = 0; a < expression.getDOFs(cell); ++a)
      expression.getBasis(a);
    EXPECT_EQ(calls, 1);
    auto transferred = std::move(points);
    EXPECT_EQ(IntegrationPoint(transferred, 0).getReferenceQuadratureFormula(), &formula);
    EXPECT_EQ(&original.getPoint(), &transferred.getPoint(0));
    const auto& replacement =
      QF::QuadratureFormula::get<QF::GaussLegendre>(cell.getGeometry(), size_t{4});
    EXPECT_EQ(original.getReferenceQuadratureFormula(), &formula);
    EXPECT_EQ(moved.getReferenceQuadratureFormula(), &formula);
    EXPECT_EQ(IntegrationPoint(transferred, 0).getReferenceQuadratureFormula(), &formula);
    const PolytopeQuadrature rebuilt(cell, replacement);
    EXPECT_EQ(IntegrationPoint(rebuilt, 0).getReferenceQuadratureFormula(), &replacement);
    EXPECT_THROW(IntegrationPoint(rebuilt, rebuilt.getSize()), Alert::Exception);
  }

  /**
   * @brief Facet derivatives and mapped averages/jumps preserve pointwise meaning.
   *
   * Boundary derivatives use the adjacent volume cell, rather than a
   * tangential derivative of the facet element. Interior averages and jumps
   * preserve origin metadata while their operands evaluate in cell charts.
   * Every result must equal the direct-point overload exactly. An interior
   * one-sided derivative without a selected trace must be rejected by both
   * overloads, rather than silently replaced by a tangential derivative.
   */
  TEST_P(IntegrationPointTest, FacetFunctionCoordinates)
  {
    auto mesh = Convergence::UniformGrid(GetParam()).makeMesh(3);
    const size_t dimension = mesh.getDimension();
    mesh.getConnectivity().compute(dimension - 1, dimension);
    const auto check = [&]<bool Directional, class Space>(Space& space) {
      using Scalar = typename FormLanguage::Traits<Space>::ScalarType;
      using Range = typename FormLanguage::Traits<Space>::RangeType;
      GridFunction field(space);
      for (size_t a = 0; a < field.getSize(); ++a)
      {
        if constexpr (std::is_same_v<Scalar, Complex>)
          field.getData()(a) = Scalar(a + 1, 2 * a + 1);
        else
          field.getData()(a) = Scalar(a + 1);
      }
      const auto compare = [&](const auto& expression, const IntegrationPoint& ip) {
        const auto actual = expression.getValue(ip);
        const auto expected = expression.getValue(ip.getPoint());
        ASSERT_TRUE(equal(actual, expected));
      };
      bool boundary = false, interior = false;
      for (auto face = mesh.getFace(); face; ++face)
      {
        const bool isBoundary = face->isBoundary();
        if (isBoundary ? boundary : interior)
          continue;
        const auto& formula = QF::QuadratureFormula::get(4, face->getGeometry());
        const PolytopeQuadrature points(*face, formula);
        const IntegrationPoint ip(points, 0);
        const auto exercise = [&](const auto& expression) {
          if (isBoundary)
            compare(expression, ip);
          else
          {
            EXPECT_THROW(expression.getValue(ip.getPoint()), Alert::Exception);
            EXPECT_THROW(expression.getValue(ip), Alert::Exception);
            compare(Average(expression), ip);
            compare(Jump(expression), ip);
          }
        };
        if constexpr (FormLanguage::IsVectorRange<Range>::Value)
        {
          exercise(Jacobian(field));
          exercise(Div(field));
        }
        else
        {
          exercise(Grad(field));
          if constexpr (Directional)
          {
            for (size_t axis = 0; axis < dimension; ++axis)
              exercise(Derivative(axis, field));
          }
        }
        if (isBoundary)
          boundary = true;
        else
          interior = true;
        if (boundary && interior)
          break;
      }
      EXPECT_TRUE(boundary);
      EXPECT_TRUE(interior);
    };
    const auto fields = [&]<class Scalar>() {
      P1<Scalar> scalar(mesh);
      P1<Math::SpatialVector<Scalar>> vector(mesh, dimension);
      H1<2, Scalar> quadratic(std::integral_constant<size_t, 2>{}, mesh);
      H1<2, Math::SpatialVector<Scalar>> quadraticVector(
        std::integral_constant<size_t, 2>{}, mesh, dimension);
      check.template operator()<false>(scalar);
      check.template operator()<false>(vector);
      check.template operator()<true>(quadratic);
      check.template operator()<false>(quadraticVector);
    };
    fields.template operator()<Real>();
    fields.template operator()<Complex>();
  }

  /// @brief A zero-dimensional owned sample has valid logical reference association.
  TEST(IntegrationPointReference, PointSample)
  {
    auto builder = LocalMesh::Builder();
    builder.initialize(3).nodes(1).vertex({0, 0, 0});
    const auto mesh = builder.finalize();
    const auto& formula = QF::QuadratureFormula::get(4, Polytope::Type::Point);
    const PolytopeQuadrature points(*mesh.getCell(), formula);
    const IntegrationPoint point(points, 0);
    EXPECT_EQ(point.getReferenceQuadratureFormula(), &formula);
    EXPECT_EQ(point.getPoint().getReferenceCoordinates().size(), 0);
  }

  INSTANTIATE_TEST_SUITE_P(AllGeometries, IntegrationPointTest,
    ::testing::Values(Polytope::Type::Segment, Polytope::Type::Triangle,
      Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron, Polytope::Type::Pyramid,
      Polytope::Type::Hexahedron, Polytope::Type::Wedge));
}
