/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <algorithm>
#include <array>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>
#include "Rodin/Alert/Exception.h"
#include "Rodin/Variational/H1/H1Element.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"

namespace Rodin::Tests::Unit
{
  using GradientGeometry = std::pair<Geometry::Polytope::Type, const char*>;

  /**
   * @brief Certifies ordered, reentrant bulk reference-gradient evaluation.
   *
   * Every entry is compared exactly with the independent public single-basis
   * derivative. Samples include interiors and reference vertices, preserving
   * historical singular-vertex conventions rather than asserting new traces.
   * Real/complex degrees one through six exercise distinct scalar
   * promotion and modal sizes. Independent caller-owned outputs and rejected invalid extents must
   * leave prior results and later invocations unchanged. No numerical
   * tolerance or coordinate-based identity is used.
   */
  class H1GradientEvaluationTest : public ::testing::TestWithParam<GradientGeometry>
  {
    protected:
      template <size_t K, class Scalar>
      static void check(Geometry::Polytope::Type geometry)
      {
        using namespace Variational;
        SCOPED_TRACE(K);
        SCOPED_TRACE((std::is_same_v<Scalar, Complex> ? "Complex" : "Real"));
        const H1Element<K, Scalar> element(geometry);
        const Geometry::Polytope::Traits traits(geometry);
        const auto& formula = QF::PolytopeQuadratureFormula::get(4, geometry);
        std::vector<Math::SpatialPoint> points;
        for (size_t qp = 0; qp < std::min(size_t(3), formula.getSize()); ++qp)
          points.push_back(formula.getPoint(qp));
        for (size_t vertex = 0; vertex < traits.getVertexCount(); ++vertex)
          points.push_back(traits.getVertex(vertex));
        ASSERT_FALSE(points.empty());
        const auto checkModalOutputs = [&]<class Modal, size_t Count>() {
          using ModalValues = typename Modal::Values;
          using ModalGradients = typename Modal::Gradients;
          typename Modal::Tabulation modalValues, modalGradients, both;
          for (auto& row : modalValues)
            row.fill(77);
          modalGradients = modalValues;
          both = modalValues;
          const Math::SpatialPoint invalid{0.25, 0.25};
          EXPECT_THROW(
            Modal::tabulate(ModalValues{modalValues}, invalid), Alert::Exception);
          EXPECT_THROW(
            Modal::tabulate(ModalGradients{modalGradients}, invalid), Alert::Exception);
          EXPECT_THROW(Modal::tabulate(both, invalid), Alert::Exception);
          EXPECT_EQ(modalValues, both);
          EXPECT_EQ(modalGradients, both);
          const Math::SpatialPoint point{0.125, 0.125, 0.25};
          Modal::tabulate(ModalValues{modalValues}, point);
          Modal::tabulate(ModalGradients{modalGradients}, point);
          Modal::tabulate(both, point);
          for (size_t mode = 0; mode < Count; ++mode)
          {
            EXPECT_EQ(modalValues[mode][0], both[mode][0]);
            EXPECT_EQ(modalGradients[mode][0], 77);
            for (size_t direction = 1; direction < 4; ++direction)
            {
              EXPECT_EQ(modalGradients[mode][direction], both[mode][direction]);
              EXPECT_EQ(modalValues[mode][direction], 77);
            }
          }
          Modal::tabulate(ModalGradients{modalValues}, point);
          Modal::tabulate(ModalValues{modalGradients}, point);
          EXPECT_EQ(modalValues, both);
          EXPECT_EQ(modalGradients, both);
        };
        if (geometry == Geometry::Polytope::Type::Tetrahedron)
        {
          checkModalOutputs
            .template operator()<DubinerTetrahedron<K>, FeketeTetrahedron<K>::Count>();
        }
        else if (geometry == Geometry::Polytope::Type::Pyramid)
        {
          checkModalOutputs
            .template operator()<PyramidModal<K>, PyramidIndex<K>::Count>();
        }
        using Values = typename H1Element<K, Scalar>::Values;
        using Gradients = typename H1Element<K, Scalar>::Gradients;
        const size_t count = element.getCount();
        const size_t dimension = traits.getDimension();
        const auto checkEvaluation = [&]<Geometry::Polytope::Type ReferenceGeometry>() {
          using Evaluation =
            typename H1Element<K, Scalar>::template Evaluation<ReferenceGeometry>;
          const Math::SpatialPoint invalid{0.25, 0.25};
          EXPECT_THROW(
            element.template tabulate<ReferenceGeometry>(invalid), Alert::Exception);
          EXPECT_THROW((Evaluation(invalid)), Alert::Exception);
          const H1Element<K, Scalar> incompatible(
            ReferenceGeometry == Geometry::Polytope::Type::Tetrahedron
              ? Geometry::Polytope::Type::Pyramid
              : Geometry::Polytope::Type::Tetrahedron);
          EXPECT_THROW(incompatible.template tabulate<ReferenceGeometry>(points.front()),
            Alert::Exception);
          for (const auto& point : points)
          {
            auto coordinates = point;
            const auto evaluation =
              element.template tabulate<ReferenceGeometry>(coordinates);
            const auto copy = evaluation;
            coordinates(0) += Real(0.03125);
            const auto other = element.template tabulate<ReferenceGeometry>(coordinates);
            ASSERT_EQ(evaluation.getCount(), count);
            ASSERT_EQ(other.getCount(), count);
            EXPECT_THROW(evaluation.getGradient(count), Alert::Exception);
            for (size_t local = 0; local < count; ++local)
            {
              const auto expected = element.getBasis(local).getGradient()(point);
              const auto actual = evaluation.getGradient(local);
              const auto copied = copy.getGradient(local);
              for (size_t direction = 0; direction < 3; ++direction)
              {
                EXPECT_EQ(actual(direction), expected(direction));
                EXPECT_EQ(copied(direction), expected(direction));
              }
            }
          }
          const auto independent = element.template tabulate<ReferenceGeometry>(
            Math::SpatialPoint{0.125, 0.125, 0.25});
          const auto expected =
            element.getBasis(0).getGradient()(Math::SpatialPoint{0.125, 0.125, 0.25});
          const auto actual = independent.getGradient(0);
          for (size_t direction = 0; direction < 3; ++direction)
            EXPECT_EQ(actual(direction), expected(direction));
        };
        if (geometry == Geometry::Polytope::Type::Tetrahedron)
          checkEvaluation.template operator()<Geometry::Polytope::Type::Tetrahedron>();
        else if (geometry == Geometry::Polytope::Type::Pyramid)
          checkEvaluation.template operator()<Geometry::Polytope::Type::Pyramid>();
        std::vector<Scalar> values(count, Scalar(77));
        std::vector<Scalar> gradients(count * dimension, Scalar(77));
        if (dimension > 0)
        {
          Math::SpatialPoint invalid(dimension - 1);
          invalid.setZero();
          try
          {
            element.tabulate(Gradients{gradients, dimension}, invalid);
            FAIL() << "Missing reference coordinates were accepted.";
          }
          catch (const Alert::Exception& error)
          {
            EXPECT_NE(std::string(error.what()).find("tabulate"), std::string::npos);
            EXPECT_NE(
              std::string(error.what()).find("In member function"), std::string::npos);
          }
          EXPECT_TRUE(std::all_of(gradients.begin(), gradients.end(),
            [](Scalar x) { return x == Scalar(77); }));
        }
        EXPECT_THROW(
          element.tabulate(Values{std::span(values).first(count - 1)}, points.front()),
          Alert::Exception);
        EXPECT_THROW(
          element.tabulate(Gradients{gradients, dimension + 1}, points.front()),
          Alert::Exception);
        Math::SpatialPoint padded(points.front());
        padded.resize(3);
        for (size_t direction = dimension; direction < 3; ++direction)
          padded(direction) = Real(7 + direction);
        points.push_back(padded);
        for (const auto& point : points)
        {
          element.tabulate(Values{values}, Gradients{gradients, dimension}, point);
          const auto originalValues = values;
          const auto originalGradients = gradients;
          std::vector<Scalar> other(count * dimension);
          element.tabulate(Gradients{other, dimension}, points.back());
          EXPECT_EQ(values, originalValues);
          EXPECT_EQ(gradients, originalGradients);
          for (size_t local = 0; local < count; ++local)
          {
            const auto& basis = element.getBasis(local);
            EXPECT_EQ(values[local], basis(point));
            const auto expected = basis.getGradient()(point);
            for (size_t direction = 0; direction < dimension; ++direction)
              EXPECT_EQ(gradients[local * dimension + direction], expected(direction));
          }
          element.tabulate(Values{values}, point);
          EXPECT_EQ(values, originalValues);
          element.tabulate(Gradients{gradients, dimension}, point);
          EXPECT_EQ(gradients, originalGradients);
        }
        typename H1Element<K, Scalar>::Tabulation table;
        element.tabulate(table, formula);
        ASSERT_EQ(table.getSize(), formula.getSize());
        ASSERT_EQ(table.getCount(), count);
        ASSERT_EQ(table.getDimension(), dimension);
        for (size_t q = 0; q < formula.getSize(); ++q)
        {
          element.tabulate(
            Values{values}, Gradients{gradients, dimension}, formula.getPoint(q));
          for (size_t local = 0; local < count; ++local)
          {
            EXPECT_EQ(table.getBasis(q, local), values[local]);
            for (size_t direction = 0; direction < dimension; ++direction)
            {
              EXPECT_EQ(table.template getDerivative<1>(q, local, direction),
                gradients[local * dimension + direction]);
            }
          }
        }
      }
  };

  TEST_P(H1GradientEvaluationTest, TypedOutputsMatchPublicBasisAndRemainIndependent)
  {
    const auto geometry = GetParam().first;
    check<1, Real>(geometry);
    check<2, Real>(geometry);
    check<3, Real>(geometry);
    check<4, Real>(geometry);
    check<5, Real>(geometry);
    check<6, Real>(geometry);
    check<1, Complex>(geometry);
    check<2, Complex>(geometry);
    check<3, Complex>(geometry);
    check<4, Complex>(geometry);
    check<5, Complex>(geometry);
    check<6, Complex>(geometry);
  }

  using G = Geometry::Polytope::Type;
  INSTANTIATE_TEST_SUITE_P(AllGeometries, H1GradientEvaluationTest,
    ::testing::Values(GradientGeometry{G::Point, "Point"},
      GradientGeometry{G::Segment, "Segment"}, GradientGeometry{G::Triangle, "Triangle"},
      GradientGeometry{G::Quadrilateral, "Quadrilateral"},
      GradientGeometry{G::Tetrahedron, "Tetrahedron"},
      GradientGeometry{G::Pyramid, "Pyramid"},
      GradientGeometry{G::Hexahedron, "Hexahedron"}, GradientGeometry{G::Wedge, "Wedge"}),
    [](const auto& info) { return info.param.second; });
}
