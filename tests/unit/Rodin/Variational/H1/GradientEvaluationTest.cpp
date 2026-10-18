/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <algorithm>
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
   * promotion and modal sizes. Nested evaluation and a throwing consumer must
   * leave subsequent entries and later invocations unchanged. No numerical
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
        if (traits.getDimension() > 0)
        {
          Math::SpatialPoint invalid(traits.getDimension() - 1);
          invalid.setZero();
          bool invoked = false;
          try
          {
            element.evaluateGradients(
              [&](size_t, const auto&) { invoked = true; }, invalid);
            FAIL() << "Missing reference coordinates were accepted.";
          }
          catch (const Alert::Exception& error)
          {
            EXPECT_NE(
              std::string(error.what()).find("evaluateGradients"), std::string::npos);
            EXPECT_NE(std::string(error.what())
                        .find("coordinate components, got " +
                          std::to_string(static_cast<size_t>(invalid.size())) + "."),
              std::string::npos);
          }
          EXPECT_FALSE(invoked);
        }
        Math::SpatialPoint padded(points.front());
        padded.resize(3);
        for (size_t direction = traits.getDimension(); direction < 3; ++direction)
          padded(direction) = Real(7 + direction);
        points.push_back(padded);
        for (const auto& point : points)
        {
          size_t next = 0;
          element.evaluateGradients(
            [&](size_t local, const auto& gradient) {
              EXPECT_EQ(local, next++);
              const auto expected = element.getBasis(local).getGradient()(point);
              if (local == 0)
              {
                size_t nestedCount = 0;
                element.evaluateGradients(
                  [&](size_t nested, const auto& value) {
                    EXPECT_EQ(nested, nestedCount++);
                    const auto reference =
                      element.getBasis(nested).getGradient()(points.back());
                    ASSERT_EQ(value.size(), reference.size());
                    for (size_t j = 0; j < traits.getDimension(); ++j)
                      EXPECT_EQ(value(j), reference(j));
                  },
                  points.back());
                EXPECT_EQ(nestedCount, element.getCount());
              }
              ASSERT_EQ(gradient.size(), expected.size());
              for (size_t j = 0; j < traits.getDimension(); ++j)
                EXPECT_EQ(gradient(j), expected(j));
            },
            point);
          EXPECT_EQ(next, element.getCount());
        }
        std::ostringstream diagnostic;
        EXPECT_THROW(element.evaluateGradients(
                       [&](size_t, const auto&) {
                         Alert::Exception(diagnostic)
                           << "Intentional consumer failure" << Alert::Raise;
                       },
                       points.front()),
          Alert::Exception);
        size_t repeated = 0;
        element.evaluateGradients(
          [&](size_t local, const auto& gradient) {
            EXPECT_EQ(local, repeated++);
            const auto expected = element.getBasis(local).getGradient()(points.front());
            for (size_t j = 0; j < traits.getDimension(); ++j)
              EXPECT_EQ(gradient(j), expected(j));
          },
          points.front());
        EXPECT_EQ(repeated, element.getCount());
      }
  };

  TEST_P(H1GradientEvaluationTest, MatchesPublicBasisAfterReentryAndFailure)
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
