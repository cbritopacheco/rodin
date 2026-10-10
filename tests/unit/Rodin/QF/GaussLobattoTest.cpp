/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 */
#include <gtest/gtest.h>
#include <array>
#include <thread>
#include "Rodin/QF/GaussLobatto.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "QuadratureInvariants.h"

using namespace Rodin;

class GaussLobattoTest : public ::testing::TestWithParam<Geometry::Polytope::Type>
{};

TEST_P(GaussLobattoTest, PositiveNormalizedBoundaryInclusive)
{
  const auto geometry = GetParam();
  const Geometry::Polytope::Traits traits(geometry);
  for (size_t count = 2; count <= 12; ++count)
  {
    const QF::GaussLobatto rule(geometry, count);
    EXPECT_TRUE(Tests::QF::allWeightsPositive(rule));
    EXPECT_TRUE(Tests::QF::allPointsInside(rule, geometry));
    EXPECT_NEAR(Tests::QF::weightSum(rule), Tests::QF::referenceMeasure(geometry), 1e-13);
    for (size_t vertex = 0; vertex < traits.getVertexCount(); ++vertex)
    {
      bool found = false;
      for (size_t q = 0; q < rule.getSize(); ++q)
      {
        bool equal = true;
        for (size_t d = 0; d < traits.getDimension(); ++d)
          equal = equal && rule.getPoint(q)[d] == traits.getVertex(vertex)[d];
        found = found || equal;
      }
      EXPECT_TRUE(found) << "vertex=" << vertex << " count=" << count;
    }
    for (size_t i = 0; i < rule.getSize(); ++i)
    {
      for (size_t j = 0; j < i; ++j)
        EXPECT_GT((rule.getPoint(i) - rule.getPoint(j)).norm(), Real(0));
    }
  }
}

TEST_P(GaussLobattoTest, ExactPolynomialMoments)
{
  for (size_t count = 2; count <= 10; ++count)
  {
    const QF::GaussLobatto rule(GetParam(), count);
    const auto report = Tests::QF::exactnessSweep(rule, GetParam(), 2 * count - 3);
    EXPECT_LT(report.worstRelativeError, Real(2e-11)) << "count=" << count;
  }
}

TEST_P(GaussLobattoTest, DirectCachePreservesRulesAndIdentity)
{
  for (size_t degree = 0; degree <= 20; ++degree)
  {
    const size_t count = std::max<size_t>(2, (degree + 4) / 2);
    const auto& rule = QF::GaussLobatto::get(GetParam(), count);
    const QF::GaussLobatto owned(GetParam(), count);
    EXPECT_EQ(&rule, &QF::GaussLobatto::get(GetParam(), count));
    EXPECT_NE(&rule, &QF::PolytopeQuadratureFormula::get(degree, GetParam()));
    ASSERT_EQ(rule.getSize(), owned.getSize());
    for (size_t q = 0; q < rule.getSize(); ++q)
    {
      EXPECT_EQ(rule.getWeight(q), owned.getWeight(q));
      EXPECT_EQ((rule.getPoint(q) - owned.getPoint(q)).norm(), Real(0));
    }
    EXPECT_TRUE(Tests::QF::allWeightsPositive(rule));
    const auto report = Tests::QF::exactnessSweep(rule, GetParam(), degree);
    EXPECT_LT(report.worstRelativeError, Real(2e-11)) << "degree=" << degree;
  }
}

TEST(GaussLobatto, CacheSurvivesHotEvictionAndConcurrentAccess)
{
  const auto geometry = Geometry::Polytope::Type::Triangle;
  const auto* first = &QF::GaussLobatto::get(geometry, 3);
  const auto identity = first->getCacheIdentity();
  for (size_t count = 4; count <= 20; ++count)
    QF::GaussLobatto::get(geometry, count);
  EXPECT_EQ(first, &QF::GaussLobatto::get(geometry, 3));
  EXPECT_EQ(identity, first->getCacheIdentity());
  EXPECT_NE(first, &QF::GaussLobatto::get(geometry, 4));
  EXPECT_NE(first, &QF::GaussLobatto::get(Geometry::Polytope::Type::Tetrahedron, 3));

  std::array<const QF::GaussLobatto*, 4> rules{};
  std::array<std::thread, 4> threads;
  for (size_t i = 0; i < threads.size(); ++i)
  {
    threads[i] = std::thread([&, i] {
      rules[i] = &QF::GaussLobatto::get(geometry, 21);
    });
  }
  for (auto& thread : threads)
    thread.join();
  for (const auto* rule : rules)
    EXPECT_EQ(rule, rules.front());
}

INSTANTIATE_TEST_SUITE_P(AllGeometries, GaussLobattoTest, ::testing::Values(
  Geometry::Polytope::Type::Point,
  Geometry::Polytope::Type::Segment,
  Geometry::Polytope::Type::Triangle,
  Geometry::Polytope::Type::Quadrilateral,
  Geometry::Polytope::Type::Tetrahedron,
  Geometry::Polytope::Type::Hexahedron,
  Geometry::Polytope::Type::Wedge,
  Geometry::Polytope::Type::Pyramid));

TEST(GaussLobatto, SeparateTriangleOrdersAreRespected)
{
  const QF::GaussLobatto rule(Geometry::Polytope::Type::Triangle, 3, 5);
  EXPECT_EQ(rule.getSize(), size_t(11));
  EXPECT_LT(Tests::QF::exactnessSweep(rule,
    Geometry::Polytope::Type::Triangle, 3).worstRelativeError, Real(1e-12));
}

TEST(GaussLobatto, WeightedGaussianInteriorMoments)
{
  for (size_t alpha = 1; alpha <= 3; ++alpha)
  {
    for (size_t count = 1; count <= 12; ++count)
    {
      std::vector<Real> x, w;
      QF::GaussLegendre::gj1dUnit(count, alpha, 1, x, w);
      for (size_t degree = 0; degree <= 2 * count - 1; ++degree)
      {
        Real value = 0;
        for (size_t q = 0; q < count; ++q)
          value += w[q] * std::pow(x[q], degree);
        const Real exact = std::exp(std::lgamma(Real(degree + 2)) +
          std::lgamma(Real(alpha + 1)) - std::lgamma(Real(degree + alpha + 3)));
        EXPECT_NEAR(value, exact, 1e-14);
      }
    }
  }
}
