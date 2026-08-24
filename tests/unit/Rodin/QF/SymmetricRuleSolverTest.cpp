/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <numeric>

#include "Rodin/QF/SymmetricRuleSolver.h"

#include "QuadratureInvariants.h"

using namespace Rodin;
using namespace Rodin::QF;
using namespace Rodin::Geometry;
using namespace Rodin::Tests::QF;

namespace
{
  /// @brief Presents a solved orbit set with the interface the exactness
  /// sweep expects, so a generated rule is measured by the same instrument
  /// that measures a shipped one.
  class OrbitRule
  {
    public:
      OrbitRule(Polytope::Type g, const std::vector<SymmetricOrbit>& orbits)
      {
        for (const auto& orbit : orbits)
        {
          for (auto& p : orbit.expandPoints(g))
          {
            m_points.push_back(std::move(p));
            m_weights.push_back(orbit.getWeight());
          }
        }
      }

      size_t getSize() const
      {
        return m_weights.size();
      }
      Real getWeight(size_t i) const
      {
        return m_weights[i];
      }
      const Math::SpatialVector<Real>& getPoint(size_t i) const
      {
        return m_points[i];
      }

    private:
      std::vector<Math::SpatialVector<Real>> m_points;
      std::vector<Real> m_weights;
  };
}

// --- Happy path: known rules, in closed form -------------------------------

/// @brief Degree 1 on the triangle is the centroid rule: one point, weight
/// equal to the area. This one is forced, so it is a genuine known-answer
/// check rather than a self-consistency check.
TEST(SymmetricRuleSolverTest, TriangleDegreeOneIsTheCentroidRule)
{
  const auto r = SymmetricRuleSolver::solve(Polytope::Type::Triangle, 1, {{3}});
  ASSERT_TRUE(r.converged) << "residual " << r.residual;
  ASSERT_EQ(r.orbits.size(), 1u);
  EXPECT_EQ(r.orbits[0].getSize(), 1u);
  EXPECT_NEAR(r.orbits[0].getWeight(), 0.5, 1e-13);
  for (const auto v : r.orbits[0].getBarycentric())
    EXPECT_NEAR(v, 1.0 / 3.0, 1e-13);
}

/// @brief Degree 2 on the triangle is the three-point rule with barycentric
/// parameter 1/6 and weight 1/6, the classical S21(1/6) rule. Its parameters
/// are determined, so this pins the solver to a value it cannot have invented.
TEST(SymmetricRuleSolverTest, TriangleDegreeTwoIsTheThreePointRule)
{
  const auto r = SymmetricRuleSolver::solve(Polytope::Type::Triangle, 2, {{2, 1}});
  ASSERT_TRUE(r.converged) << "residual " << r.residual;
  ASSERT_EQ(r.orbits.size(), 1u);
  EXPECT_EQ(r.orbits[0].getSize(), 3u);
  EXPECT_NEAR(r.orbits[0].getWeight(), 1.0 / 6.0, 1e-12);

  auto b = r.orbits[0].getBarycentric();
  std::sort(b.begin(), b.end());
  EXPECT_NEAR(b[0], 1.0 / 6.0, 1e-10);
  EXPECT_NEAR(b[1], 1.0 / 6.0, 1e-10);
  EXPECT_NEAR(b[2], 2.0 / 3.0, 1e-10);
}

// --- Invariant: whatever it returns as converged really is exact ------------

/// @brief Any rule the solver reports as converged integrates every monomial
/// of its degree exactly, measured by the independent moment oracle. This is
/// the property that lets generated coefficients be trusted without a
/// published table to compare against.
TEST(SymmetricRuleSolverTest, ConvergedRulesPassTheExactnessSweep)
{
  struct Case
  {
      Polytope::Type g;
      size_t degree;
      SymmetricRuleSolver::Configuration config;
  };

  const std::vector<Case> cases = {{Polytope::Type::Triangle, 1, {{3}}},
    {Polytope::Type::Triangle, 2, {{2, 1}}},
    {Polytope::Type::Triangle, 3, {{2, 1}, {2, 1}}},
    {Polytope::Type::Triangle, 4, {{2, 1}, {2, 1}}},
    {Polytope::Type::Tetrahedron, 2, {{3, 1}}},
    {Polytope::Type::Tetrahedron, 3, {{3, 1}, {3, 1}}}};

  for (const auto& c : cases)
  {
    const auto r = SymmetricRuleSolver::solve(c.g, c.degree, c.config);
    ASSERT_TRUE(r.converged) << "degree " << c.degree << " residual " << r.residual;
    const OrbitRule rule(c.g, r.orbits);
    const auto report = exactnessSweep(rule, c.g, c.degree);
    EXPECT_LT(report.worstRelativeError, 1e-11)
      << "degree " << c.degree << " worst x^" << report.worstExponents[0] << " y^"
      << report.worstExponents[1] << " z^" << report.worstExponents[2];
  }
}

/// @brief A converged and admissible rule has positive weights and interior
/// points, which is what makes it usable for a barrier or a metric.
TEST(SymmetricRuleSolverTest, AdmissibleRulesArePositiveAndInterior)
{
  const auto r =
    SymmetricRuleSolver::solve(Polytope::Type::Triangle, 4, {{2, 1}, {2, 1}});
  ASSERT_TRUE(r.converged) << "residual " << r.residual;
  ASSERT_TRUE(r.admissible);
  const OrbitRule rule(Polytope::Type::Triangle, r.orbits);
  EXPECT_TRUE(allWeightsPositive(rule));
  EXPECT_TRUE(allPointsInside(rule, Polytope::Type::Triangle));
  EXPECT_NEAR(weightAmplification(rule), 1.0, 1e-14);
}

// --- Determinism: regenerated coefficients must not drift -------------------

/// @brief The same request yields bitwise identical coefficients. Coefficients
/// are inlined from this solver, so a regeneration that drifted would silently
/// change every assembled integral in the library.
TEST(SymmetricRuleSolverTest, SolveIsDeterministic)
{
  const auto a =
    SymmetricRuleSolver::solve(Polytope::Type::Triangle, 4, {{2, 1}, {2, 1}});
  const auto b =
    SymmetricRuleSolver::solve(Polytope::Type::Triangle, 4, {{2, 1}, {2, 1}});
  ASSERT_EQ(a.orbits.size(), b.orbits.size());
  for (size_t i = 0; i < a.orbits.size(); ++i)
  {
    EXPECT_EQ(a.orbits[i].getWeight(), b.orbits[i].getWeight());
    EXPECT_EQ(a.orbits[i].getBarycentric(), b.orbits[i].getBarycentric());
  }
}

// --- Negative and boundary -------------------------------------------------

/// @brief An under-parameterised configuration cannot reach the tolerance.
/// A single centroid orbit has one unknown and cannot satisfy degree 2, so
/// the solver must report failure rather than a plausible wrong rule.
TEST(SymmetricRuleSolverTest, UnderparameterisedConfigurationDoesNotConverge)
{
  const auto r =
    SymmetricRuleSolver::solve(Polytope::Type::Triangle, 2, {{3}}, 8 /* restarts */);
  EXPECT_FALSE(r.converged);
  EXPECT_GT(r.residual, 1e-6);
}

/// @brief The moment oracle used by the solver agrees with the one used by the
/// tests, which were written independently.
TEST(SymmetricRuleSolverTest, SolverMomentsAgreeWithTestOracle)
{
  for (size_t a = 0; a <= 4; ++a)
  {
    for (size_t b = 0; a + b <= 4; ++b)
    {
      EXPECT_NEAR(SymmetricRuleSolver::simplexMoment({a, b}),
        exactMoment(Polytope::Type::Triangle, a, b, 0), 1e-16)
        << "triangle x^" << a << " y^" << b;
      for (size_t c = 0; a + b + c <= 4; ++c)
      {
        EXPECT_NEAR(SymmetricRuleSolver::simplexMoment({a, b, c}),
          exactMoment(Polytope::Type::Tetrahedron, a, b, c), 1e-16)
          << "tet x^" << a << " y^" << b << " z^" << c;
      }
    }
  }
}

/// @brief Unknown counting matches the layout toOrbits consumes.
TEST(SymmetricRuleSolverTest, UnknownCountMatchesParameterLayout)
{
  EXPECT_EQ(SymmetricRuleSolver::unknownCount({{3}}), 1u);
  EXPECT_EQ(SymmetricRuleSolver::unknownCount({{2, 1}}), 2u);
  EXPECT_EQ(SymmetricRuleSolver::unknownCount({{1, 1, 1}}), 3u);
  EXPECT_EQ(SymmetricRuleSolver::unknownCount({{2, 1}, {2, 1}}), 4u);
}

// --- Configuration search --------------------------------------------------

/// @brief patternsFor enumerates the orbit classes of each simplex: the
/// partitions of d+1. Triangle {3},{2,1},{1,1,1}; tetrahedron {4},{3,1},
/// {2,2},{2,1,1},{1,1,1,1}.
TEST(SymmetricRuleSolverTest, PatternsForEnumeratesTheOrbitClasses)
{
  const auto tri = SymmetricRuleSolver::patternsFor(Polytope::Type::Triangle);
  EXPECT_EQ(tri.size(), 3u);
  const auto tet = SymmetricRuleSolver::patternsFor(Polytope::Type::Tetrahedron);
  EXPECT_EQ(tet.size(), 5u);
  for (const auto& p : tri)
    EXPECT_EQ(std::accumulate(p.begin(), p.end(), size_t{0}), 3u);
  for (const auto& p : tet)
    EXPECT_EQ(std::accumulate(p.begin(), p.end(), size_t{0}), 4u);
}

/// @brief patternSize is the multinomial coefficient of the pattern, i.e. the
/// cardinality of the orbit it names.
TEST(SymmetricRuleSolverTest, PatternSizeMatchesOrbitCardinality)
{
  EXPECT_EQ(SymmetricRuleSolver::patternSize({3}), 1u);
  EXPECT_EQ(SymmetricRuleSolver::patternSize({2, 1}), 3u);
  EXPECT_EQ(SymmetricRuleSolver::patternSize({1, 1, 1}), 6u);
  EXPECT_EQ(SymmetricRuleSolver::patternSize({4}), 1u);
  EXPECT_EQ(SymmetricRuleSolver::patternSize({3, 1}), 4u);
  EXPECT_EQ(SymmetricRuleSolver::patternSize({2, 2}), 6u);
  EXPECT_EQ(SymmetricRuleSolver::patternSize({2, 1, 1}), 12u);
  EXPECT_EQ(SymmetricRuleSolver::patternSize({1, 1, 1, 1}), 24u);
  EXPECT_EQ(SymmetricRuleSolver::configurationSize({{3}, {2, 1}, {2, 1}}), 7u);
}

/// @brief The search recovers the published minimal point counts on the
/// triangle. These are the figure of merit of the Xiao-Gimbutas and
/// Witherden-Vincent families, so agreeing with them is evidence that the
/// search finds the same rules and not merely valid ones.
TEST(SymmetricRuleSolverTest, SearchRecoversPublishedTrianglePointCounts)
{
  const std::vector<std::pair<size_t, size_t>> expected = {
    {1, 1}, {2, 3}, {3, 6}, {4, 6}, {5, 7}};
  for (const auto& [degree, points] : expected)
  {
    SymmetricRuleSolver::Configuration config;
    const auto r =
      SymmetricRuleSolver::search(Polytope::Type::Triangle, degree, 24, 96, &config);
    ASSERT_TRUE(r.converged && r.admissible)
      << "degree " << degree << " residual " << r.residual;
    EXPECT_EQ(SymmetricRuleSolver::configurationSize(config), points)
      << "degree " << degree;
  }
}

/// @brief Known-answer: the degree-5 seven-point triangle rule carries
/// centroid weight 9/80, a determined value the search cannot have invented.
TEST(SymmetricRuleSolverTest, TriangleDegreeFiveCentroidWeightIsNineEightieths)
{
  const auto r = SymmetricRuleSolver::search(Polytope::Type::Triangle, 5, 24);
  ASSERT_TRUE(r.converged && r.admissible);
  bool found = false;
  for (const auto& orbit : r.orbits)
  {
    if (orbit.getSize() == 1u)
    {
      EXPECT_NEAR(orbit.getWeight(), 9.0 / 80.0, 1e-9);
      found = true;
    }
  }
  EXPECT_TRUE(found) << "the degree-5 rule must contain a centroid orbit";
}

/// @brief Known-answer: the degree-4 six-point triangle rule has the
/// classical barycentric parameters 0.445948490915965 and
/// 0.091576213509771, with weights 0.111690794839005 and 0.054975871827661.
TEST(SymmetricRuleSolverTest, TriangleDegreeFourMatchesClassicalParameters)
{
  const auto r = SymmetricRuleSolver::search(Polytope::Type::Triangle, 4, 24);
  ASSERT_TRUE(r.converged && r.admissible);
  ASSERT_EQ(r.orbits.size(), 2u);

  std::vector<std::pair<Real, Real>> got;
  for (const auto& orbit : r.orbits)
  {
    auto b = orbit.getBarycentric();
    std::sort(b.begin(), b.end());
    got.emplace_back(b.front(), orbit.getWeight());
  }
  std::sort(got.begin(), got.end());

  EXPECT_NEAR(got[0].first, 0.091576213509771, 1e-9);
  EXPECT_NEAR(got[0].second, 0.054975871827661, 1e-9);
  EXPECT_NEAR(got[1].first, 0.108103018168070, 1e-9);
  EXPECT_NEAR(got[1].second, 0.111690794839005, 1e-9);
}

/// @brief Every rule the search returns as admissible passes the independent
/// exactness sweep, on both simplices.
TEST(SymmetricRuleSolverTest, SearchedRulesAreExactOnBothSimplices)
{
  for (const auto g : {Polytope::Type::Triangle, Polytope::Type::Tetrahedron})
  {
    for (size_t degree = 1; degree <= 4; ++degree)
    {
      const auto r = SymmetricRuleSolver::search(g, degree, 24);
      ASSERT_TRUE(r.converged && r.admissible)
        << "degree " << degree << " residual " << r.residual;
      const OrbitRule rule(g, r.orbits);
      EXPECT_LT(exactnessSweep(rule, g, degree).worstRelativeError, 1e-11)
        << "degree " << degree;
      EXPECT_TRUE(allWeightsPositive(rule));
      EXPECT_TRUE(allPointsInside(rule, g));
      EXPECT_NEAR(weightSum(rule), referenceMeasure(g), 1e-12);
    }
  }
}

/// @brief Boundary: a point budget below the minimal rule yields no
/// admissible configuration rather than a wrong one.
TEST(SymmetricRuleSolverTest, SearchRespectsThePointBudget)
{
  const auto r = SymmetricRuleSolver::search(Polytope::Type::Triangle, 5, 4);
  EXPECT_FALSE(r.converged && r.admissible);
}

/// @brief The search is deterministic.
TEST(SymmetricRuleSolverTest, SearchIsDeterministic)
{
  SymmetricRuleSolver::Configuration ca, cb;
  const auto a = SymmetricRuleSolver::search(Polytope::Type::Triangle, 4, 24, 96, &ca);
  const auto b = SymmetricRuleSolver::search(Polytope::Type::Triangle, 4, 24, 96, &cb);
  EXPECT_EQ(ca, cb);
  ASSERT_EQ(a.orbits.size(), b.orbits.size());
  for (size_t i = 0; i < a.orbits.size(); ++i)
  {
    EXPECT_EQ(a.orbits[i].getWeight(), b.orbits[i].getWeight());
    EXPECT_EQ(a.orbits[i].getBarycentric(), b.orbits[i].getBarycentric());
  }
}

// --- Wedge: the product orbit model ----------------------------------------

/// @brief A wedge orbit is the product of a triangle orbit with its tensor
/// part, so its cardinality multiplies: S21(a) on the mid-plane is 3 points,
/// and the same orbit on a reflected pair is 6.
TEST(SymmetricRuleSolverTest, WedgeOrbitCardinalitiesMultiply)
{
  const SymmetricOrbit mid({0.2, 0.2, 0.6}, {0.5}, 1);
  EXPECT_EQ(mid.getSize(), 3u);
  const SymmetricOrbit pair({0.2, 0.2, 0.6}, {0.3, 0.7}, 1);
  EXPECT_EQ(pair.getSize(), 6u);
  const SymmetricOrbit centroidPair({1. / 3, 1. / 3, 1. / 3}, {0.3, 0.7}, 1);
  EXPECT_EQ(centroidPair.getSize(), 2u);

  using OrbitClass = SymmetricRuleSolver::OrbitClass;
  EXPECT_EQ(SymmetricRuleSolver::classSize(OrbitClass({2, 1}, true)), 3u);
  EXPECT_EQ(SymmetricRuleSolver::classSize(OrbitClass({2, 1}, false)), 6u);
  EXPECT_EQ(SymmetricRuleSolver::classSize(OrbitClass({3}, false)), 2u);
}

/// @brief expandPoints places a wedge orbit inside the reference wedge, which
/// is the unit triangle crossed with [0,1] as published by Polytope::Traits.
TEST(SymmetricRuleSolverTest, WedgeOrbitPointsLieInsideReferenceWedge)
{
  const SymmetricOrbit orbit({0.1, 0.3, 0.6}, {0.25, 0.75}, 1);
  const auto points = orbit.expandPoints(Polytope::Type::Wedge);
  EXPECT_EQ(points.size(), 12u);
  for (const auto& p : points)
  {
    EXPECT_EQ(p.size(), 3);
    EXPECT_TRUE(isInside(Polytope::Type::Wedge, p, 1e-13));
  }
}

/// @brief The wedge moment used by the solver agrees with the independently
/// written test oracle. The wedge is a product, so its moment is the triangle
/// moment times 1/(c+1); getting that factor wrong would still produce a rule
/// that looked exact against itself.
TEST(SymmetricRuleSolverTest, WedgeMomentAgreesWithTestOracle)
{
  for (size_t a = 0; a <= 4; ++a)
    for (size_t b = 0; a + b <= 4; ++b)
      for (size_t c = 0; a + b + c <= 4; ++c)
        EXPECT_NEAR(
          SymmetricRuleSolver::referenceMoment(Polytope::Type::Wedge, {a, b, c}),
          exactMoment(Polytope::Type::Wedge, a, b, c), 1e-16)
          << "x^" << a << " y^" << b << " z^" << c;
}

/// @brief The search recovers the published Witherden-Vincent prism point
/// counts. These are the minima their paper reports, so agreeing with them is
/// evidence the product orbit model is the same one they use.
TEST(SymmetricRuleSolverTest, SearchRecoversPublishedWedgePointCounts)
{
  const std::vector<std::pair<size_t, size_t>> expected = {
    {1, 1}, {2, 5}, {3, 8}, {4, 11}};
  for (const auto& [degree, points] : expected)
  {
    SymmetricRuleSolver::Configuration config;
    const auto r =
      SymmetricRuleSolver::search(Polytope::Type::Wedge, degree, 16, 96, &config);
    ASSERT_TRUE(r.converged && r.admissible)
      << "degree " << degree << " residual " << r.residual;
    EXPECT_EQ(SymmetricRuleSolver::configurationSize(config), points)
      << "degree " << degree;
  }
}

/// @brief Wedge rules the search returns are exact, positive and interior,
/// measured by the same instrument used for the simplices.
TEST(SymmetricRuleSolverTest, SearchedWedgeRulesAreExact)
{
  for (size_t degree = 1; degree <= 4; ++degree)
  {
    const auto r = SymmetricRuleSolver::search(Polytope::Type::Wedge, degree, 16);
    ASSERT_TRUE(r.converged && r.admissible)
      << "degree " << degree << " residual " << r.residual;
    const OrbitRule rule(Polytope::Type::Wedge, r.orbits);
    const auto report = exactnessSweep(rule, Polytope::Type::Wedge, degree);
    EXPECT_LT(report.worstRelativeError, 1e-11)
      << "degree " << degree << " worst x^" << report.worstExponents[0] << " y^"
      << report.worstExponents[1] << " z^" << report.worstExponents[2];
    EXPECT_TRUE(allWeightsPositive(rule));
    EXPECT_TRUE(allPointsInside(rule, Polytope::Type::Wedge));
    EXPECT_NEAR(weightSum(rule), referenceMeasure(Polytope::Type::Wedge), 1e-12);
  }
}

/// @brief Wedge orbits placed on the mid-plane carry the reflection-invariant
/// value exactly, and reflected pairs are symmetric about it.
TEST(SymmetricRuleSolverTest, WedgeTensorPartsRespectTheReflection)
{
  const auto r = SymmetricRuleSolver::search(Polytope::Type::Wedge, 4, 16);
  ASSERT_TRUE(r.converged && r.admissible);
  for (const auto& orbit : r.orbits)
  {
    const auto& t = orbit.getTensor();
    ASSERT_FALSE(t.empty()) << "every wedge orbit carries a tensor part";
    if (t.size() == 1)
      EXPECT_NEAR(t[0], 0.5, 1e-15);
    else
    {
      ASSERT_EQ(t.size(), 2u);
      EXPECT_NEAR(t[0] + t[1], 1.0, 1e-14);
    }
  }
}

/// @brief The wedge search is deterministic.
TEST(SymmetricRuleSolverTest, WedgeSearchIsDeterministic)
{
  SymmetricRuleSolver::Configuration ca, cb;
  const auto a = SymmetricRuleSolver::search(Polytope::Type::Wedge, 3, 16, 96, &ca);
  const auto b = SymmetricRuleSolver::search(Polytope::Type::Wedge, 3, 16, 96, &cb);
  EXPECT_EQ(ca, cb);
  ASSERT_EQ(a.orbits.size(), b.orbits.size());
  for (size_t i = 0; i < a.orbits.size(); ++i)
  {
    EXPECT_EQ(a.orbits[i].getWeight(), b.orbits[i].getWeight());
    EXPECT_EQ(a.orbits[i].getTensor(), b.orbits[i].getTensor());
  }
}

/**
 * @brief The analytic Jacobian of the eliminated system agrees with a
 * difference quotient.
 *
 * A Jacobian that is quietly wrong does not announce itself: the search still
 * runs, converges more slowly, and returns rules a few digits short of exact.
 * The comparison sweeps the difference step and takes the best agreement,
 * because a single step proves nothing --- too large and truncation dominates,
 * too small and roundoff does.
 */
TEST(SymmetricRuleSolverTest, AnalyticJacobianMatchesDifferenceQuotient)
{
  using Pattern = SymmetricRuleSolver::Pattern;
  const Pattern t3{3}, t21{2, 1}, t111{1, 1, 1};
  const Pattern q4{4}, q31{3, 1}, q22{2, 2}, q211{2, 1, 1};

  struct Case
  {
      const char* name;
      Polytope::Type g;
      size_t degree;
      SymmetricRuleSolver::Configuration config;
  };
  const std::vector<Case> cases = {
    {"triangle, two edge orbits", Polytope::Type::Triangle, 4, {t21, t21}},
    {"triangle, general orbit", Polytope::Type::Triangle, 8, {t21, t21, t111}},
    {"triangle, high degree", Polytope::Type::Triangle, 12, {t3, t21, t111, t111}},
    {"tetrahedron, face orbits", Polytope::Type::Tetrahedron, 5, {q31, q31, q22}},
    {"tetrahedron, general orbit", Polytope::Type::Tetrahedron, 7, {q4, q31, q211}},
    {"wedge, mid-plane orbit", Polytope::Type::Wedge, 4, {{t21, true}}},
    {"wedge, reflected pair", Polytope::Type::Wedge, 4, {{t21, false}}},
    {"wedge, mixed orbits", Polytope::Type::Wedge, 5, {{t111, false}, {t3, true}}},
  };

  for (const auto& c : cases)
  {
    const Real worst =
      SymmetricRuleSolver::verifyJacobian(c.g, c.degree, c.config);
    EXPECT_LT(worst, 1e-6) << "analytic Jacobian disagrees with the difference "
                              "quotient on " << c.name;
  }
}

/**
 * @brief Searched triangle rules are exact, positive and interior, and reach
 * the published point counts.
 *
 * The counts are those of Xiao and Gimbutas. Matching them is the point of the
 * exercise; the exactness sweep is what says the match is real rather than a
 * rule that merely has the right number of points.
 */
TEST(SymmetricRuleSolverTest, SearchedTriangleRulesMatchPublishedCounts)
{
  const std::vector<size_t> published = {0, 1, 3, 6, 6, 7, 12, 15, 16};
  for (size_t degree = 1; degree < published.size(); ++degree)
  {
    const auto result =
      SymmetricRuleSolver::search(Polytope::Type::Triangle, degree, 40, 256);
    ASSERT_TRUE(result.converged) << "degree " << degree;
    ASSERT_TRUE(result.admissible) << "degree " << degree;

    const OrbitRule rule(Polytope::Type::Triangle, result.orbits);
    EXPECT_EQ(rule.getSize(), published[degree])
      << "degree " << degree << " should reach the published count";

    const auto report = exactnessSweep(rule, Polytope::Type::Triangle, degree);
    EXPECT_LT(report.worstRelativeError, 1e-12) << "degree " << degree;
    EXPECT_TRUE(allWeightsPositive(rule)) << "degree " << degree;
    EXPECT_TRUE(allPointsInside(rule, Polytope::Type::Triangle))
      << "degree " << degree;
  }
}
