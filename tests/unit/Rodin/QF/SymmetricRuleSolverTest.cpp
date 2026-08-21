/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

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
          for (const auto& b : orbit.expand())
          {
            m_points.push_back(SymmetricOrbit::toReference(g, b));
            m_weights.push_back(orbit.getWeight());
          }
        }
      }

      size_t getSize() const { return m_weights.size(); }
      Real getWeight(size_t i) const { return m_weights[i]; }
      const Math::SpatialVector<Real>& getPoint(size_t i) const { return m_points[i]; }

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

  const std::vector<Case> cases = {
    {Polytope::Type::Triangle,    1, {{3}}},
    {Polytope::Type::Triangle,    2, {{2, 1}}},
    {Polytope::Type::Triangle,    3, {{2, 1}, {2, 1}}},
    {Polytope::Type::Triangle,    4, {{2, 1}, {2, 1}}},
    {Polytope::Type::Tetrahedron, 2, {{3, 1}}},
    {Polytope::Type::Tetrahedron, 3, {{3, 1}, {3, 1}}}
  };

  for (const auto& c : cases)
  {
    const auto r = SymmetricRuleSolver::solve(c.g, c.degree, c.config);
    ASSERT_TRUE(r.converged)
      << "degree " << c.degree << " residual " << r.residual;
    const OrbitRule rule(c.g, r.orbits);
    const auto report = exactnessSweep(rule, c.g, c.degree);
    EXPECT_LT(report.worstRelativeError, 1e-11)
      << "degree " << c.degree << " worst x^" << report.worstExponents[0]
      << " y^" << report.worstExponents[1] << " z^" << report.worstExponents[2];
  }
}

/// @brief A converged and admissible rule has positive weights and interior
/// points, which is what makes it usable for a barrier or a metric.
TEST(SymmetricRuleSolverTest, AdmissibleRulesArePositiveAndInterior)
{
  const auto r = SymmetricRuleSolver::solve(Polytope::Type::Triangle, 4,
    {{2, 1}, {2, 1}});
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
  const auto a = SymmetricRuleSolver::solve(Polytope::Type::Triangle, 4, {{2, 1}, {2, 1}});
  const auto b = SymmetricRuleSolver::solve(Polytope::Type::Triangle, 4, {{2, 1}, {2, 1}});
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
  const auto r = SymmetricRuleSolver::solve(Polytope::Type::Triangle, 2, {{3}},
    8 /* restarts */);
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
