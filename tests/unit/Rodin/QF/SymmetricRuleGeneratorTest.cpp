/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Tests SymmetricRuleGenerator: the generated rules, and the reasoning
 * that produces them.
 *
 * The rules are checked as objects in their own right --- exact, positive,
 * interior, symmetric, no larger than the published counts --- and separately
 * the machinery that decides what to search is checked against an independent
 * computation. Both matter: a generator can be wrong by returning a bad rule,
 * or by never looking where the good one is.
 */
#include <gtest/gtest.h>

#include <algorithm>
#include <random>

#include "Rodin/QF/SymmetricRuleGenerator.h"

#include "PublishedCounts.h"
#include "QuadratureInvariants.h"

using namespace Rodin;
using namespace Rodin::QF;
using namespace Rodin::Geometry;
using namespace Rodin::Tests::QF;

namespace
{
  struct Element
  {
      const char* name;
      Polytope::Type type;
      size_t degrees;   ///< Strengths to check, from one upward.
  };

  /// @brief Kept modest so the suite stays a unit test; the higher strengths
  /// are exercised by the generation driver.
  const std::vector<Element>& elements()
  {
    static const std::vector<Element> all = {
      {"triangle", Polytope::Type::Triangle, 8},
      {"quadrilateral", Polytope::Type::Quadrilateral, 8},
      {"tetrahedron", Polytope::Type::Tetrahedron, 5},
      {"wedge", Polytope::Type::Wedge, 5},
      {"pyramid", Polytope::Type::Pyramid, 5},
      {"hexahedron", Polytope::Type::Hexahedron, 5},
    };
    return all;
  }
}

/**
 * @brief The invariant dimension from Molien's series agrees with a direct
 * count.
 *
 * This number decides which decompositions the search will even attempt, so
 * getting it wrong does not produce a wrong rule --- it produces no rule, or a
 * needlessly large one, with nothing to indicate why. The direct count builds
 * the symmetrisation of every basis mode at many points and takes the rank of
 * the result, which is the dimension of the invariant subspace by definition.
 */
TEST(SymmetricRuleGeneratorTest, InvariantDimensionMatchesADirectCount)
{
  std::mt19937 rng(20260824u);
  std::uniform_real_distribution<Real> uniform(0, 1);

  for (const auto& element : elements())
  {
    const size_t d = Polytope::Traits(element.type).getDimension();
    const auto& group = SymmetryGroup::maps(element.type);
    for (const size_t degree : {1u, 2u, 3u, 4u})
    {
      const auto basis = CollapsedBasis::indices(element.type, degree);
      const size_t samples = 4 * basis.size() + 16;

      Math::Matrix<Real> symmetrised(
        static_cast<Eigen::Index>(basis.size()), static_cast<Eigen::Index>(samples));
      for (size_t q = 0; q < samples; ++q)
      {
        Math::SpatialVector<Real> x;
        x.resize(static_cast<Eigen::Index>(d));
        for (size_t attempt = 0;; ++attempt)
        {
          for (size_t k = 0; k < d; ++k)
            x[static_cast<Eigen::Index>(k)] = uniform(rng);
          if (isInside(element.type, x, 0.0))
            break;
          ASSERT_LT(attempt, 10000u) << element.name;
        }

        // The orbit of the sample, then each mode summed over it.
        std::vector<Math::SpatialVector<Real>> images;
        Math::Vector<Real> point(static_cast<Eigen::Index>(d));
        for (size_t k = 0; k < d; ++k)
          point(static_cast<Eigen::Index>(k)) = x[static_cast<Eigen::Index>(k)];
        for (const auto& map : group)
        {
          const Math::Vector<Real> mapped = map(point);
          Math::SpatialVector<Real> image;
          image.resize(static_cast<Eigen::Index>(d));
          for (size_t k = 0; k < d; ++k)
            image[static_cast<Eigen::Index>(k)] = mapped(static_cast<Eigen::Index>(k));
          images.push_back(std::move(image));
        }

        for (size_t e = 0; e < basis.size(); ++e)
        {
          Real sum = 0;
          for (const auto& image : images)
            sum += CollapsedBasis::evaluate(element.type, basis[e], image);
          symmetrised(static_cast<Eigen::Index>(e), static_cast<Eigen::Index>(q)) = sum;
        }
      }

      Eigen::ColPivHouseholderQR<Math::Matrix<Real>> qr(symmetrised);
      qr.setThreshold(1e-9);
      EXPECT_EQ(static_cast<size_t>(qr.rank()),
        SymmetricRuleGenerator::invariantDimension(element.type, degree))
        << element.name << ", strength " << degree;
    }
  }
}

/**
 * @brief Generated rules are exact, positive, and inside the element, and no
 * larger than the published counts.
 */
TEST(SymmetricRuleGeneratorTest, GeneratedRulesAreValidAndNoLargerThanPublished)
{
  for (const auto& element : elements())
  {
    const size_t d = Polytope::Traits(element.type).getDimension();
    for (size_t degree = 1; degree <= element.degrees; ++degree)
    {
      const auto rule = SymmetricRuleGenerator::search(element.type, degree, 64);
      ASSERT_TRUE(rule.converged)
        << element.name << ", strength " << degree << ": no rule found";
      ASSERT_TRUE(rule.admissible) << element.name << ", strength " << degree;

      const auto report = exactnessSweep(rule, element.type, degree);
      EXPECT_LT(report.worstRelativeError, 1e-12)
        << element.name << ", strength " << degree;
      EXPECT_TRUE(allWeightsPositive(rule)) << element.name;
      EXPECT_TRUE(allPointsInside(rule, element.type, 1e-12)) << element.name;

      // No rule can have fewer points than the counting bound; one that does
      // is not a discovery but a moment system that has stopped spanning.
      EXPECT_GE(rule.getSize(), countingBound(d, degree))
        << element.name << ", strength " << degree << ": below the counting bound";

      const size_t published =
        publishedCount(witherdenVincentCounts(element.type), degree);
      if (published > 0)
        EXPECT_LE(rule.getSize(), published) << element.name << ", strength " << degree
                                             << ": more points than the published rule";
    }
  }
}

/**
 * @brief Generated rules really are symmetric.
 *
 * Symmetry is what the whole construction assumes, and it is a property of the
 * points themselves rather than of the search: applying any symmetry of the
 * element must permute the rule, carrying each point to a point of equal
 * weight. A rule failing this would still integrate correctly on the
 * polynomials tested and be wrong in general.
 */
TEST(SymmetricRuleGeneratorTest, GeneratedRulesAreInvariantUnderTheGroup)
{
  for (const auto& element : elements())
  {
    const size_t d = Polytope::Traits(element.type).getDimension();
    for (size_t degree = 1; degree <= std::min<size_t>(element.degrees, 5); ++degree)
    {
      const auto rule = SymmetricRuleGenerator::search(element.type, degree, 64);
      ASSERT_TRUE(rule.converged && rule.admissible) << element.name;

      for (const auto& map : SymmetryGroup::maps(element.type))
      {
        for (size_t i = 0; i < rule.getSize(); ++i)
        {
          Math::Vector<Real> point(static_cast<Eigen::Index>(d));
          for (size_t k = 0; k < d; ++k)
            point(static_cast<Eigen::Index>(k)) =
              rule.getPoint(i)[static_cast<Eigen::Index>(k)];
          const Math::Vector<Real> image = map(point);

          bool matched = false;
          for (size_t j = 0; j < rule.getSize() && !matched; ++j)
          {
            Real distance = 0;
            for (size_t k = 0; k < d; ++k)
              distance = std::max(distance,
                std::abs(image(static_cast<Eigen::Index>(k)) -
                  rule.getPoint(j)[static_cast<Eigen::Index>(k)]));
            if (distance < 1e-9)
            {
              matched = true;
              EXPECT_NEAR(rule.getWeight(i), rule.getWeight(j), 1e-12)
                << element.name << ": symmetric points carry different weights";
            }
          }
          EXPECT_TRUE(matched) << element.name << ", strength " << degree
                               << ": a symmetry carried a point off the rule";
        }
      }
    }
  }
}

/// @brief The same request yields the same rule, so inlined coefficients do
/// not drift when regenerated.
TEST(SymmetricRuleGeneratorTest, GenerationIsDeterministic)
{
  for (const auto& element : elements())
  {
    const auto first = SymmetricRuleGenerator::search(element.type, 4, 64);
    const auto second = SymmetricRuleGenerator::search(element.type, 4, 64);
    ASSERT_EQ(first.getSize(), second.getSize()) << element.name;
    for (size_t i = 0; i < first.getSize(); ++i)
    {
      EXPECT_EQ(first.getWeight(i), second.getWeight(i)) << element.name;
      for (Eigen::Index k = 0; k < first.getPoint(i).size(); ++k)
        EXPECT_EQ(first.getPoint(i)[k], second.getPoint(i)[k]) << element.name;
    }
  }
}
