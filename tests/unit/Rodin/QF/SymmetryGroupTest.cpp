/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Tests SymmetryGroup independently of any rule generator.
 *
 * The group and its orbit types are derived rather than tabulated, so they are
 * checked against things known independently: the order of each element's
 * symmetry group, the group axioms, and the orbit sizes Witherden and Vincent
 * tabulate domain by domain.
 */
#include <gtest/gtest.h>

#include <algorithm>
#include <random>

#include "Rodin/QF/SymmetryGroup.h"

#include "QuadratureInvariants.h"

using namespace Rodin;
using namespace Rodin::QF;
using namespace Rodin::Geometry;
using namespace Rodin::Tests::QF;

namespace
{
  /**
   * @brief Each element, the order of its symmetry group, and the orbit sizes
   * that group admits.
   *
   * The orders are those of the familiar groups --- the symmetric group of a
   * simplex's vertices, the dihedral group of a square, the octahedral group
   * of a cube --- and the orbit sizes are those tabulated in arXiv:1409.1865.
   */
  struct Expectation
  {
      const char* name;
      Polytope::Type type;
      size_t order;
      std::vector<size_t> orbitSizes;
  };

  const std::vector<Expectation>& expectations()
  {
    static const std::vector<Expectation> all = {
      {"segment", Polytope::Type::Segment, 2, {1, 2}},
      {"triangle", Polytope::Type::Triangle, 6, {1, 3, 6}},
      {"quadrilateral", Polytope::Type::Quadrilateral, 8, {1, 4, 4, 8}},
      {"tetrahedron", Polytope::Type::Tetrahedron, 24, {1, 4, 6, 12, 24}},
      {"wedge", Polytope::Type::Wedge, 12, {1, 2, 3, 6, 6, 12}},
      {"pyramid", Polytope::Type::Pyramid, 8, {1, 4, 4, 8}},
      {"hexahedron", Polytope::Type::Hexahedron, 48, {1, 6, 8, 12, 24, 24, 48}},
    };
    return all;
  }

  bool sameMap(const SymmetryGroup::Map& a, const SymmetryGroup::Map& b)
  {
    return (a.linear - b.linear).cwiseAbs().maxCoeff() < 1e-9 &&
      (a.translation - b.translation).cwiseAbs().maxCoeff() < 1e-9;
  }
}

/// @brief The derived group has the order the element's geometry requires.
TEST(SymmetryGroupTest, GroupOrderMatchesTheElement)
{
  for (const auto& expected : expectations())
  {
    const auto& group = SymmetryGroup::maps(expected.type);
    EXPECT_EQ(group.size(), expected.order) << expected.name;
  }
}

/**
 * @brief The derived set really is a group.
 *
 * It is built by a search over candidate maps, not by composing generators, so
 * closure and inverses are properties to be checked rather than assumed.
 */
TEST(SymmetryGroupTest, DerivedSetSatisfiesTheGroupAxioms)
{
  for (const auto& expected : expectations())
  {
    const auto& group = SymmetryGroup::maps(expected.type);
    ASSERT_FALSE(group.empty()) << expected.name;
    const Eigen::Index d = group.front().linear.rows();

    size_t identities = 0;
    for (const auto& a : group)
    {
      const bool isIdentity =
        (a.linear - Math::Matrix<Real>::Identity(d, d)).cwiseAbs().maxCoeff() < 1e-9 &&
        a.translation.cwiseAbs().maxCoeff() < 1e-9;
      identities += isIdentity;
    }
    EXPECT_EQ(identities, 1u) << expected.name << ": exactly one identity";

    for (const auto& a : group)
    {
      // Closure.
      for (const auto& b : group)
      {
        SymmetryGroup::Map composed;
        composed.linear = a.linear * b.linear;
        composed.translation = a.linear * b.translation + a.translation;
        const bool present = std::any_of(group.begin(), group.end(),
          [&](const SymmetryGroup::Map& c) { return sameMap(c, composed); });
        EXPECT_TRUE(present) << expected.name << ": not closed under composition";
      }

      // Inverses.
      SymmetryGroup::Map inverse;
      inverse.linear = a.linear.inverse();
      inverse.translation = -inverse.linear * a.translation;
      const bool present = std::any_of(group.begin(), group.end(),
        [&](const SymmetryGroup::Map& c) { return sameMap(c, inverse); });
      EXPECT_TRUE(present) << expected.name << ": missing an inverse";
    }
  }
}

/**
 * @brief Every symmetry maps the element onto itself.
 *
 * Checked on interior points rather than vertices: a map carrying the vertex
 * set onto itself but not the body would still pass a vertex-only check.
 */
TEST(SymmetryGroupTest, SymmetriesPreserveTheElement)
{
  std::mt19937 rng(20260824u);
  std::uniform_real_distribution<Real> uniform(0, 1);

  for (const auto& expected : expectations())
  {
    const size_t d = Polytope::Traits(expected.type).getDimension();
    for (size_t trial = 0; trial < 64; ++trial)
    {
      Math::SpatialVector<Real> x;
      x.resize(static_cast<Eigen::Index>(d));
      for (size_t attempt = 0;; ++attempt)
      {
        for (size_t k = 0; k < d; ++k)
          x[static_cast<Eigen::Index>(k)] = uniform(rng);
        if (isInside(expected.type, x, 0.0))
          break;
        ASSERT_LT(attempt, 10000u) << expected.name;
      }

      Math::Vector<Real> point(static_cast<Eigen::Index>(d));
      for (size_t k = 0; k < d; ++k)
        point(static_cast<Eigen::Index>(k)) = x[static_cast<Eigen::Index>(k)];

      for (const auto& map : SymmetryGroup::maps(expected.type))
      {
        const Math::Vector<Real> mapped = map(point);
        Math::SpatialVector<Real> image;
        image.resize(static_cast<Eigen::Index>(d));
        for (size_t k = 0; k < d; ++k)
          image[static_cast<Eigen::Index>(k)] = mapped(static_cast<Eigen::Index>(k));
        EXPECT_TRUE(isInside(expected.type, image, 1e-10))
          << expected.name << ": a symmetry carried an interior point outside";
      }
    }
  }
}

/**
 * @brief The derived orbit types are those Witherden and Vincent tabulate.
 *
 * This is the check that the derivation reproduces what the literature states
 * per domain, and in particular that the pyramid --- whose group fixes the
 * apex and permutes the base, and which the barycentric description of a
 * simplex cannot express --- comes out right without a path of its own.
 */
TEST(SymmetryGroupTest, OrbitTypesMatchThePublishedOnes)
{
  for (const auto& expected : expectations())
  {
    std::vector<size_t> sizes;
    for (const auto& stratum : SymmetryGroup::strata(expected.type))
      sizes.push_back(stratum.orbitSize);
    EXPECT_EQ(sizes, expected.orbitSizes)
      << expected.name << ": derived orbit sizes differ from the published ones";
  }
}

/**
 * @brief A generic point of each orbit type has exactly the stated orbit size,
 * and that size divides the group order.
 */
TEST(SymmetryGroupTest, OrbitSizesAreConsistentWithTheGroup)
{
  std::mt19937 rng(20260824u);
  std::uniform_real_distribution<Real> uniform(0.13, 0.29);

  for (const auto& expected : expectations())
  {
    const auto& group = SymmetryGroup::maps(expected.type);
    for (const auto& stratum : SymmetryGroup::strata(expected.type))
    {
      EXPECT_EQ(group.size() % stratum.orbitSize, 0u)
        << expected.name << ": orbit size must divide the group order";

      for (size_t trial = 0; trial < 8; ++trial)
      {
        Math::Vector<Real> probe = stratum.point;
        for (Eigen::Index k = 0; k < stratum.basis.cols(); ++k)
          probe += uniform(rng) * stratum.basis.col(k);
        EXPECT_EQ(SymmetryGroup::orbit(expected.type, probe).size(), stratum.orbitSize)
          << expected.name
          << ": a generic point of the stratum has the wrong "
             "orbit size";
      }
    }
  }
}
