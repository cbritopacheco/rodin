/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <set>

#include "Rodin/QF/SymmetricOrbit.h"

#include "QuadratureInvariants.h"

using namespace Rodin;
using namespace Rodin::QF;
using namespace Rodin::Geometry;
using namespace Rodin::Tests::QF;

namespace
{
  /// @brief Multiset of a tuple, rounded, so orbits can be compared as sets.
  std::multiset<long long> quantise(const SymmetricOrbit::Barycentric& b)
  {
    std::multiset<long long> s;
    for (const auto v : b)
      s.insert(static_cast<long long>(std::llround(v * 1e12)));
    return s;
  }
}

// --- Happy path: the conventional orbit classes and their cardinalities -----

/// @brief Triangle orbits S3, S21(a) and S111(a,b) have 1, 3 and 6 points.
/// The classes are not enumerated in the code; they are the multiplicity
/// patterns of the stored tuple, so this pins that expansion recovers them.
TEST(SymmetricOrbitTest, TriangleOrbitCardinalities)
{
  EXPECT_EQ(SymmetricOrbit({1. / 3, 1. / 3, 1. / 3}, 1).getSize(), 1u);
  EXPECT_EQ(SymmetricOrbit({0.2, 0.2, 0.6}, 1).getSize(), 3u);
  EXPECT_EQ(SymmetricOrbit({0.1, 0.3, 0.6}, 1).getSize(), 6u);
}

/// @brief Tetrahedron orbits S4, S31, S22, S211 and S1111 have 1, 4, 6, 12
/// and 24 points.
TEST(SymmetricOrbitTest, TetrahedronOrbitCardinalities)
{
  EXPECT_EQ(SymmetricOrbit({0.25, 0.25, 0.25, 0.25}, 1).getSize(), 1u);
  EXPECT_EQ(SymmetricOrbit({0.1, 0.1, 0.1, 0.7}, 1).getSize(), 4u);
  EXPECT_EQ(SymmetricOrbit({0.2, 0.2, 0.3, 0.3}, 1).getSize(), 6u);
  EXPECT_EQ(SymmetricOrbit({0.1, 0.1, 0.3, 0.5}, 1).getSize(), 12u);
  EXPECT_EQ(SymmetricOrbit({0.1, 0.2, 0.3, 0.4}, 1).getSize(), 24u);
}

// --- Invariant: the orbit is what it claims to be --------------------------

/// @brief The expansion is closed under permutation: permuting any member
/// yields another member. This is the property that makes the rule symmetric,
/// and it is what a taxonomy of named orbit classes would otherwise assert by
/// construction rather than prove.
TEST(SymmetricOrbitTest, ExpansionIsClosedUnderPermutation)
{
  for (const SymmetricOrbit& orbit : {
    SymmetricOrbit({0.1, 0.3, 0.6}, 1),
    SymmetricOrbit({0.1, 0.1, 0.3, 0.5}, 1),
    SymmetricOrbit({0.1, 0.2, 0.3, 0.4}, 1)})
  {
    const auto points = orbit.expand();
    std::set<std::multiset<long long>> patterns;
    for (const auto& p : points)
      patterns.insert(quantise(p));
    EXPECT_EQ(patterns.size(), 1u)
      << "every member must be a permutation of the same multiset";

    std::set<std::vector<long long>> distinct;
    for (const auto& p : points)
    {
      std::vector<long long> q;
      for (const auto v : p)
        q.push_back(static_cast<long long>(std::llround(v * 1e12)));
      distinct.insert(q);
    }
    EXPECT_EQ(distinct.size(), points.size()) << "members must be distinct";
  }
}

/// @brief Barycentric coordinates of every member sum to one.
TEST(SymmetricOrbitTest, ExpansionPreservesBarycentricSum)
{
  const SymmetricOrbit orbit({0.1, 0.2, 0.3, 0.4}, 1);
  for (const auto& p : orbit.expand())
  {
    Real s = 0;
    for (const auto v : p)
      s += v;
    EXPECT_NEAR(s, 1.0, 1e-15);
  }
}

// --- Interaction: the mapping uses Rodin's reference elements --------------

/// @brief A barycentric tuple maps onto the reference element published by
/// Polytope::Traits, and the unit tuples land exactly on its vertices.
TEST(SymmetricOrbitTest, ToReferenceLandsOnPolytopeTraitsVertices)
{
  for (const auto g : {Polytope::Type::Triangle, Polytope::Type::Tetrahedron})
  {
    const Polytope::Traits traits(g);
    const size_t nv = traits.getVertexCount();
    for (size_t i = 0; i < nv; ++i)
    {
      SymmetricOrbit::Barycentric b(nv, 0.0);
      b[i] = 1.0;
      const auto x = SymmetricOrbit::toReference(g, b);
      const auto& v = traits.getVertex(i);
      for (size_t k = 0; k < traits.getDimension(); ++k)
      {
        EXPECT_NEAR(x[static_cast<Eigen::Index>(k)],
                    v[static_cast<Eigen::Index>(k)], 1e-15)
          << "vertex " << i << " component " << k;
      }
    }
  }
}

/// @brief Every expanded point of an interior orbit lies inside the reference
/// element, checked against the half-space system Polytope::Traits publishes.
TEST(SymmetricOrbitTest, ExpandedPointsLieInsideReferenceElement)
{
  const SymmetricOrbit tri({0.1, 0.3, 0.6}, 1);
  for (const auto& b : tri.expand())
    EXPECT_TRUE(isInside(Polytope::Type::Triangle,
      SymmetricOrbit::toReference(Polytope::Type::Triangle, b), 1e-13));

  const SymmetricOrbit tet({0.1, 0.2, 0.3, 0.4}, 1);
  for (const auto& b : tet.expand())
    EXPECT_TRUE(isInside(Polytope::Type::Tetrahedron,
      SymmetricOrbit::toReference(Polytope::Type::Tetrahedron, b), 1e-13));
}

/// @brief The centroid orbit maps to the centroid Polytope::Traits publishes.
TEST(SymmetricOrbitTest, CentroidOrbitMapsToTraitsCentroid)
{
  for (const auto g : {Polytope::Type::Triangle, Polytope::Type::Tetrahedron})
  {
    const Polytope::Traits traits(g);
    const size_t nv = traits.getVertexCount();
    const SymmetricOrbit::Barycentric b(nv, 1.0 / static_cast<Real>(nv));
    const auto x = SymmetricOrbit::toReference(g, b);
    const auto& c = traits.getCentroid();
    for (size_t k = 0; k < traits.getDimension(); ++k)
      EXPECT_NEAR(x[static_cast<Eigen::Index>(k)],
                  c[static_cast<Eigen::Index>(k)], 1e-15);
  }
}

// --- Determinism -----------------------------------------------------------

/// @brief Expansion order is deterministic across calls.
TEST(SymmetricOrbitTest, ExpansionIsDeterministic)
{
  const SymmetricOrbit orbit({0.1, 0.2, 0.3, 0.4}, 1);
  EXPECT_EQ(orbit.expand(), orbit.expand());
}
