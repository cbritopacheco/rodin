/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Constraint-map state and mesh-region iteration contracts. */
#include <gtest/gtest.h>
#include <array>
#include <set>
#include <type_traits>
#include <vector>
#include "Rodin/Assembly/ConstraintMap.h"
#include "Rodin/Assembly/Sequential.h"

using namespace Rodin;

template <class Scalar>
class ConstraintContracts : public ::testing::Test
{};
using ConstraintScalars = ::testing::Types<Real, Complex>;
TYPED_TEST_SUITE(ConstraintContracts, ConstraintScalars);

/// @brief Rejected constraints preserve existing fixed and identified rows.
TYPED_TEST(ConstraintContracts, RejectionsPreserveRows)
{
  using Scalar = TypeParam;
  using Map = Assembly::ConstraintMap<Scalar>;
  using Entries = std::vector<typename Map::Entry>;
  Map map(4);
  const Entries row{{0, Scalar(2)}};
  map.setIdentification(1, row, Scalar(3));
  map.setFixed(2, Scalar(4));
  EXPECT_THROW(map.setFixed(1, Scalar(9)), Alert::Exception);
  EXPECT_THROW(map.setIdentification(1, row), Alert::Exception);
  EXPECT_THROW(map.setIdentification(2, row), Alert::Exception);
  for (const Index index : {-1, 4})
  {
    EXPECT_THROW(map.setFixed(index, Scalar(1)), Alert::Exception);
    EXPECT_THROW(map.setIdentification(index, row), Alert::Exception);
    EXPECT_THROW(map.getIdentificationValue(index), Alert::Exception);
    EXPECT_FALSE(map.isFixed(index));
    EXPECT_FALSE(map.isIdentified(index));
  }
  EXPECT_THROW(map.setIdentification(3, Entries{{4, Scalar(1)}}), Alert::Exception);
  EXPECT_THROW(
    map.setIdentification(3, Entries{{3, Scalar(1)}}, Scalar(1)), Alert::Exception);
  EXPECT_THROW(
    map.setIdentification(3, Entries{{3, Scalar(1)}, {0, Scalar(2)}}), Alert::Exception);
  ASSERT_EQ(map.expand(1).size(), 1u);
  EXPECT_EQ(map.expand(1)[0].index, 0);
  EXPECT_EQ(map.expand(1)[0].coefficient, Scalar(2));
  EXPECT_EQ(map.getIdentificationValue(1), Scalar(3));
  EXPECT_EQ(map.getFixedValue(2), Scalar(4));
  EXPECT_FALSE(map.isIdentified(3));
  EXPECT_FALSE(map.isFixed(3));
  Map cycle(2);
  cycle.setIdentification(0, Entries{{1, Scalar(1)}});
  cycle.setIdentification(1, Entries{{0, Scalar(1)}});
  EXPECT_THROW(cycle.finalize(), Alert::Exception);
  cycle.reset(2);
  EXPECT_TRUE(cycle.getIdentifiedRows().empty());
  cycle.setFixed(0, Scalar(2));
  cycle.finalize();
  EXPECT_EQ(cycle.getFixedValue(0), Scalar(2));
}

/// @brief Flattening carries affine offsets through chains and fixed masters.
TYPED_TEST(ConstraintContracts, AffineFlatteningAndIndependentCopies)
{
  using Scalar = TypeParam;
  using Map = Assembly::ConstraintMap<Scalar>;
  using Entries = std::vector<typename Map::Entry>;
  const Scalar factor = [] {
    if constexpr (std::is_same_v<Scalar, Complex>)
      return Scalar(1, -2);
    else
      return Scalar(1);
  }();
  const auto [alpha, beta, expectedCoefficient, expectedOffset, expectedFirst,
    expectedSecond] = [] {
    if constexpr (std::is_same_v<Scalar, Complex>)
      return std::array<Scalar, 6>{Scalar(2, 1), Scalar(-0.5, 0.5), Scalar(-1.5, 0.5),
        Scalar(5.5, -3.5), Scalar(31, -27), Scalar(2, 21)};
    else
      return std::array<Scalar, 6>{2, -0.5, -1, 2.5, 17, -4.5};
  }();
  Map original(4);
  original.setIdentification(
    1, Entries{{0, alpha + Scalar(1)}, {0, Scalar(-1)}}, Scalar(3) * factor);
  original.setIdentification(2, Entries{{1, beta}}, Scalar(4) * factor);
  original.setIdentification(3, Entries{{3, Scalar(2)}}, Scalar(6) * factor);
  for (size_t repeat = 0; repeat < 2; ++repeat)
  {
    original.finalize();
    ASSERT_EQ(original.expand(2).size(), 1u);
    EXPECT_EQ(original.expand(2)[0].index, 0);
    EXPECT_EQ(original.expand(2)[0].coefficient, expectedCoefficient);
    EXPECT_EQ(original.getIdentificationValue(2), expectedOffset);
    EXPECT_EQ(original.getFixedValue(3), Scalar(-6) * factor);
  }
  auto copied = original;
  copied.setFixed(0, Scalar(7) * factor);
  copied.finalize();
  ASSERT_TRUE(copied.isFixed(1));
  ASSERT_TRUE(copied.isFixed(2));
  EXPECT_EQ(copied.getFixedValue(1), expectedFirst);
  EXPECT_EQ(copied.getFixedValue(2), expectedSecond);
  EXPECT_FALSE(original.isFixed(0));
  EXPECT_TRUE(original.isIdentified(1));
  EXPECT_TRUE(original.isIdentified(2));
  Map cancelled(2);
  cancelled.setIdentification(
    1, Entries{{0, Scalar(1)}, {0, Scalar(-1)}}, Scalar(5) * factor);
  ASSERT_TRUE(cancelled.isFixed(1));
  EXPECT_EQ(cancelled.getFixedValue(1), Scalar(5) * factor);
}

/// @brief Iterator and indexed traversal select the same cells, faces and boundary/interface partition.
TEST(AssemblyRegionContracts, IteratorMatchesIndexedSelection)
{
  using namespace Geometry;
  auto mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3});
  mesh.getConnectivity().compute(1, 2);
  for (const auto region :
    {Region::Cells, Region::Faces, Region::Boundary, Region::Interface})
  {
    Assembly::SequentialIteration traversal(mesh, region);
    const size_t dimension = region == Region::Cells ? 2 : 1;
    EXPECT_EQ(traversal.getDimension(), dimension);
    EXPECT_EQ(traversal.getCount(), mesh.getPolytopeCount(dimension));
    std::set<Index> iterated, indexed;
    for (auto it = traversal.getIterator(); it; ++it)
      iterated.insert(it->getIndex());
    for (Index index = 0; index < static_cast<Index>(traversal.getCount()); ++index)
    {
      if (traversal.filter(index))
      {
        const auto polytope = traversal.getPolytope(index);
        EXPECT_EQ(polytope.getDimension(), dimension);
        EXPECT_EQ(polytope.getIndex(), index);
        indexed.insert(index);
      }
    }
    EXPECT_EQ(iterated, indexed);
    EXPECT_EQ(iterated.size(), region == Region::Faces ? 16u : 8u);
  }
}
