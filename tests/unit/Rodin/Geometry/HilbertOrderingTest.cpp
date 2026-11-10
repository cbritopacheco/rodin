/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <Rodin/Geometry/HilbertOrdering.h>

using namespace Rodin;
using namespace Rodin::Geometry;

TEST(Geometry_HilbertOrdering, KnownTwoDimensionalCurve)
{
  LocalMesh::Builder builder;
  builder.initialize(2).nodes(16);
  for (size_t x = 0; x < 4; ++x)
  {
    for (size_t y = 0; y < 4; ++y)
      builder.vertex({Real(x), Real(y)});
  }
  HilbertOrdering ordering;
  ordering.setBitsPerCoordinate(2);
  EXPECT_EQ(ordering(builder.finalize()),
    (IndexVector{0, 4, 5, 1, 2, 3, 7, 6, 10, 11, 15, 14, 13, 9, 8, 12}));
}

TEST(Geometry_HilbertOrdering, CompleteGridsHaveFaceAdjacentSuccessors)
{
  for (size_t dim = 1; dim <= 3; ++dim)
  {
    const size_t count = size_t(1) << (2 * dim);
    LocalMesh::Builder builder;
    builder.initialize(dim).nodes(count);
    for (Index v = 0; v < count; ++v)
    {
      Math::SpatialPoint point(dim);
      for (size_t d = 0; d < dim; ++d)
        point(d) = Real((v >> (2 * d)) & 3);
      builder.vertex(point);
    }
    const auto mesh = builder.finalize();
    HilbertOrdering ordering;
    ordering.setBitsPerCoordinate(2);
    const auto order = ordering(mesh);
    auto sorted = order;
    std::sort(sorted.begin(), sorted.end());
    ASSERT_EQ(sorted.size(), count);
    for (Index v = 0; v < count; ++v)
      EXPECT_EQ(sorted[v], v);
    for (size_t i = 1; i < order.size(); ++i)
    {
      const auto delta =
        mesh.getVertexCoordinates(order[i]) - mesh.getVertexCoordinates(order[i - 1]);
      Real distance = 0;
      for (size_t d = 0; d < dim; ++d)
        distance += std::abs(delta(d));
      EXPECT_EQ(distance, 1) << "dimension " << dim << ", step " << i;
    }
  }
}

TEST(Geometry_HilbertOrdering, PrecisionEmptyAndCoincidentVertices)
{
  HilbertOrdering ordering;
  EXPECT_FALSE(ordering.getBitsPerCoordinate());
  EXPECT_THROW(ordering.setBitsPerCoordinate(0), Alert::Exception);
  EXPECT_THROW(ordering.setBitsPerCoordinate(53), Alert::Exception);
  LocalMesh::Builder empty;
  empty.initialize(3).nodes(0);
  const auto emptyMesh = empty.finalize();
  EXPECT_TRUE(ordering(emptyMesh).empty());
  ordering.setBitsPerCoordinate(22);
  EXPECT_THROW(ordering(emptyMesh), Alert::Exception);
  ordering.setBitsPerCoordinate(21);
  LocalMesh::Builder builder;
  builder.initialize(3).nodes(3).vertex({2, 2, 2}).vertex({2, 2, 2}).vertex({2, 2, 2});
  EXPECT_EQ(ordering(builder.finalize()), (IndexVector{0, 1, 2}));
  ordering.setBitsPerCoordinate(std::nullopt);
  EXPECT_FALSE(ordering.getBitsPerCoordinate());
}
