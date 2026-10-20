/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include "Rodin/Assembly/OpenMP.h"

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  /// @brief Indexed assembly candidates match interval topology and each region's filter.
  TEST(OpenMPIteration, IndexedDescriptorsAndFiltersMatchConnectedInterval)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {4});
    mesh.getConnectivity().compute(0, 1);
    for (const auto region :
      {Region::Cells, Region::Faces, Region::Boundary, Region::Interface})
    {
      SCOPED_TRACE(static_cast<int>(region));
      Assembly::OpenMPIteration iteration(mesh, region);
      const bool cells = region == Region::Cells;
      const size_t expectedCount = cells ? 3 : 4;
      EXPECT_EQ(iteration.getCount(), expectedCount);
      EXPECT_EQ(iteration.getDimension(), cells ? 1u : 0u);
      size_t selected = 0;
      for (Index i = 0; i < expectedCount; ++i)
      {
        SCOPED_TRACE(i);
        const bool boundary = i == 0 || i == 3;
        const bool expected = region == Region::Boundary ? boundary
          : region == Region::Interface                  ? !boundary
                                                         : true;
        EXPECT_EQ(iteration.filter(i), expected);
        selected += iteration.filter(i);
        const auto value = iteration.getPolytope(i);
        auto iterator = iteration.getIterator(i);
        ASSERT_TRUE(iterator);
        EXPECT_EQ(value.getDimension(), cells ? 1u : 0u);
        EXPECT_EQ(iterator->getDimension(), cells ? 1u : 0u);
        EXPECT_EQ(value.getIndex(), i);
        EXPECT_EQ(iterator->getIndex(), i);
        EXPECT_EQ(&value.getMesh(), &mesh);
        EXPECT_EQ(&iterator->getMesh(), &mesh);
        EXPECT_EQ(
          value.getGeometry(), cells ? Polytope::Type::Segment : Polytope::Type::Point);
        const auto vertices = value.getVertices();
        ASSERT_EQ(vertices.size(), cells ? 2 : 1);
        EXPECT_EQ(vertices[0], i);
        if (cells)
          EXPECT_EQ(vertices[1], i + 1);
        else
        {
          EXPECT_EQ(mesh.getVertexCoordinates(iterator->getVertices()[0])(0),
            static_cast<Real>(i));
        }
      }
      EXPECT_EQ(selected, cells ? 3u : region == Region::Faces ? 4u : 2u);
    }
  }

  /// @brief A finalized empty mesh has zero-dimensional cell traversal with no candidates.
  TEST(OpenMPIteration, EmptyCellsHaveNoCandidates)
  {
    auto mesh = LocalMesh::Builder().initialize(1).nodes(0).finalize();
    ASSERT_EQ(mesh.getDimension(), 0u);
    Assembly::OpenMPIteration iteration(mesh, Region::Cells);
    EXPECT_EQ(iteration.getCount(), 0u);
    EXPECT_EQ(iteration.getDimension(), 0u);
  }
}
