/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include "Rodin/Geometry.h"
#include "Rodin/MMG.h"

using namespace Rodin;
using namespace Rodin::Geometry;

/// @brief Native copies retain the capacity and free lists needed for insertion.
TEST(MMG5Copy, InsertionAfterSourceDestruction)
{
  using Native = MMG::MMG5;
  using Owner = std::unique_ptr<MMG5_Mesh, decltype(&Native::destroyMesh)>;
  for (const bool volume : {false, true})
  {
    SCOPED_TRACE(volume);
    auto mesh = volume ? LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {2, 2, 2})
                       : LocalMesh::UniformGrid(Polytope::Type::Triangle, {2, 2});
    Owner source(Native::rodinToMesh(mesh), Native::destroyMesh);
    Owner copy(
      static_cast<MMG5_pMesh>(mycalloc(1, sizeof(MMG5_Mesh))), Native::destroyMesh);
    Native::copyMesh(source.get(), copy.get());
    ASSERT_EQ(copy->npmax, source->npmax);
    ASSERT_GT(copy->npmax, copy->np);
    const auto oldCount = copy->np;
    source.reset();
    double coordinates[3] = {0.25, 0.25, 0.25};
    const auto vertex = volume ? MMG3D_newPt(copy.get(), coordinates, 0, 0)
                               : MMG2D_newPt(copy.get(), coordinates, 0);
    ASSERT_GT(vertex, oldCount);
    EXPECT_DOUBLE_EQ(copy->point[vertex].c[0], 0.25);
    EXPECT_DOUBLE_EQ(copy->point[vertex].c[1], 0.25);
    if (volume)
    {
      EXPECT_GT(MMG3D_newElt(copy.get()), mesh.getCellCount());
    }
    else
    {
      EXPECT_GT(MMG2D_newElt(copy.get()), mesh.getCellCount());
    }
  }
}
