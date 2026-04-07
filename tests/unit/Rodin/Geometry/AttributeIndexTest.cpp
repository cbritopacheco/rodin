/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Configure.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  // ---- Mesh-level attribute tests (tests AttributeIndex indirectly) ----

  TEST(Geometry_AttributeIndex, SetAndGetCellAttribute)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    const size_t D = mesh.getDimension();

    // Set attribute on cell 0
    mesh.setAttribute({D, 0}, 42);
    auto attr = mesh.getPolytope(D, 0)->getAttribute();
    EXPECT_TRUE(attr.has_value());
    EXPECT_EQ(*attr, 42);
  }

  TEST(Geometry_AttributeIndex, SetAndGetMultipleCells)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    const size_t D = mesh.getDimension();

    for (Index i = 0; i < mesh.getCellCount(); ++i)
      mesh.setAttribute({D, i}, i % 3 + 1);

    for (Index i = 0; i < mesh.getCellCount(); ++i)
    {
      auto attr = mesh.getPolytope(D, i)->getAttribute();
      EXPECT_TRUE(attr.has_value());
      EXPECT_EQ(*attr, i % 3 + 1);
    }
  }

  TEST(Geometry_AttributeIndex, GetMultipleAttributes)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    const size_t D = mesh.getDimension();

    // Set different attributes
    for (Index i = 0; i < mesh.getCellCount(); ++i)
      mesh.setAttribute({D, i}, i % 3 + 1);

    // Verify we can read back all the individual attributes
    FlatSet<Attribute> uniqueAttrs;
    for (Index i = 0; i < mesh.getCellCount(); ++i)
    {
      auto attr = mesh.getAttribute(D, i);
      EXPECT_TRUE(attr.has_value());
      uniqueAttrs.insert(*attr);
    }
    // Should have 3 unique attributes (1, 2, 3)
    EXPECT_EQ(uniqueAttrs.size(), 3);
  }

  TEST(Geometry_AttributeIndex, SetBoundaryAttributes)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(1, 2);
    const size_t D = mesh.getDimension();

    // Set attributes on boundary faces
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
      mesh.setAttribute({D - 1, it->getIndex()}, 10);

    // Verify they are set
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
    {
      auto attr = it->getAttribute();
      EXPECT_TRUE(attr.has_value());
      EXPECT_EQ(*attr, 10);
    }
  }

  TEST(Geometry_AttributeIndex, VertexAttributes)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});

    // Set attribute on vertex 0
    mesh.setAttribute({0, 0}, 99);
    auto attr = mesh.getPolytope(0, 0)->getAttribute();
    EXPECT_TRUE(attr.has_value());
    EXPECT_EQ(*attr, 99);
  }
}
