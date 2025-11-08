/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <Rodin/Geometry.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  // ============================================================================
  // Comprehensive Polytope Class Tests
  // ============================================================================

  TEST(Geometry_Polytope, BasicPolytopeAccess)
  {
    constexpr const size_t mdim = 2;
    
    Mesh mesh =
      Mesh<Context::Local>::Builder()
      .initialize(mdim)
      .nodes(3)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .finalize();

    mesh.getConnectivity().compute(mdim - 1, mdim);
    
    auto it = mesh.getBoundary();
    ASSERT_FALSE(it.end());
    
    // Access polytope via iterator
    Polytope edge = *it;
    EXPECT_EQ(edge.getDimension(), 1);
  }

  TEST(Geometry_Polytope, PolytopeIndexAccess)
  {
    constexpr const size_t mdim = 2;
    
    Mesh mesh =
      Mesh<Context::Local>::Builder()
      .initialize(mdim)
      .nodes(3)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .finalize();

    mesh.getConnectivity().compute(mdim - 1, mdim);
    
    auto it = mesh.getBoundary();
    ASSERT_FALSE(it.end());
    
    Polytope edge = *it;
    Index idx = edge.getIndex();
    
    EXPECT_GE(idx, 0);
  }

  TEST(Geometry_Polytope, PolytopeAttributeAccess)
  {
    constexpr const size_t mdim = 2;
    
    Mesh mesh =
      Mesh<Context::Local>::Builder()
      .initialize(mdim)
      .nodes(3)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .finalize();

    mesh.getConnectivity().compute(mdim - 1, mdim);
    
    auto it = mesh.getBoundary();
    ASSERT_FALSE(it.end());
    
    Polytope edge = *it;
    Attribute attr = edge.getAttribute();
    
    EXPECT_GE(attr, 0);
  }

  TEST(Geometry_Polytope, PolytopeDimension)
  {
    constexpr const size_t mdim = 2;
    
    Mesh mesh =
      Mesh<Context::Local>::Builder()
      .initialize(mdim)
      .nodes(3)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .finalize();

    mesh.getConnectivity().compute(mdim - 1, mdim);
    
    // Boundary polytopes should be 1D (edges)
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
    {
      Polytope edge = *it;
      EXPECT_EQ(edge.getDimension(), 1);
    }
  }

  TEST(Geometry_Polytope, PolytopeMeshReference)
  {
    constexpr const size_t mdim = 2;
    
    Mesh mesh =
      Mesh<Context::Local>::Builder()
      .initialize(mdim)
      .nodes(3)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .finalize();

    mesh.getConnectivity().compute(mdim - 1, mdim);
    
    auto it = mesh.getBoundary();
    ASSERT_FALSE(it.end());
    
    Polytope edge = *it;
    const MeshBase& edgeMesh = edge.getMesh();
    
    EXPECT_EQ(&edgeMesh, &mesh);
  }

  TEST(Geometry_Polytope, PolytopeCopyConstruction)
  {
    constexpr const size_t mdim = 2;
    
    Mesh mesh =
      Mesh<Context::Local>::Builder()
      .initialize(mdim)
      .nodes(3)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .finalize();

    mesh.getConnectivity().compute(mdim - 1, mdim);
    
    auto it = mesh.getBoundary();
    ASSERT_FALSE(it.end());
    
    Polytope p1 = *it;
    Polytope p2(p1);
    
    EXPECT_EQ(p1.getIndex(), p2.getIndex());
    EXPECT_EQ(p1.getDimension(), p2.getDimension());
  }

  TEST(Geometry_Polytope, PolytopeCopyAssignment)
  {
    constexpr const size_t mdim = 2;
    
    Mesh mesh =
      Mesh<Context::Local>::Builder()
      .initialize(mdim)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .polytope(Polytope::Type::Triangle, {1, 3, 2})
      .finalize();

    mesh.getConnectivity().compute(mdim - 1, mdim);
    
    auto it = mesh.getBoundary();
    ASSERT_FALSE(it.end());
    
    Polytope p1 = *it;
    ++it;
    
    if (!it.end())
    {
      Polytope p2 = *it;
      p2 = p1;  // Copy assignment
      
      EXPECT_EQ(p1.getIndex(), p2.getIndex());
    }
  }

  TEST(Geometry_Polytope, PolytopeMoveConstruction)
  {
    constexpr const size_t mdim = 2;
    
    Mesh mesh =
      Mesh<Context::Local>::Builder()
      .initialize(mdim)
      .nodes(3)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .finalize();

    mesh.getConnectivity().compute(mdim - 1, mdim);
    
    auto it = mesh.getBoundary();
    ASSERT_FALSE(it.end());
    
    Polytope p1 = *it;
    Index idx = p1.getIndex();
    
    Polytope p2(std::move(p1));
    EXPECT_EQ(p2.getIndex(), idx);
  }

  TEST(Geometry_Polytope, PolytopeMoveAssignment)
  {
    constexpr const size_t mdim = 2;
    
    Mesh mesh =
      Mesh<Context::Local>::Builder()
      .initialize(mdim)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .polytope(Polytope::Type::Triangle, {1, 3, 2})
      .finalize();

    mesh.getConnectivity().compute(mdim - 1, mdim);
    
    auto it = mesh.getBoundary();
    ASSERT_FALSE(it.end());
    
    Polytope p1 = *it;
    Index idx = p1.getIndex();
    
    ++it;
    if (!it.end())
    {
      Polytope p2 = *it;
      p2 = std::move(p1);  // Move assignment
      
      EXPECT_EQ(p2.getIndex(), idx);
    }
  }

  TEST(Geometry_Polytope, PolytopeEquality)
  {
    constexpr const size_t mdim = 2;
    
    Mesh mesh =
      Mesh<Context::Local>::Builder()
      .initialize(mdim)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .polytope(Polytope::Type::Triangle, {1, 3, 2})
      .finalize();

    mesh.getConnectivity().compute(mdim - 1, mdim);
    
    auto it1 = mesh.getBoundary();
    ASSERT_FALSE(it1.end());
    
    Polytope p1 = *it1;
    Polytope p1_copy = p1;
    
    EXPECT_TRUE(p1 == p1_copy);
    
    ++it1;
    if (!it1.end())
    {
      Polytope p2 = *it1;
      EXPECT_FALSE(p1 == p2);
    }
  }

  TEST(Geometry_Polytope, PolytopeInequality)
  {
    constexpr const size_t mdim = 2;
    
    Mesh mesh =
      Mesh<Context::Local>::Builder()
      .initialize(mdim)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .polytope(Polytope::Type::Triangle, {1, 3, 2})
      .finalize();

    mesh.getConnectivity().compute(mdim - 1, mdim);
    
    auto it1 = mesh.getBoundary();
    ASSERT_FALSE(it1.end());
    
    Polytope p1 = *it1;
    Polytope p1_copy = p1;
    
    EXPECT_FALSE(p1 != p1_copy);
    
    ++it1;
    if (!it1.end())
    {
      Polytope p2 = *it1;
      EXPECT_TRUE(p1 != p2);
    }
  }

  TEST(Geometry_Polytope, Polytope3D_Boundary)
  {
    constexpr const size_t mdim = 3;
    
    Mesh mesh =
      Mesh<Context::Local>::Builder()
      .initialize(mdim)
      .nodes(4)
      .vertex({0, 0, 0})
      .vertex({1, 0, 0})
      .vertex({0, 1, 0})
      .vertex({0, 0, 1})
      .polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 3})
      .finalize();

    mesh.getConnectivity().compute(mdim - 1, mdim);
    
    // Boundary polytopes should be 2D (faces)
    size_t count = 0;
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
    {
      Polytope face = *it;
      EXPECT_EQ(face.getDimension(), 2);
      ++count;
    }
    
    EXPECT_EQ(count, 4);  // Tetrahedron has 4 faces
  }

  TEST(Geometry_Polytope, Polytope3D_Interface)
  {
    constexpr const size_t mdim = 3;
    
    // Create two tetrahedra sharing a face
    Mesh mesh =
      Mesh<Context::Local>::Builder()
      .initialize(mdim)
      .nodes(5)
      .vertex({0, 0, 0})
      .vertex({1, 0, 0})
      .vertex({0, 1, 0})
      .vertex({0, 0, 1})
      .vertex({0, 0, -1})
      .polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 3})
      .polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 4})
      .finalize();

    mesh.getConnectivity().compute(mdim - 1, mdim);
    
    // Interface polytopes should be 2D (shared faces)
    size_t count = 0;
    for (auto it = mesh.getInterface(); !it.end(); ++it)
    {
      Polytope face = *it;
      EXPECT_EQ(face.getDimension(), 2);
      ++count;
    }
    
    EXPECT_EQ(count, 1);  // One shared face
  }

  TEST(Geometry_Polytope, PolytopeBoundaryCount_Triangle)
  {
    constexpr const size_t mdim = 2;
    
    Mesh mesh =
      Mesh<Context::Local>::Builder()
      .initialize(mdim)
      .nodes(3)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .finalize();

    mesh.getConnectivity().compute(mdim - 1, mdim);
    
    size_t count = 0;
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
    {
      ++count;
    }
    
    EXPECT_EQ(count, 3);  // Triangle has 3 edges
  }

  TEST(Geometry_Polytope, PolytopeBoundaryCount_Quadrilateral)
  {
    constexpr const size_t mdim = 2;
    
    Mesh mesh =
      Mesh<Context::Local>::Builder()
      .initialize(mdim)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({1, 1})
      .vertex({0, 1})
      .polytope(Polytope::Type::Quadrilateral, {0, 1, 2, 3})
      .finalize();

    mesh.getConnectivity().compute(mdim - 1, mdim);
    
    size_t count = 0;
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
    {
      ++count;
    }
    
    EXPECT_EQ(count, 4);  // Quadrilateral has 4 edges
  }

  TEST(Geometry_Polytope, PolytopeInterfaceCount_TwoTriangles)
  {
    constexpr const size_t mdim = 2;
    
    Mesh mesh =
      Mesh<Context::Local>::Builder()
      .initialize(mdim)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .polytope(Polytope::Type::Triangle, {1, 3, 2})
      .finalize();

    mesh.getConnectivity().compute(mdim - 1, mdim);
    
    size_t count = 0;
    for (auto it = mesh.getInterface(); !it.end(); ++it)
    {
      ++count;
    }
    
    EXPECT_EQ(count, 1);  // One shared edge
  }

  TEST(Geometry_Polytope, PolytopeConsistentDimension_BoundaryVsInterface)
  {
    constexpr const size_t mdim = 2;
    
    Mesh mesh =
      Mesh<Context::Local>::Builder()
      .initialize(mdim)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .polytope(Polytope::Type::Triangle, {1, 3, 2})
      .finalize();

    mesh.getConnectivity().compute(mdim - 1, mdim);
    
    // All boundary and interface polytopes should have same dimension
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
    {
      EXPECT_EQ(it->getDimension(), 1);
    }
    
    for (auto it = mesh.getInterface(); !it.end(); ++it)
    {
      EXPECT_EQ(it->getDimension(), 1);
    }
  }

  TEST(Geometry_Polytope, PolytopeIndexUniqueness)
  {
    constexpr const size_t mdim = 2;
    
    Mesh mesh =
      Mesh<Context::Local>::Builder()
      .initialize(mdim)
      .nodes(3)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .finalize();

    mesh.getConnectivity().compute(mdim - 1, mdim);
    
    std::set<Index> indices;
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
    {
      Polytope edge = *it;
      indices.insert(edge.getIndex());
    }
    
    EXPECT_EQ(indices.size(), 3);  // All indices should be unique
  }

  TEST(Geometry_Polytope, PolytopeComplexMesh)
  {
    constexpr const size_t mdim = 2;
    
    // Create a more complex mesh with multiple triangles
    Mesh mesh =
      Mesh<Context::Local>::Builder()
      .initialize(mdim)
      .nodes(5)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .vertex({0.5, 0.5})
      .polytope(Polytope::Type::Triangle, {0, 1, 4})
      .polytope(Polytope::Type::Triangle, {1, 3, 4})
      .polytope(Polytope::Type::Triangle, {3, 2, 4})
      .polytope(Polytope::Type::Triangle, {2, 0, 4})
      .finalize();

    mesh.getConnectivity().compute(mdim - 1, mdim);
    
    // Count boundary edges (exterior edges)
    size_t boundary_count = 0;
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
    {
      EXPECT_EQ(it->getDimension(), 1);
      ++boundary_count;
    }
    
    // Count interface edges (shared between cells)
    size_t interface_count = 0;
    for (auto it = mesh.getInterface(); !it.end(); ++it)
    {
      EXPECT_EQ(it->getDimension(), 1);
      ++interface_count;
    }
    
    EXPECT_GT(boundary_count, 0);
    EXPECT_GT(interface_count, 0);
  }

  TEST(Geometry_Polytope, PolytopeAccessMultipleTimes)
  {
    constexpr const size_t mdim = 2;
    
    Mesh mesh =
      Mesh<Context::Local>::Builder()
      .initialize(mdim)
      .nodes(3)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .finalize();

    mesh.getConnectivity().compute(mdim - 1, mdim);
    
    // Access same polytope multiple times
    auto it1 = mesh.getBoundary();
    auto it2 = mesh.getBoundary();
    
    ASSERT_FALSE(it1.end());
    ASSERT_FALSE(it2.end());
    
    Polytope p1 = *it1;
    Polytope p2 = *it2;
    
    EXPECT_EQ(p1.getIndex(), p2.getIndex());
    EXPECT_TRUE(p1 == p2);
  }
}
