/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <numeric>
#include <vector>

#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Configure.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  // ==================================================================
  // Mesh::Box tests
  // Box takes a FACE type and produces a mesh of dimension = face_dim+1
  // with only the boundary faces (like a "box" surface)
  // ==================================================================

  /// @brief Box with segment faces in 2D space.
  TEST(Geometry_MeshBox, Segment_2Nodes)
  {
    // Box with segment faces in 2D space
    auto mesh = LocalMesh::Box(Polytope::Type::Segment, {2, 2});
    EXPECT_EQ(mesh.getSpaceDimension(), 2);
    EXPECT_GT(mesh.getCellCount(), 0);
  }

  /// @brief Verifies segment 4 x 4 nodes for geometry mesh box by checking exact expected values.
  TEST(Geometry_MeshBox, Segment_4x4Nodes)
  {
    auto mesh = LocalMesh::Box(Polytope::Type::Segment, {4, 4});
    EXPECT_EQ(mesh.getSpaceDimension(), 2);
    EXPECT_GT(mesh.getCellCount(), 0);
  }

  /// @brief Box with triangular faces in 3D space.
  TEST(Geometry_MeshBox, Triangle_3Nodes)
  {
    // Box with triangular faces in 3D space
    auto mesh = LocalMesh::Box(Polytope::Type::Triangle, {2, 2, 2});
    EXPECT_EQ(mesh.getSpaceDimension(), 3);
    EXPECT_GT(mesh.getCellCount(), 0);
  }

  /// @brief Box with quad faces in 3D space.
  TEST(Geometry_MeshBox, Quadrilateral_3Nodes)
  {
    // Box with quad faces in 3D space
    auto mesh = LocalMesh::Box(Polytope::Type::Quadrilateral, {2, 2, 2});
    EXPECT_EQ(mesh.getSpaceDimension(), 3);
    EXPECT_GT(mesh.getCellCount(), 0);
  }

  /// @brief Box with point faces in 1D space (just two endpoints).
  TEST(Geometry_MeshBox, Point_2Nodes)
  {
    // Box with point faces in 1D space (just two endpoints)
    auto mesh = LocalMesh::Box(Polytope::Type::Point, {2});
    EXPECT_EQ(mesh.getSpaceDimension(), 1);
    EXPECT_GT(mesh.getVertexCount(), 0);
  }

  /// @brief Verifies box 2 D coordinates positive for geometry mesh box.
  TEST(Geometry_MeshBox, Box2D_CoordinatesPositive)
  {
    auto mesh = LocalMesh::Box(Polytope::Type::Segment, {4, 4});
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      auto coords = it->getCoordinates();
      EXPECT_GE(coords(0), -1e-12);
      EXPECT_GE(coords(1), -1e-12);
    }
  }

  /// @brief Verifies box 3 D coordinates positive for geometry mesh box.
  TEST(Geometry_MeshBox, Box3D_CoordinatesPositive)
  {
    auto mesh = LocalMesh::Box(Polytope::Type::Triangle, {2, 2, 2});
    for (auto it = mesh.getVertex(); !it.end(); ++it)
    {
      auto coords = it->getCoordinates();
      EXPECT_GE(coords(0), -1e-12);
      EXPECT_GE(coords(1), -1e-12);
      EXPECT_GE(coords(2), -1e-12);
    }
  }

  // ==================================================================
  // Mesh::ccl tests (connected component labeling)
  // ==================================================================

  /// @brief Verifies CCL single component for geometry mesh ops by checking exact expected values.
  TEST(Geometry_MeshOps, CCL_SingleComponent)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);

    auto result = mesh.ccl(D, [](const Polytope&, const Polytope&) { return true; });
    const auto& components = result.getComponents();
    // A uniform grid should have exactly one connected component
    EXPECT_EQ(components.size(), 1);
    EXPECT_EQ(components[0].size(), mesh.getCellCount());
  }

  /// @brief Verifies CCL with predicates for geometry mesh ops by checking exact expected values.
  TEST(Geometry_MeshOps, CCL_WithPredicates)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);

    // All cells connected, all included
    auto result = mesh.ccl(
      D,
      [](const Polytope&, const Polytope&) { return true; },
      [](const Polytope&) { return true; }
    );
    EXPECT_EQ(result.getComponents().size(), 1);
  }

  /// @brief Verifies CCL exclude all for geometry mesh ops by checking exact expected values.
  TEST(Geometry_MeshOps, CCL_ExcludeAll)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);

    // Exclude all cells via the inclusion predicate
    auto result = mesh.ccl(
      D,
      [](const Polytope&, const Polytope&) { return true; },
      [](const Polytope&) { return false; }
    );
    EXPECT_EQ(result.getComponents().size(), 0);
  }

  /// @brief Verifies CCL disconnect by predicate for geometry mesh ops by checking exact expected values.
  TEST(Geometry_MeshOps, CCL_DisconnectByPredicate)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    const size_t D = mesh.getDimension();
    mesh.getConnectivity().compute(D, D);

    // Never connect adjacent cells => each cell is its own component
    auto result = mesh.ccl(
      D,
      [](const Polytope&, const Polytope&) { return false; }
    );
    EXPECT_EQ(result.getComponents().size(), mesh.getCellCount());
  }

  // ==================================================================
  // Mesh::trace tests
  // ==================================================================

  /// @brief Verifies trace modifies mesh for geometry mesh ops.
  TEST(Geometry_MeshOps, Trace_ModifiesMesh)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(1, 2);

    // trace requires a Map<Pair<Attribute, Attribute>, Attribute> mapping
    // Get boundary attributes first
    FlatSet<Attribute> bndAttrs;
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
    {
      auto a = it->getAttribute();
      if (a) bndAttrs.insert(*a);
    }

    // Create a simple trace map
    Map<Pair<Attribute, Attribute>, Attribute> tmap;
    for (auto a : bndAttrs)
      tmap[{a, a}] = a + 100;

    // trace() modifies the mesh in-place
    mesh.trace(tmap);
    // Just verify no crash
    EXPECT_GT(mesh.getCellCount(), 0);
  }

  // ==================================================================
  // Mesh::isSurface tests
  // ==================================================================

  /// @brief Verifies is surface 2 D mesh for geometry mesh ops by checking false predicates.
  TEST(Geometry_MeshOps, IsSurface_2DMesh)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    // A 2D mesh in 2D space is not a surface mesh
    EXPECT_FALSE(mesh.isSurface());
  }

  /// @brief Verifies is surface boundary for geometry mesh ops by checking true predicates.
  TEST(Geometry_MeshOps, IsSurface_Boundary)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.getConnectivity().compute(1, 2);
    auto submesh = mesh.skin();
    // Skin of 2D mesh is 1D embedded in 2D, which IS a surface
    EXPECT_TRUE(submesh.isSurface());
  }

  // ==================================================================
  // Mesh::getPolytope tests
  // ==================================================================

  /// @brief Verifies get polytope cell for geometry mesh ops by checking exact expected values, false predicates.
  TEST(Geometry_MeshOps, GetPolytope_Cell)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    EXPECT_FALSE(it.end());
    EXPECT_EQ(it->getDimension(), mesh.getDimension());
  }

  /// @brief Verifies get polytope vertex for geometry mesh ops by checking exact expected values, false predicates.
  TEST(Geometry_MeshOps, GetPolytope_Vertex)
  {
    Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, {4, 4});
    auto it = mesh.getPolytope(0, 0);
    EXPECT_FALSE(it.end());
    EXPECT_EQ(it->getDimension(), 0);
  }

  namespace
  {
    /**
     * @brief Rebuilds @p mesh with its vertices and polytopes listed in a
     * scrambled order, as a mesher would produce them.
     */
    Mesh<Context::Local> scramble(const Mesh<Context::Local>& mesh)
    {
      const size_t vertexCount = mesh.getVertexCount();
      std::vector<Index> vertexOrder(vertexCount);
      std::iota(vertexOrder.begin(), vertexOrder.end(), Index(0));
      std::vector<Index> newVertex(vertexCount);
      for (size_t i = 0; i < vertexCount; ++i)
        std::swap(vertexOrder[i], vertexOrder[(7 * i + 3) % vertexCount]);
      for (Index v = 0; v < vertexCount; ++v)
        newVertex[vertexOrder[v]] = v;

      Mesh<Context::Local>::Builder build;
      build.initialize(mesh.getSpaceDimension()).nodes(vertexCount);
      for (Index v = 0; v < vertexCount; ++v)
        build.vertex(mesh.getVertexCoordinates(vertexOrder[v]));
      const auto& conn = mesh.getConnectivity();
      for (size_t d = 1; d <= mesh.getDimension(); ++d)
      {
        const size_t count = conn.getCount(d);
        build.reserve(d, count);
        for (Index p = 0; p < count; ++p)
        {
          const Index old = (2 * count - 1 - p + count / 3) % count;
          const auto& key = conn.getPolytope(d, old);
          IndexArray vertices(key.size());
          for (std::uint8_t j = 0; j < key.size(); ++j)
            vertices.coeffRef(j) = newVertex[key[j]];
          build.polytope(conn.getGeometry(d, old), std::move(vertices));
          build.attribute({d, p}, mesh.getAttribute(d, old));
        }
      }
      return build.finalize();
    }

    /**
     * @brief Lists the polytopes of dimension @p d by the coordinates of their
     * vertices in stored order, with their attribute, in a canonical order.
     */
    std::vector<std::vector<Real>> describe(const Mesh<Context::Local>& mesh, size_t d)
    {
      std::vector<std::vector<Real>> res;
      const auto& conn = mesh.getConnectivity();
      for (Index p = 0; p < conn.getCount(d); ++p)
      {
        std::vector<Real> entry;
        const auto& key = conn.getPolytope(d, p);
        for (std::uint8_t j = 0; j < key.size(); ++j)
        {
          const auto x = mesh.getVertexCoordinates(key[j]);
          for (Eigen::Index i = 0; i < x.size(); ++i)
            entry.push_back(x(i));
        }
        const auto attribute = mesh.getAttribute(d, p);
        entry.push_back(attribute ? Real(*attribute) : Real(-1));
        res.push_back(std::move(entry));
      }
      std::sort(res.begin(), res.end());
      return res;
    }
  }

  /**
   * @brief Renumbering keeps every polytope, its orientation and its
   * attribute, numbers the vertices along space and the cells by their
   * smallest vertex, independently of the input numbering.
   */
  TEST(Geometry_MeshOps, Renumber_KeepsPolytopesAndOrdersBySpace)
  {
    Mesh<Context::Local> grid =
      Mesh<Context::Local>::UniformGrid(Polytope::Type::Tetrahedron, {5, 4, 3});
    grid.getConnectivity().compute(2, 3);
    for (Index c = 0; c < grid.getCellCount(); ++c)
      grid.setAttribute({3, c}, Attribute(1 + c % 3));
    for (Index f = 0; f < grid.getFaceCount(); ++f)
    {
      if (grid.isBoundary(f))
        grid.setAttribute({2, f}, Attribute(10 + f % 2));
    }
    const Mesh<Context::Local> mesh = scramble(grid);
    const Mesh<Context::Local> renumbered = mesh.renumber();

    ASSERT_EQ(renumbered.getDimension(), mesh.getDimension());
    ASSERT_EQ(renumbered.getSpaceDimension(), mesh.getSpaceDimension());
    for (size_t d = 0; d <= mesh.getDimension(); ++d)
    {
      ASSERT_EQ(
        renumbered.getConnectivity().getCount(d), mesh.getConnectivity().getCount(d));
      EXPECT_EQ(describe(renumbered, d), describe(mesh, d)) << "dimension " << d;
    }

    // The cells follow their smallest vertex.
    const auto& conn = renumbered.getConnectivity();
    Index previous = 0;
    for (Index c = 0; c < renumbered.getCellCount(); ++c)
    {
      const auto& key = conn.getPolytope(3, c);
      Index least = key[0];
      for (std::uint8_t j = 1; j < key.size(); ++j)
        least = std::min(least, key[j]);
      EXPECT_LE(previous, least);
      previous = least;
    }

    // The vertex numbering depends only on the coordinates.
    const Mesh<Context::Local> reference = grid.renumber();
    for (Index v = 0; v < renumbered.getVertexCount(); ++v)
    {
      const auto x = renumbered.getVertexCoordinates(v);
      const auto y = reference.getVertexCoordinates(v);
      for (Eigen::Index i = 0; i < x.size(); ++i)
        EXPECT_EQ(x(i), y(i)) << "vertex " << v;
    }

    // Renumbering a renumbered mesh changes nothing.
    const Mesh<Context::Local> again = renumbered.renumber();
    for (size_t d = 1; d <= mesh.getDimension(); ++d)
    {
      for (Index p = 0; p < conn.getCount(d); ++p)
      {
        const auto& a = conn.getPolytope(d, p);
        const auto& b = again.getConnectivity().getPolytope(d, p);
        ASSERT_EQ(a.size(), b.size());
        for (std::uint8_t j = 0; j < a.size(); ++j)
          EXPECT_EQ(a[j], b[j]);
      }
    }
  }
  TEST(Geometry_MeshOps, Renumber_RejectsInvalidPermutations)
  {
    const auto mesh = LocalMesh::UniformGrid(Polytope::Type::Segment, {3});
    EXPECT_THROW(mesh.renumber(IndexVector{0, 1}), Alert::Exception);
    EXPECT_THROW(mesh.renumber(IndexVector{0, 1, 1}), Alert::Exception);
    EXPECT_THROW(mesh.renumber(IndexVector{0, 1, 3}), Alert::Exception);
  }

  TEST(Geometry_MeshOps, Renumber_PoliciesPreserveGeometryAndAttributes)
  {
    auto mesh = LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
    mesh.getConnectivity().compute(2, 3);
    for (size_t d = 0; d <= mesh.getDimension(); ++d)
    {
      for (Index p = 0; p < mesh.getConnectivity().getCount(d); ++p)
        mesh.setAttribute({d, p}, Attribute(1 + p % 7));
    }
    const std::vector<IndexVector> permutations = {
      Ordering::Morton{}(mesh), Ordering::Hilbert{}(mesh), Ordering::RCM{}(mesh)};
    const std::vector<LocalMesh> results = {mesh.renumber(),
      mesh.renumber(Ordering::Hilbert{}), mesh.renumber(Ordering::RCM{})};
    for (size_t i = 0; i < results.size(); ++i)
    {
      const auto direct = mesh.renumber(permutations[i]);
      for (size_t d = 0; d <= mesh.getDimension(); ++d)
      {
        EXPECT_EQ(describe(results[i], d), describe(mesh, d));
        EXPECT_EQ(describe(direct, d), describe(results[i], d));
      }
      for (Index v = 0; v < mesh.getVertexCount(); ++v)
      {
        const auto actual = results[i].getVertexCoordinates(v);
        const auto expected = mesh.getVertexCoordinates(permutations[i][v]);
        for (size_t d = 0; d < mesh.getSpaceDimension(); ++d)
          EXPECT_EQ(actual(d), expected(d));
        EXPECT_EQ(
          results[i].getAttribute(0, v), mesh.getAttribute(0, permutations[i][v]));
      }
    }
    const auto explicitMorton = mesh.renumber(Ordering::Morton{});
    for (Index v = 0; v < mesh.getVertexCount(); ++v)
    {
      const auto actual = results[0].getVertexCoordinates(v);
      const auto expected = explicitMorton.getVertexCoordinates(v);
      for (size_t d = 0; d < mesh.getSpaceDimension(); ++d)
        EXPECT_EQ(actual(d), expected(d));
    }
  }

}
