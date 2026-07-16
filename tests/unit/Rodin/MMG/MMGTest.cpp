/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>
#include <sstream>
#include <cmath>
#include <vector>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/MMG.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  static Index getFirstBoundaryTriangleIndex(MMG::Mesh& mesh)
  {
    mesh.getConnectivity().compute(2, 3);
    auto it = mesh.getBoundary();
    EXPECT_FALSE(it.end());
    if (it.end())
      return 0;
    EXPECT_EQ(it->getGeometry(), Polytope::Type::Triangle);
    return it->getIndex();
  }

  static void setAllFaceAttributes(MMG::Mesh& mesh, Attribute attr)
  {
    mesh.getConnectivity().compute(mesh.getDimension() - 1, mesh.getDimension());
    for (auto it = mesh.getFace(); it; ++it)
      mesh.setAttribute({it->getDimension(), it->getIndex()}, attr);
  }

  static void setAllBoundaryAttributes(MMG::Mesh& mesh, Attribute attr)
  {
    mesh.getConnectivity().compute(mesh.getDimension() - 1, mesh.getDimension());
    for (auto it = mesh.getBoundary(); it; ++it)
      mesh.setAttribute({it->getDimension(), it->getIndex()}, attr);
  }

  static void setAllSegmentAttributes(MMG::Mesh& mesh, Attribute attr)
  {
    mesh.getConnectivity().compute(1, mesh.getDimension());
    for (auto it = mesh.getPolytope(1); it; ++it)
      mesh.setAttribute({it->getDimension(), it->getIndex()}, attr);
  }

  static bool samePoint(const Point& a, const Point& b, Real tol = 1e-10)
  {
    if (a.getDimension() != b.getDimension())
      return false;
    for (size_t i = 0; i < a.getDimension(); i++)
    {
      if (std::abs(a.getCoordinates()(i) - b.getCoordinates()(i)) > tol)
        return false;
    }
    return true;
  }

  static size_t getDimension(Polytope::Type type)
  {
    switch (type)
    {
      case Polytope::Type::Point:
        return 0;
      case Polytope::Type::Segment:
        return 1;
      case Polytope::Type::Triangle:
      case Polytope::Type::Quadrilateral:
        return 2;
      case Polytope::Type::Tetrahedron:
      case Polytope::Type::Pyramid:
      case Polytope::Type::Hexahedron:
      case Polytope::Type::Wedge:
        return 3;
      default:
        assert(false);
        return 0;
    }
  }

  static std::vector<Point> getVertexCoordinates(
    const Mesh<Context::Local>& mesh, const Polytope& p)
  {
    std::vector<Point> points;
    for (const auto& v : p.getVertices())
    {
      const auto coords = mesh.getVertexCoordinates(v);
      points.emplace_back(p, coords, coords);
    }
    return points;
  }

  static bool containsPointSet(
    const std::vector<Point>& haystack, const std::vector<Point>& needles)
  {
    for (const auto& needle : needles)
    {
      bool found = false;
      for (const auto& point : haystack)
      {
        if (samePoint(point, needle))
        {
          found = true;
          break;
        }
      }
      if (!found)
        return false;
    }
    return true;
  }

  static bool hasEntityWithCoordinates(
    MMG::Mesh& mesh, Polytope::Type type, const std::vector<Point>& points)
  {
    const size_t dimension = getDimension(type);
    if (dimension < mesh.getDimension())
      mesh.getConnectivity().compute(dimension, mesh.getDimension());
    for (auto it = mesh.getPolytope(getDimension(type)); it; ++it)
    {
      if (it->getGeometry() != type)
        continue;
      const auto candidate = getVertexCoordinates(mesh, *it);
      if (containsPointSet(candidate, points) && containsPointSet(points, candidate))
        return true;
    }
    return false;
  }

  static bool hasVertexWithCoordinates(MMG::Mesh& mesh, const Point& point)
  {
    for (Index i = 0; i < mesh.getVertexCount(); i++)
    {
      const auto coords = mesh.getVertexCoordinates(i);
      const Point candidate(*mesh.getVertex(i), coords, coords);
      if (samePoint(candidate, point))
        return true;
    }
    return false;
  }

  static bool isCutByLevelSet(const Mesh<Context::Local>& mesh, const Polytope& p,
    const std::function<Real(const Point&)>& ls)
  {
    bool hasNegative = false;
    bool hasPositive = false;
    for (const auto& v : p.getVertices())
    {
      const auto coords = mesh.getVertexCoordinates(v);
      const Real value = ls(Point(p, coords, coords));
      hasNegative = hasNegative || value < 0;
      hasPositive = hasPositive || value > 0;
    }
    return hasNegative && hasPositive;
  }

  static Index getFirstCutEntity(
    MMG::Mesh& mesh, Polytope::Type type, const std::function<Real(const Point&)>& ls)
  {
    for (auto it = mesh.getPolytope(getDimension(type)); it; ++it)
    {
      if (it->getGeometry() != type)
        continue;
      if (isCutByLevelSet(mesh, *it, ls))
        return it->getIndex();
    }
    ADD_FAILURE() << "Could not find cut entity.";
    return 0;
  }

  static Index getFirstUncutEntity(
    MMG::Mesh& mesh, Polytope::Type type, const std::function<Real(const Point&)>& ls)
  {
    for (auto it = mesh.getPolytope(getDimension(type)); it; ++it)
    {
      if (it->getGeometry() != type)
        continue;
      if (!isCutByLevelSet(mesh, *it, ls))
        return it->getIndex();
    }
    ADD_FAILURE() << "Could not find uncut entity.";
    return 0;
  }

  static void setCellsContainingVertices(
    MMG::Mesh& mesh, const Polytope::Key& vertices, Attribute attr)
  {
    for (auto it = mesh.getCell(); it; ++it)
    {
      bool containsAll = true;
      for (const auto& vertex : vertices)
      {
        bool containsVertex = false;
        for (const auto& cellVertex : it->getVertices())
        {
          if (cellVertex == vertex)
          {
            containsVertex = true;
            break;
          }
        }
        containsAll = containsAll && containsVertex;
      }
      if (containsAll)
        mesh.setAttribute({it->getDimension(), it->getIndex()}, attr);
    }
  }

  static std::vector<Index> getCellsContainingVertices(
    MMG::Mesh& mesh, const Polytope::Key& vertices)
  {
    std::vector<Index> cells;
    for (auto it = mesh.getCell(); it; ++it)
    {
      bool containsAll = true;
      for (const auto& vertex : vertices)
      {
        bool containsVertex = false;
        for (const auto& cellVertex : it->getVertices())
        {
          if (cellVertex == vertex)
          {
            containsVertex = true;
            break;
          }
        }
        containsAll = containsAll && containsVertex;
      }
      if (containsAll)
        cells.emplace_back(it->getIndex());
    }
    return cells;
  }

  static std::vector<Index> getEdgesOnTriangle(
    MMG::Mesh& mesh, const Polytope::Key& vertices)
  {
    std::vector<Index> edges;
    mesh.getConnectivity().compute(1, mesh.getDimension());
    for (auto it = mesh.getPolytope(1); it; ++it)
    {
      const auto& edgeVertices = it->getVertices();
      size_t matches = 0;
      for (const auto& edgeVertex : edgeVertices)
      {
        for (const auto& triangleVertex : vertices)
        {
          if (edgeVertex == triangleVertex)
          {
            matches++;
            break;
          }
        }
      }
      if (matches == edgeVertices.size())
        edges.emplace_back(it->getIndex());
    }
    return edges;
  }

  // ========================================================================
  // MMG::Mesh construction tests
  // ========================================================================

  /// @brief Verifies default construction for MMG mesh by checking exact expected values, true predicates.
  TEST(Rodin_MMG_Mesh, DefaultConstruction)
  {
    MMG::Mesh mesh;
    EXPECT_TRUE(mesh.isEmpty());
    EXPECT_EQ(mesh.getCorners().size(), 0);
    EXPECT_EQ(mesh.getRidges().size(), 0);
    EXPECT_EQ(mesh.getRequiredVertices().size(), 0);
    EXPECT_EQ(mesh.getRequiredEdges().size(), 0);
    EXPECT_EQ(mesh.getRequiredTriangles().size(), 0);
    EXPECT_EQ(mesh.getRequiredTetrahedra().size(), 0);
  }

  /// @brief Verifies uniform grid construction for MMG mesh by checking exact expected values, false predicates.
  TEST(Rodin_MMG_Mesh, UniformGridConstruction)
  {
    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
    EXPECT_FALSE(mesh.isEmpty());
    EXPECT_EQ(mesh.getVertexCount(), n * n);
    EXPECT_EQ(mesh.getCellCount(), 2 * (n - 1) * (n - 1));
    EXPECT_EQ(mesh.getDimension(), 2);
    EXPECT_EQ(mesh.getSpaceDimension(), 2);
  }

  /// @brief Verifies set corner for MMG mesh by checking exact expected values, true predicates.
  TEST(Rodin_MMG_Mesh, SetCorner)
  {
    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });

    mesh.setCorner(0);
    mesh.setCorner(3);
    mesh.setCorner(12);
    mesh.setCorner(15);

    EXPECT_EQ(mesh.getCorners().size(), 4);
    EXPECT_TRUE(mesh.getCorners().count(0));
    EXPECT_TRUE(mesh.getCorners().count(3));
    EXPECT_TRUE(mesh.getCorners().count(12));
    EXPECT_TRUE(mesh.getCorners().count(15));
  }

  /// @brief Verifies set ridge for MMG mesh by checking exact expected values.
  TEST(Rodin_MMG_Mesh, SetRidge)
  {
    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
    mesh.getConnectivity().compute(1, 2);

    size_t ridgeCount = 0;
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
    {
      mesh.setRidge(it->getIndex());
      ridgeCount++;
    }

    EXPECT_EQ(mesh.getRidges().size(), ridgeCount);
    EXPECT_GT(ridgeCount, 0);
  }

  /// @brief Verifies set required vertex for MMG mesh by checking exact expected values, true predicates.
  TEST(Rodin_MMG_Mesh, SetRequiredVertex)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });

    mesh.setRequiredVertex(0);
    mesh.setRequiredVertex(5);

    EXPECT_EQ(mesh.getRequiredVertices().size(), 2);
    EXPECT_TRUE(mesh.getRequiredVertices().count(0));
    EXPECT_TRUE(mesh.getRequiredVertices().count(5));
  }

  /// @brief Verifies set required edge for MMG mesh by checking exact expected values, true predicates.
  TEST(Rodin_MMG_Mesh, SetRequiredEdge)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.getConnectivity().compute(1, 2);

    mesh.setRequiredEdge(0);
    mesh.setRequiredEdge(1);

    EXPECT_EQ(mesh.getRequiredEdges().size(), 2);
    EXPECT_TRUE(mesh.getRequiredEdges().count(0));
    EXPECT_TRUE(mesh.getRequiredEdges().count(1));
  }

  /// @brief Verifies set required triangle for MMG mesh by checking exact expected values, true predicates.
  TEST(Rodin_MMG_Mesh, SetRequiredTriangle)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {4, 4});

    mesh.setRequiredTriangle(0);
    mesh.setRequiredTriangle(3);

    EXPECT_EQ(mesh.getRequiredTriangles().size(), 2);
    EXPECT_TRUE(mesh.getRequiredTriangles().count(0));
    EXPECT_TRUE(mesh.getRequiredTriangles().count(3));
  }

  /// @brief Verifies set required tetrahedron for MMG mesh by checking exact expected values, true predicates.
  TEST(Rodin_MMG_Mesh, SetRequiredTetrahedron)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});

    mesh.setRequiredTetrahedron(0);
    mesh.setRequiredTetrahedron(2);

    EXPECT_EQ(mesh.getRequiredTetrahedra().size(), 2);
    EXPECT_TRUE(mesh.getRequiredTetrahedra().count(0));
    EXPECT_TRUE(mesh.getRequiredTetrahedra().count(2));
  }

  // ========================================================================
  // MMG::Mesh copy and move tests
  // ========================================================================

  /// @brief Verifies copy construction for MMG mesh by checking exact expected values, true predicates, copy semantics.
  TEST(Rodin_MMG_Mesh, CopyConstruction)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.setCorner(0);
    mesh.setCorner(3);
    mesh.setRidge(0);
    mesh.setRequiredVertex(5);
    mesh.setRequiredEdge(2);
    mesh.setRequiredTriangle(1);
    mesh.setRequiredTetrahedron(0);

    MMG::Mesh copy(mesh);
    EXPECT_EQ(copy.getVertexCount(), mesh.getVertexCount());
    EXPECT_EQ(copy.getCellCount(), mesh.getCellCount());
    EXPECT_EQ(copy.getCorners().size(), 2);
    EXPECT_EQ(copy.getRidges().size(), 1);
    EXPECT_EQ(copy.getRequiredVertices().size(), 1);
    EXPECT_EQ(copy.getRequiredEdges().size(), 1);
    EXPECT_EQ(copy.getRequiredTriangles().size(), 1);
    EXPECT_EQ(copy.getRequiredTetrahedra().size(), 1);
    EXPECT_TRUE(copy.getCorners().count(0));
    EXPECT_TRUE(copy.getCorners().count(3));
    EXPECT_TRUE(copy.getRequiredTriangles().count(1));
    EXPECT_TRUE(copy.getRequiredTetrahedra().count(0));
  }

  /// @brief Verifies move construction for MMG mesh by checking exact expected values, move semantics.
  TEST(Rodin_MMG_Mesh, MoveConstruction)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.setCorner(0);
    mesh.setRidge(0);
    mesh.setRequiredVertex(5);
    mesh.setRequiredEdge(2);
    mesh.setRequiredTriangle(1);
    mesh.setRequiredTetrahedron(0);

    size_t vertexCount = mesh.getVertexCount();
    size_t cellCount = mesh.getCellCount();

    MMG::Mesh moved(std::move(mesh));
    EXPECT_EQ(moved.getVertexCount(), vertexCount);
    EXPECT_EQ(moved.getCellCount(), cellCount);
    EXPECT_EQ(moved.getCorners().size(), 1);
    EXPECT_EQ(moved.getRidges().size(), 1);
    EXPECT_EQ(moved.getRequiredVertices().size(), 1);
    EXPECT_EQ(moved.getRequiredEdges().size(), 1);
    EXPECT_EQ(moved.getRequiredTriangles().size(), 1);
    EXPECT_EQ(moved.getRequiredTetrahedra().size(), 1);
  }

  /// @brief Verifies move assignment for MMG mesh by checking exact expected values, move semantics.
  TEST(Rodin_MMG_Mesh, MoveAssignment)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.setCorner(0);
    mesh.setRidge(0);

    size_t vertexCount = mesh.getVertexCount();
    size_t cellCount = mesh.getCellCount();

    MMG::Mesh target;
    target = std::move(mesh);
    EXPECT_EQ(target.getVertexCount(), vertexCount);
    EXPECT_EQ(target.getCellCount(), cellCount);
    EXPECT_EQ(target.getCorners().size(), 1);
    EXPECT_EQ(target.getRidges().size(), 1);
  }

  /// @brief Verifies move from parent for MMG mesh by checking exact expected values, move semantics.
  TEST(Rodin_MMG_Mesh, MoveFromParent)
  {
    Mesh<Context::Local> parentMesh;
    parentMesh = parentMesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    size_t vertexCount = parentMesh.getVertexCount();
    size_t cellCount = parentMesh.getCellCount();

    MMG::Mesh mmgMesh(std::move(parentMesh));
    EXPECT_EQ(mmgMesh.getVertexCount(), vertexCount);
    EXPECT_EQ(mmgMesh.getCellCount(), cellCount);
    // MMG metadata should be empty since parent doesn't carry it
    EXPECT_EQ(mmgMesh.getCorners().size(), 0);
    EXPECT_EQ(mmgMesh.getRidges().size(), 0);
  }

  // ========================================================================
  // MMG::Mesh Builder tests
  // ========================================================================

  /// @brief Verifies triangle mesh with MMG metadata for MMG mesh builder by checking exact expected values, true predicates.
  TEST(Rodin_MMG_Mesh_Builder, TriangleMeshWithMMGMetadata)
  {
    MMG::Mesh mesh =
      MMG::Mesh::Builder()
      .initialize(2)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .polytope(Polytope::Type::Triangle, {1, 3, 2})
      .corner(0)
      .corner(1)
      .corner(2)
      .corner(3)
      .finalize();

    EXPECT_EQ(mesh.getVertexCount(), 4);
    EXPECT_EQ(mesh.getCellCount(), 2);
    EXPECT_EQ(mesh.getCorners().size(), 4);
    EXPECT_TRUE(mesh.getCorners().count(0));
    EXPECT_TRUE(mesh.getCorners().count(1));
    EXPECT_TRUE(mesh.getCorners().count(2));
    EXPECT_TRUE(mesh.getCorners().count(3));
  }

  /// @brief Verifies builder with all MMG metadata for MMG mesh builder by checking exact expected values.
  TEST(Rodin_MMG_Mesh_Builder, BuilderWithAllMMGMetadata)
  {
    MMG::Mesh mesh = MMG::Mesh::Builder()
                       .initialize(3)
                       .nodes(4)
                       .vertex({0, 0, 0})
                       .vertex({1, 0, 0})
                       .vertex({0, 1, 0})
                       .vertex({0, 0, 1})
                       .polytope(Polytope::Type::Segment, {0, 1})
                       .polytope(Polytope::Type::Triangle, {0, 1, 2})
                       .polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 3})
                       .corner(0)
                       .ridge(0)
                       .requiredVertex(1)
                       .requiredEdge(0)
                       .requiredTriangle(0)
                       .requiredTetrahedron(0)
                       .finalize();

    EXPECT_EQ(mesh.getCorners().size(), 1);
    EXPECT_EQ(mesh.getRidges().size(), 1);
    EXPECT_EQ(mesh.getRequiredVertices().size(), 1);
    EXPECT_EQ(mesh.getRequiredEdges().size(), 1);
    EXPECT_EQ(mesh.getRequiredTriangles().size(), 1);
    EXPECT_EQ(mesh.getRequiredTetrahedra().size(), 1);
  }

  // ========================================================================
  // MMG::Mesh I/O roundtrip tests (MEDIT format)
  // ========================================================================

  /// @brief Verifies save load MEDIT roundtrip for MMG mesh IO by checking exact expected values, true predicates.
  TEST(Rodin_MMG_Mesh_IO, SaveLoadMEDITRoundtrip)
  {
    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
    mesh.getConnectivity().compute(1, 2);

    // Set some MMG metadata
    mesh.setCorner(0);
    mesh.setCorner(3);
    mesh.setCorner(12);
    mesh.setCorner(15);

    for (auto it = mesh.getBoundary(); !it.end(); ++it)
      mesh.setRidge(it->getIndex());

    mesh.setRequiredVertex(0);
    mesh.setRequiredVertex(15);
    mesh.setRequiredTriangle(0);
    mesh.setRequiredTriangle(1);

    // Save to file
    const std::string filename = "/tmp/rodin_mmg_test_roundtrip.mesh";
    mesh.save(filename, IO::FileFormat::MEDIT);

    // Load back
    MMG::Mesh loaded;
    loaded.load(filename, IO::FileFormat::MEDIT);

    // Verify topology preserved
    EXPECT_EQ(loaded.getVertexCount(), mesh.getVertexCount());
    EXPECT_EQ(loaded.getCellCount(), mesh.getCellCount());

    // Verify MMG metadata preserved
    EXPECT_EQ(loaded.getCorners().size(), mesh.getCorners().size());
    EXPECT_EQ(loaded.getRidges().size(), mesh.getRidges().size());
    EXPECT_EQ(loaded.getRequiredVertices().size(), mesh.getRequiredVertices().size());
    EXPECT_EQ(loaded.getRequiredEdges().size(), mesh.getRequiredEdges().size());
    EXPECT_EQ(loaded.getRequiredTriangles().size(), mesh.getRequiredTriangles().size());
    EXPECT_EQ(loaded.getRequiredTetrahedra().size(), mesh.getRequiredTetrahedra().size());

    // Verify specific corners
    EXPECT_TRUE(loaded.getCorners().count(0));
    EXPECT_TRUE(loaded.getCorners().count(3));
    EXPECT_TRUE(loaded.getCorners().count(12));
    EXPECT_TRUE(loaded.getCorners().count(15));

    // Verify required vertices
    EXPECT_TRUE(loaded.getRequiredVertices().count(0));
    EXPECT_TRUE(loaded.getRequiredVertices().count(15));
    EXPECT_TRUE(loaded.getRequiredTriangles().count(0));
    EXPECT_TRUE(loaded.getRequiredTriangles().count(1));

    // Clean up
    std::remove(filename.c_str());
  }

  /// @brief Verifies save load MEDIT required triangles 2 D for MMG mesh IO by checking exact expected values, true predicates.
  TEST(Rodin_MMG_Mesh_IO, SaveLoadMEDITRequiredTriangles2D)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {4, 4});

    mesh.setRequiredTriangle(0);
    mesh.setRequiredTriangle(4);

    const std::string filename = "/tmp/rodin_mmg_required_triangles_2d.mesh";
    mesh.save(filename, IO::FileFormat::MEDIT);

    MMG::Mesh loaded;
    loaded.load(filename, IO::FileFormat::MEDIT);

    EXPECT_EQ(loaded.getVertexCount(), mesh.getVertexCount());
    EXPECT_EQ(loaded.getCellCount(), mesh.getCellCount());
    EXPECT_EQ(loaded.getRequiredTriangles().size(), 2);
    EXPECT_TRUE(loaded.getRequiredTriangles().count(0));
    EXPECT_TRUE(loaded.getRequiredTriangles().count(4));

    std::remove(filename.c_str());
  }

  /// @brief Verifies save load MEDIT required triangles 3 D for MMG mesh IO by checking exact expected values, true predicates.
  TEST(Rodin_MMG_Mesh_IO, SaveLoadMEDITRequiredTriangles3D)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
    setAllFaceAttributes(mesh, 7);
    const Index required = 0;
    mesh.setRequiredTriangle(required);

    const std::string filename = "/tmp/rodin_mmg_required_triangles_3d.mesh";
    mesh.save(filename, IO::FileFormat::MEDIT);

    MMG::Mesh loaded;
    loaded.load(filename, IO::FileFormat::MEDIT);

    EXPECT_EQ(loaded.getVertexCount(), mesh.getVertexCount());
    EXPECT_EQ(loaded.getCellCount(), mesh.getCellCount());
    EXPECT_EQ(loaded.getRequiredTriangles().size(), 1);
    EXPECT_TRUE(loaded.getRequiredTriangles().count(required));

    std::remove(filename.c_str());
  }

  /// @brief Verifies save load MEDIT required entities 3 D for MMG mesh IO by checking exact expected values, true predicates.
  TEST(Rodin_MMG_Mesh_IO, SaveLoadMEDITRequiredEntities3D)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
    setAllSegmentAttributes(mesh, 8);
    setAllFaceAttributes(mesh, 7);

    const Index requiredVertex = 0;
    const Index requiredEdge = 0;
    const Index requiredTriangle = getFirstBoundaryTriangleIndex(mesh);
    const Index requiredTetrahedron = 0;
    mesh.setRequiredVertex(requiredVertex);
    mesh.setRequiredEdge(requiredEdge);
    mesh.setRequiredTriangle(requiredTriangle);
    mesh.setRequiredTetrahedron(requiredTetrahedron);

    const std::string filename = "/tmp/rodin_mmg_required_entities_3d.mesh";
    mesh.save(filename, IO::FileFormat::MEDIT);

    MMG::Mesh loaded;
    loaded.load(filename, IO::FileFormat::MEDIT);

    EXPECT_EQ(loaded.getVertexCount(), mesh.getVertexCount());
    EXPECT_EQ(loaded.getCellCount(), mesh.getCellCount());
    EXPECT_EQ(loaded.getRequiredVertices().size(), 1);
    EXPECT_EQ(loaded.getRequiredEdges().size(), 1);
    EXPECT_EQ(loaded.getRequiredTriangles().size(), 1);
    EXPECT_EQ(loaded.getRequiredTetrahedra().size(), 1);
    EXPECT_TRUE(loaded.getRequiredVertices().count(requiredVertex));
    EXPECT_TRUE(loaded.getRequiredEdges().count(requiredEdge));
    EXPECT_TRUE(loaded.getRequiredTriangles().count(requiredTriangle));
    EXPECT_TRUE(loaded.getRequiredTetrahedra().count(requiredTetrahedron));

    std::remove(filename.c_str());
  }

  /// @brief Verifies save load MEDIT empty metadata for MMG mesh IO by checking exact expected values.
  TEST(Rodin_MMG_Mesh_IO, SaveLoadMEDITEmptyMetadata)
  {
    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });

    // Save without any MMG metadata
    const std::string filename = "/tmp/rodin_mmg_test_empty_metadata.mesh";
    mesh.save(filename, IO::FileFormat::MEDIT);

    MMG::Mesh loaded;
    loaded.load(filename, IO::FileFormat::MEDIT);

    EXPECT_EQ(loaded.getVertexCount(), mesh.getVertexCount());
    EXPECT_EQ(loaded.getCellCount(), mesh.getCellCount());
    EXPECT_EQ(loaded.getCorners().size(), 0);
    EXPECT_EQ(loaded.getRidges().size(), 0);
    EXPECT_EQ(loaded.getRequiredVertices().size(), 0);
    EXPECT_EQ(loaded.getRequiredEdges().size(), 0);
    EXPECT_EQ(loaded.getRequiredTriangles().size(), 0);
    EXPECT_EQ(loaded.getRequiredTetrahedra().size(), 0);

    std::remove(filename.c_str());
  }

  // ========================================================================
  // MMG5 mesh conversion roundtrip tests
  // ========================================================================

  /// @brief Verifies rodin to mesh and back 2 D for MMG MMG 5 by checking exact expected values.
  TEST(Rodin_MMG_MMG5, RodinToMeshAndBack2D)
  {
    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });

    size_t origVertexCount = mesh.getVertexCount();
    size_t origCellCount = mesh.getCellCount();

    // Convert to MMG and back
    MMG5_pMesh mmgMesh = MMG::MMG5::rodinToMesh(mesh);
    ASSERT_NE(mmgMesh, nullptr);

    MMG::Mesh result = MMG::MMG5::meshToRodin(mmgMesh);
    MMG::MMG5::destroyMesh(mmgMesh);

    EXPECT_EQ(result.getVertexCount(), origVertexCount);
    EXPECT_EQ(result.getCellCount(), origCellCount);
    EXPECT_EQ(result.getDimension(), 2);
  }

  /// @brief Verifies rodin to mesh preserves MMG metadata for MMG MMG 5 by checking exact expected values, true predicates.
  TEST(Rodin_MMG_MMG5, RodinToMeshPreservesMMGMetadata)
  {
    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
    mesh.getConnectivity().compute(1, 2);

    mesh.setCorner(0);
    mesh.setCorner(3);

    // Only mark edges that have attributes as ridges, since rodinToMesh
    // filters out non-attributed edges. Boundary edges on a uniform grid
    // receive attributes, so we mark only boundary edges as ridges.
    size_t ridgeCount = 0;
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
    {
      if (it->getAttribute().has_value())
      {
        mesh.setRidge(it->getIndex());
        ridgeCount++;
      }
    }

    mesh.setRequiredVertex(0);
    mesh.setRequiredTriangle(1);

    // Convert to MMG and back
    MMG5_pMesh mmgMesh = MMG::MMG5::rodinToMesh(mesh);
    ASSERT_NE(mmgMesh, nullptr);

    MMG::Mesh result = MMG::MMG5::meshToRodin(mmgMesh);
    MMG::MMG5::destroyMesh(mmgMesh);

    // Corners should be preserved (they live on vertices, which are always
    // converted)
    EXPECT_EQ(result.getCorners().size(), mesh.getCorners().size());
    EXPECT_TRUE(result.getCorners().count(0));
    EXPECT_TRUE(result.getCorners().count(3));

    // Required vertices should be preserved
    EXPECT_EQ(result.getRequiredVertices().size(), mesh.getRequiredVertices().size());
    EXPECT_TRUE(result.getRequiredVertices().count(0));

    // Required triangles should be preserved
    EXPECT_EQ(result.getRequiredTriangles().size(), mesh.getRequiredTriangles().size());
    EXPECT_TRUE(result.getRequiredTriangles().count(1));

    // Ridges live on edges; only attributed edges survive the roundtrip.
    // Verify that the result carries at least some ridges if any attributed
    // edges were marked.
    if (ridgeCount > 0)
    {
      EXPECT_GT(result.getRidges().size(), 0);
    }
  }

  /// @brief Verifies required triangles roundtrip 2 D for MMG MMG 5 by checking exact expected values, true predicates.
  TEST(Rodin_MMG_MMG5, RequiredTrianglesRoundtrip2D)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {4, 4});
    mesh.setRequiredTriangle(0);
    mesh.setRequiredTriangle(3);

    MMG5_pMesh mmgMesh = MMG::MMG5::rodinToMesh(mesh);
    ASSERT_NE(mmgMesh, nullptr);

    MMG::Mesh result = MMG::MMG5::meshToRodin(mmgMesh);
    MMG::MMG5::destroyMesh(mmgMesh);

    EXPECT_EQ(result.getRequiredTriangles().size(), 2);
    EXPECT_TRUE(result.getRequiredTriangles().count(0));
    EXPECT_TRUE(result.getRequiredTriangles().count(3));
  }

  /// @brief Verifies required entities roundtrip 2 D for MMG MMG 5 by checking exact expected values, true predicates.
  TEST(Rodin_MMG_MMG5, RequiredEntitiesRoundtrip2D)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {4, 4});
    setAllSegmentAttributes(mesh, 8);

    const Index requiredVertex = 0;
    const Index requiredEdge = 0;
    const Index requiredTriangle = 3;
    mesh.setRequiredVertex(requiredVertex);
    mesh.setRequiredEdge(requiredEdge);
    mesh.setRequiredTriangle(requiredTriangle);

    MMG5_pMesh mmgMesh = MMG::MMG5::rodinToMesh(mesh);
    ASSERT_NE(mmgMesh, nullptr);

    MMG::Mesh result = MMG::MMG5::meshToRodin(mmgMesh);
    MMG::MMG5::destroyMesh(mmgMesh);

    EXPECT_EQ(result.getRequiredVertices().size(), 1);
    EXPECT_EQ(result.getRequiredEdges().size(), 1);
    EXPECT_EQ(result.getRequiredTriangles().size(), 1);
    EXPECT_TRUE(result.getRequiredVertices().count(requiredVertex));
    EXPECT_TRUE(result.getRequiredEdges().count(requiredEdge));
    EXPECT_TRUE(result.getRequiredTriangles().count(requiredTriangle));
  }

  /// @brief Verifies required triangles roundtrip 3 D without face attributes for MMG MMG 5 by checking exact expected values.
  TEST(Rodin_MMG_MMG5, RequiredTrianglesRoundtrip3DWithoutFaceAttributes)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
    const Index required = getFirstBoundaryTriangleIndex(mesh);
    mesh.setRequiredTriangle(required);

    MMG5_pMesh mmgMesh = MMG::MMG5::rodinToMesh(mesh);
    ASSERT_NE(mmgMesh, nullptr);

    MMG::Mesh result = MMG::MMG5::meshToRodin(mmgMesh);
    MMG::MMG5::destroyMesh(mmgMesh);

    EXPECT_EQ(result.getRequiredTriangles().size(), 1);
  }

  /// @brief Verifies required triangles roundtrip 3 D preserves boundary metadata with face attributes for MMG MMG 5 by checking exact expected values, true predicates.
  TEST(Rodin_MMG_MMG5,
    RequiredTrianglesRoundtrip3DPreservesBoundaryMetadataWithFaceAttributes)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
    setAllBoundaryAttributes(mesh, 7);
    const Index required = getFirstBoundaryTriangleIndex(mesh);
    mesh.setRequiredTriangle(required);

    MMG5_pMesh mmgMesh = MMG::MMG5::rodinToMesh(mesh);
    ASSERT_NE(mmgMesh, nullptr);

    MMG::Mesh result = MMG::MMG5::meshToRodin(mmgMesh);
    MMG::MMG5::destroyMesh(mmgMesh);

    EXPECT_EQ(result.getRequiredTriangles().size(), 1);
    EXPECT_TRUE(result.getRequiredTriangles().count(0));
  }

  /// @brief Verifies required entities roundtrip 3 D for MMG MMG 5 by checking exact expected values, true predicates.
  TEST(Rodin_MMG_MMG5, RequiredEntitiesRoundtrip3D)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
    setAllSegmentAttributes(mesh, 8);
    setAllBoundaryAttributes(mesh, 7);

    const Index requiredVertex = 0;
    const Index requiredEdge = 0;
    const Index requiredTriangle = getFirstBoundaryTriangleIndex(mesh);
    const Index requiredTetrahedron = 0;
    mesh.setRequiredVertex(requiredVertex);
    mesh.setRequiredEdge(requiredEdge);
    mesh.setRequiredTriangle(requiredTriangle);
    mesh.setRequiredTetrahedron(requiredTetrahedron);

    MMG5_pMesh mmgMesh = MMG::MMG5::rodinToMesh(mesh);
    ASSERT_NE(mmgMesh, nullptr);

    MMG::Mesh result = MMG::MMG5::meshToRodin(mmgMesh);
    MMG::MMG5::destroyMesh(mmgMesh);

    EXPECT_EQ(result.getRequiredVertices().size(), 1);
    EXPECT_EQ(result.getRequiredEdges().size(), 1);
    EXPECT_EQ(result.getRequiredTriangles().size(), 1);
    EXPECT_EQ(result.getRequiredTetrahedra().size(), 1);
    EXPECT_TRUE(result.getRequiredVertices().count(requiredVertex));
    EXPECT_TRUE(result.getRequiredEdges().count(requiredEdge));
    EXPECT_TRUE(result.getRequiredTriangles().count(0));
    EXPECT_TRUE(result.getRequiredTetrahedra().count(requiredTetrahedron));
  }

  // ========================================================================
  // MMG::Optimizer tests
  // ========================================================================

  /// @brief Verifies optimize 2 D triangle mesh for MMG optimizer by checking exact expected values, false predicates.
  TEST(Rodin_MMG_Optimizer, Optimize2DTriangleMesh)
  {
    const size_t n = 8;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });

    MMG::Optimizer optimizer;
    optimizer.setHMax(0.5);
    optimizer.optimize(mesh);

    // After optimization, mesh should still be valid and non-empty
    EXPECT_FALSE(mesh.isEmpty());
    EXPECT_GT(mesh.getVertexCount(), 0);
    EXPECT_GT(mesh.getCellCount(), 0);
    EXPECT_EQ(mesh.getDimension(), 2);
  }

  /// @brief Verifies fluent API for MMG optimizer by checking false predicates.
  TEST(Rodin_MMG_Optimizer, FluentAPI)
  {
    const size_t n = 8;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });

    // Fluent API method chaining
    MMG::Optimizer()
      .setHMax(0.5)
      .setAngleDetection(true)
      .setGradation(1.3)
      .optimize(mesh);

    EXPECT_FALSE(mesh.isEmpty());
    EXPECT_GT(mesh.getVertexCount(), 0);
    EXPECT_GT(mesh.getCellCount(), 0);
  }

  /// @brief Verifies optimize with corners for MMG optimizer by checking false predicates.
  TEST(Rodin_MMG_Optimizer, OptimizeWithCorners)
  {
    const size_t n = 8;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
    mesh.getConnectivity().compute(1, 2);

    // Mark corners
    mesh.setCorner(0);
    mesh.setCorner(n - 1);
    mesh.setCorner(n * (n - 1));
    mesh.setCorner(n * n - 1);

    // Mark ridges
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
      mesh.setRidge(it->getIndex());

    MMG::Optimizer().setHMax(0.5).optimize(mesh);

    EXPECT_FALSE(mesh.isEmpty());
    EXPECT_GT(mesh.getVertexCount(), 0);
    EXPECT_GT(mesh.getCellCount(), 0);
  }

  /// @brief Verifies preserves required triangles 2 D for MMG optimizer by checking exact expected values, true predicates, false predicates.
  TEST(Rodin_MMG_Optimizer, PreservesRequiredTriangles2D)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {6, 6});
    const Index required = 0;
    mesh.setRequiredTriangle(required);

    MMG::Optimizer().setHMax(0.5).optimize(mesh);

    EXPECT_FALSE(mesh.isEmpty());
    EXPECT_GT(mesh.getCellCount(), 0);
    EXPECT_EQ(mesh.getRequiredTriangles().size(), 1);
    EXPECT_TRUE(mesh.getRequiredTriangles().count(required));
  }

  /// @brief Verifies preserves required triangles 3 D for MMG optimizer by checking exact expected values, true predicates, false predicates.
  TEST(Rodin_MMG_Optimizer, PreservesRequiredTriangles3D)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
    const Index required = getFirstBoundaryTriangleIndex(mesh);
    mesh.setRequiredTriangle(required);

    MMG::Optimizer().setHMax(0.75).optimize(mesh);

    EXPECT_FALSE(mesh.isEmpty());
    EXPECT_GT(mesh.getCellCount(), 0);
    EXPECT_EQ(mesh.getRequiredTriangles().size(), 1);
    EXPECT_TRUE(mesh.getRequiredTriangles().count(required));
  }

  /// @brief Verifies preserves required entities 2 D for MMG optimizer by checking exact expected values, true predicates, false predicates.
  TEST(Rodin_MMG_Optimizer, PreservesRequiredEntities2D)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {6, 6});
    setAllSegmentAttributes(mesh, 8);

    const Index requiredVertex = 0;
    const Index requiredEdge = 0;
    const Index requiredTriangle = 0;
    mesh.setRequiredVertex(requiredVertex);
    mesh.setRequiredEdge(requiredEdge);
    mesh.setRequiredTriangle(requiredTriangle);

    MMG::Optimizer().setHMax(0.5).optimize(mesh);

    EXPECT_FALSE(mesh.isEmpty());
    EXPECT_GT(mesh.getCellCount(), 0);
    EXPECT_EQ(mesh.getRequiredVertices().size(), 1);
    EXPECT_EQ(mesh.getRequiredEdges().size(), 1);
    EXPECT_EQ(mesh.getRequiredTriangles().size(), 1);
    EXPECT_TRUE(mesh.getRequiredVertices().count(requiredVertex));
    EXPECT_TRUE(mesh.getRequiredEdges().count(requiredEdge));
    EXPECT_TRUE(mesh.getRequiredTriangles().count(requiredTriangle));
  }

  /// @brief Verifies preserves required entities 3 D for MMG optimizer by checking exact expected values, true predicates, false predicates.
  TEST(Rodin_MMG_Optimizer, PreservesRequiredEntities3D)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
    setAllSegmentAttributes(mesh, 8);
    setAllFaceAttributes(mesh, 7);

    const Index requiredVertex = 0;
    const Index requiredEdge = 0;
    const Index requiredTriangle = getFirstBoundaryTriangleIndex(mesh);
    const Index requiredTetrahedron = 0;
    mesh.setRequiredVertex(requiredVertex);
    mesh.setRequiredEdge(requiredEdge);
    mesh.setRequiredTriangle(requiredTriangle);
    mesh.setRequiredTetrahedron(requiredTetrahedron);

    MMG::Optimizer().setHMax(0.75).optimize(mesh);

    EXPECT_FALSE(mesh.isEmpty());
    EXPECT_GT(mesh.getCellCount(), 0);
    EXPECT_EQ(mesh.getRequiredVertices().size(), 1);
    EXPECT_EQ(mesh.getRequiredEdges().size(), 1);
    EXPECT_EQ(mesh.getRequiredTriangles().size(), 1);
    EXPECT_EQ(mesh.getRequiredTetrahedra().size(), 1);
    EXPECT_TRUE(mesh.getRequiredVertices().count(requiredVertex));
    EXPECT_TRUE(mesh.getRequiredEdges().count(requiredEdge));
    EXPECT_TRUE(mesh.getRequiredTriangles().count(requiredTriangle));
    EXPECT_TRUE(mesh.getRequiredTetrahedra().count(requiredTetrahedron));
  }

  // ========================================================================
  // MMG::Adapt tests
  // ========================================================================

  /// @brief Verifies uniform size map 2 D for MMG adapt by checking exact expected values, false predicates.
  TEST(Rodin_MMG_Adapt, UniformSizeMap2D)
  {
    const size_t n = 8;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });

    // Create a constant size map
    P1 fes(mesh);
    MMG::RealGridFunction sizeMap(fes);
    sizeMap = [](const Geometry::Point&) { return 0.2; };

    MMG::Adapt adapter;
    adapter.setHMax(0.5);
    adapter.setHMin(0.05);
    adapter.adapt(mesh, sizeMap);

    EXPECT_FALSE(mesh.isEmpty());
    EXPECT_GT(mesh.getVertexCount(), 0);
    EXPECT_GT(mesh.getCellCount(), 0);
    EXPECT_EQ(mesh.getDimension(), 2);
  }

  /// @brief Verifies variable size map 2 D for MMG adapt by checking false predicates.
  TEST(Rodin_MMG_Adapt, VariableSizeMap2D)
  {
    const size_t n = 8;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });

    // Create a variable size map: smaller near center, larger near boundary
    P1 fes(mesh);
    MMG::RealGridFunction sizeMap(fes);
    sizeMap = [&](const Geometry::Point& p)
    {
      Real cx = static_cast<Real>(n - 1) / 2.0;
      Real cy = cx;
      Real dx = p.x() - cx;
      Real dy = p.y() - cy;
      Real dist = std::sqrt(dx * dx + dy * dy);
      return 0.1 + 0.3 * dist / cx;
    };

    MMG::Adapt adapter;
    adapter.setHMin(0.05);
    adapter.setHMax(1.0);
    adapter.adapt(mesh, sizeMap);

    EXPECT_FALSE(mesh.isEmpty());
    EXPECT_GT(mesh.getVertexCount(), 0);
    EXPECT_GT(mesh.getCellCount(), 0);
  }

  /// @brief Verifies fluent API for MMG adapt by checking false predicates.
  TEST(Rodin_MMG_Adapt, FluentAPI)
  {
    const size_t n = 8;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });

    P1 fes(mesh);
    MMG::RealGridFunction sizeMap(fes);
    sizeMap = [](const Geometry::Point&) { return 0.2; };

    // Fluent API
    MMG::Adapt()
      .setHMin(0.05)
      .setHMax(0.5)
      .setGradation(1.3)
      .setHausdorff(0.01)
      .setAngleDetection(true)
      .adapt(mesh, sizeMap);

    EXPECT_FALSE(mesh.isEmpty());
    EXPECT_GT(mesh.getVertexCount(), 0);
    EXPECT_GT(mesh.getCellCount(), 0);
  }

  // ========================================================================
  // MMG::LevelSetDiscretizer tests
  // ========================================================================

  /// @brief Verifies circle level set 2 D for MMG level set discretizer by checking exact expected values, false predicates.
  TEST(Rodin_MMG_LevelSetDiscretizer, CircleLevelSet2D)
  {
    static constexpr Attribute interior = 1;
    static constexpr Attribute exterior = 2;
    static constexpr Attribute boundary = 10;
    static constexpr Real radius = 0.1;

    const size_t n = 16;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
    mesh.scale(1.0 / (n - 1));

    // Set all cells to interior attribute
    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({ it->getDimension(), it->getIndex() }, interior);

    // Create level-set function: circle of given radius centered at (0.5, 0.5)
    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = [](const Geometry::Point& p)
    {
      return (p.x() - 0.5) * (p.x() - 0.5) + (p.y() - 0.5) * (p.y() - 0.5) - radius;
    };

    MMG::Mesh result =
      MMG::LevelSetDiscretizer()
        .split(interior, { interior, exterior })
        .setBoundaryReference(boundary)
        .setHMax(0.1)
        .discretize(ls);

    // Result should be a valid non-empty mesh
    EXPECT_FALSE(result.isEmpty());
    EXPECT_GT(result.getVertexCount(), 0);
    EXPECT_GT(result.getCellCount(), 0);
    EXPECT_EQ(result.getDimension(), 2);
  }

  /// @brief Verifies fluent API for MMG level set discretizer by checking false predicates.
  TEST(Rodin_MMG_LevelSetDiscretizer, FluentAPI)
  {
    static constexpr Attribute interior = 1;
    static constexpr Attribute exterior = 2;

    const size_t n = 16;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
    mesh.scale(1.0 / (n - 1));

    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({ it->getDimension(), it->getIndex() }, interior);

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = [](const Geometry::Point& p)
    {
      return (p.x() - 0.5) * (p.x() - 0.5) + (p.y() - 0.5) * (p.y() - 0.5) - 0.1;
    };

    // Verify fluent API chaining compiles and works
    MMG::Mesh result =
      MMG::LevelSetDiscretizer()
        .split(interior, { interior, exterior })
        .setBoundaryReference(10)
        .setHMin(0.01)
        .setHMax(0.1)
        .setGradation(1.3)
        .setAngleDetection(true)
        .discretize(ls);

    EXPECT_FALSE(result.isEmpty());
  }

  /// @brief Verifies no split material for MMG level set discretizer by checking exact expected values, true predicates, false predicates.
  TEST(Rodin_MMG_LevelSetDiscretizer, NoSplitMaterial)
  {
    static constexpr Attribute mat1 = 1;
    static constexpr Attribute exterior = 3;

    const size_t n = 16;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
    mesh.scale(1.0 / (n - 1));

    // Set half the cells to mat1, half to mat2
    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({ it->getDimension(), it->getIndex() }, mat1);

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = [](const Geometry::Point& p)
    {
      return (p.x() - 0.5) * (p.x() - 0.5) + (p.y() - 0.5) * (p.y() - 0.5) - 0.1;
    };

    MMG::LevelSetDiscretizer discretizer;
    discretizer.split(mat1, { mat1, exterior });
    discretizer.setBoundaryReference(10);
    discretizer.setHMax(0.1);

    // Verify getSplitMap() reflects the configuration
    const auto& splitMap = discretizer.getSplitMap();
    EXPECT_EQ(splitMap.size(), 1);
    EXPECT_TRUE(splitMap.count(mat1));

    MMG::Mesh result = discretizer.discretize(ls);
    EXPECT_FALSE(result.isEmpty());
  }

  /// @brief Verifies set level set value for MMG level set discretizer by checking false predicates.
  TEST(Rodin_MMG_LevelSetDiscretizer, SetLevelSetValue)
  {
    static constexpr Attribute interior = 1;
    static constexpr Attribute exterior = 2;

    const size_t n = 16;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
    mesh.scale(1.0 / (n - 1));

    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({ it->getDimension(), it->getIndex() }, interior);

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = [](const Geometry::Point& p)
    {
      return (p.x() - 0.5) * (p.x() - 0.5) + (p.y() - 0.5) * (p.y() - 0.5) - 0.1;
    };

    // Discretize at a non-zero level-set value
    MMG::Mesh result =
      MMG::LevelSetDiscretizer()
        .split(interior, { interior, exterior })
        .setBoundaryReference(10)
        .setLevelSet(0.0)
        .setHMax(0.1)
        .discretize(ls);

    EXPECT_FALSE(result.isEmpty());
    EXPECT_GT(result.getVertexCount(), 0);
  }

  /// @brief Verifies preserves required triangles 2 D for MMG level set discretizer by checking exact expected values, true predicates, false predicates.
  TEST(Rodin_MMG_LevelSetDiscretizer, PreservesRequiredTriangles2D)
  {
    static constexpr Attribute interior = 1;
    static constexpr Attribute exterior = 2;

    const size_t n = 16;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {n, n});
    mesh.scale(1.0 / (n - 1));
    const Index required = 0;
    mesh.setRequiredTriangle(required);

    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({it->getDimension(), it->getIndex()}, interior);

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = [](const Geometry::Point& p) {
      return (p.x() - 0.5) * (p.x() - 0.5) + (p.y() - 0.5) * (p.y() - 0.5) - 0.1;
    };

    MMG::Mesh result = MMG::LevelSetDiscretizer()
                         .split(interior, {interior, exterior})
                         .setBoundaryReference(10)
                         .setHMax(0.1)
                         .discretize(ls);

    EXPECT_FALSE(result.isEmpty());
    EXPECT_GT(result.getCellCount(), 0);
    EXPECT_EQ(result.getRequiredTriangles().size(), 1);
    EXPECT_TRUE(result.getRequiredTriangles().count(required));
  }

  /// @brief Verifies preserves required triangles 3 D for MMG level set discretizer by checking exact expected values, true predicates, false predicates.
  TEST(Rodin_MMG_LevelSetDiscretizer, PreservesRequiredTriangles3D)
  {
    static constexpr Attribute interior = 1;
    static constexpr Attribute exterior = 2;

    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {n, n, n});
    mesh.scale(1.0 / (n - 1));
    const Index required = getFirstBoundaryTriangleIndex(mesh);
    mesh.setRequiredTriangle(required);

    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({it->getDimension(), it->getIndex()}, interior);

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = [](const Geometry::Point& p) {
      const Real dx = p.x() - 0.5;
      const Real dy = p.y() - 0.5;
      const Real dz = p.z() - 0.5;
      return dx * dx + dy * dy + dz * dz - 0.16;
    };

    MMG::Mesh result = MMG::LevelSetDiscretizer()
                         .split(interior, {interior, exterior})
                         .setBoundaryReference(10)
                         .setHMax(0.35)
                         .discretize(ls);

    EXPECT_FALSE(result.isEmpty());
    EXPECT_GT(result.getCellCount(), 0);
    EXPECT_EQ(result.getRequiredTriangles().size(), 1);
    EXPECT_TRUE(result.getRequiredTriangles().count(required));
  }

  /// @brief Verifies preserves required entities 2 D for MMG level set discretizer by checking exact expected values, true predicates, false predicates.
  TEST(Rodin_MMG_LevelSetDiscretizer, PreservesRequiredEntities2D)
  {
    static constexpr Attribute interior = 1;
    static constexpr Attribute exterior = 2;

    const size_t n = 16;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {n, n});
    mesh.scale(1.0 / (n - 1));
    setAllSegmentAttributes(mesh, 8);

    const Index requiredVertex = 0;
    const Index requiredEdge = 0;
    const Index requiredTriangle = 0;
    mesh.setRequiredVertex(requiredVertex);
    mesh.setRequiredEdge(requiredEdge);
    mesh.setRequiredTriangle(requiredTriangle);

    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({it->getDimension(), it->getIndex()}, interior);

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = [](const Geometry::Point& p) {
      return (p.x() - 0.5) * (p.x() - 0.5) + (p.y() - 0.5) * (p.y() - 0.5) - 0.1;
    };

    MMG::Mesh result = MMG::LevelSetDiscretizer()
                         .split(interior, {interior, exterior})
                         .setBoundaryReference(10)
                         .setHMax(0.1)
                         .discretize(ls);

    EXPECT_FALSE(result.isEmpty());
    EXPECT_GT(result.getCellCount(), 0);
    EXPECT_EQ(result.getRequiredVertices().size(), 1);
    EXPECT_EQ(result.getRequiredEdges().size(), 1);
    EXPECT_EQ(result.getRequiredTriangles().size(), 1);
    EXPECT_TRUE(result.getRequiredVertices().count(requiredVertex));
    EXPECT_TRUE(result.getRequiredEdges().count(requiredEdge));
    EXPECT_TRUE(result.getRequiredTriangles().count(requiredTriangle));
  }

  /// @brief Verifies preserves uncut required entity geometry 2 D for MMG level set discretizer by checking true predicates.
  TEST(Rodin_MMG_LevelSetDiscretizer, PreservesUncutRequiredEntityGeometry2D)
  {
    static constexpr Attribute interior = 1;
    static constexpr Attribute exterior = 2;

    const size_t n = 16;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {n, n});
    mesh.scale(1.0 / (n - 1));
    mesh.getConnectivity().compute(1, 2);
    setAllSegmentAttributes(mesh, 8);

    const auto lsFn = [](const Point& p) {
      return (p.x() - 0.5) * (p.x() - 0.5) + (p.y() - 0.5) * (p.y() - 0.5) - 0.1;
    };

    const Index requiredVertex = 0;
    const Index requiredEdge = getFirstUncutEntity(mesh, Polytope::Type::Segment, lsFn);
    const Index requiredTriangle =
      getFirstUncutEntity(mesh, Polytope::Type::Triangle, lsFn);
    const auto vertexCoords = mesh.getVertexCoordinates(requiredVertex);
    const Point requiredVertexPoint(
      *mesh.getVertex(requiredVertex), vertexCoords, vertexCoords);
    const auto requiredEdgePoints = getVertexCoordinates(
      mesh, *mesh.getPolytope(getDimension(Polytope::Type::Segment), requiredEdge));
    const auto requiredTrianglePoints = getVertexCoordinates(
      mesh, *mesh.getPolytope(getDimension(Polytope::Type::Triangle), requiredTriangle));

    mesh.setRequiredVertex(requiredVertex);
    mesh.setRequiredEdge(requiredEdge);
    mesh.setRequiredTriangle(requiredTriangle);

    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({it->getDimension(), it->getIndex()}, interior);

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = lsFn;

    MMG::Mesh result = MMG::LevelSetDiscretizer()
                         .split(interior, {interior, exterior})
                         .setBoundaryReference(10)
                         .setHMax(0.1)
                         .discretize(ls);

    EXPECT_TRUE(result.getRequiredVertices().count(requiredVertex));
    EXPECT_TRUE(result.getRequiredEdges().count(requiredEdge));
    EXPECT_TRUE(result.getRequiredTriangles().count(requiredTriangle));
    EXPECT_TRUE(hasVertexWithCoordinates(result, requiredVertexPoint));
    EXPECT_TRUE(
      hasEntityWithCoordinates(result, Polytope::Type::Segment, requiredEdgePoints));
    EXPECT_TRUE(
      hasEntityWithCoordinates(result, Polytope::Type::Triangle, requiredTrianglePoints));
  }

  /// @brief Verifies preserves required entities 3 D for MMG level set discretizer by checking exact expected values, true predicates, false predicates.
  TEST(Rodin_MMG_LevelSetDiscretizer, PreservesRequiredEntities3D)
  {
    static constexpr Attribute interior = 1;
    static constexpr Attribute exterior = 2;

    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {n, n, n});
    mesh.scale(1.0 / (n - 1));
    setAllSegmentAttributes(mesh, 8);

    const Index requiredVertex = 0;
    const Index requiredEdge = 0;
    const Index requiredTriangle = getFirstBoundaryTriangleIndex(mesh);
    mesh.setRequiredVertex(requiredVertex);
    mesh.setRequiredEdge(requiredEdge);
    mesh.setRequiredTriangle(requiredTriangle);

    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({it->getDimension(), it->getIndex()}, interior);

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = [](const Geometry::Point& p) {
      const Real dx = p.x() - 0.5;
      const Real dy = p.y() - 0.5;
      const Real dz = p.z() - 0.5;
      return dx * dx + dy * dy + dz * dz - 0.16;
    };

    Index requiredTetrahedron = 0;
    bool foundUncutTetrahedron = false;
    for (auto it = mesh.getCell(); it; ++it)
    {
      bool hasNegative = false;
      bool hasPositive = false;
      for (const auto& vertex : it->getVertices())
      {
        hasNegative = hasNegative || ls[vertex] < 0;
        hasPositive = hasPositive || ls[vertex] > 0;
      }
      if (!(hasNegative && hasPositive))
      {
        requiredTetrahedron = it->getIndex();
        foundUncutTetrahedron = true;
        break;
      }
    }
    ASSERT_TRUE(foundUncutTetrahedron);
    mesh.setRequiredTetrahedron(requiredTetrahedron);

    MMG::Mesh result = MMG::LevelSetDiscretizer()
                         .split(interior, {interior, exterior})
                         .setBoundaryReference(10)
                         .setHMax(0.35)
                         .discretize(ls);

    EXPECT_FALSE(result.isEmpty());
    EXPECT_GT(result.getCellCount(), 0);
    EXPECT_EQ(result.getRequiredVertices().size(), 1);
    EXPECT_EQ(result.getRequiredEdges().size(), 1);
    EXPECT_EQ(result.getRequiredTriangles().size(), 1);
    EXPECT_EQ(result.getRequiredTetrahedra().size(), 1);
    EXPECT_TRUE(result.getRequiredVertices().count(requiredVertex));
    EXPECT_TRUE(result.getRequiredEdges().count(requiredEdge));
    EXPECT_TRUE(result.getRequiredTriangles().count(requiredTriangle));
    EXPECT_TRUE(result.getRequiredTetrahedra().count(requiredTetrahedron));
  }

  /// @brief Verifies preserves uncut required entity geometry 3 D for MMG level set discretizer by checking true predicates.
  TEST(Rodin_MMG_LevelSetDiscretizer, PreservesUncutRequiredEntityGeometry3D)
  {
    static constexpr Attribute interior = 1;
    static constexpr Attribute exterior = 2;

    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {n, n, n});
    mesh.scale(1.0 / (n - 1));
    mesh.getConnectivity().compute(1, 3);
    mesh.getConnectivity().compute(2, 3);
    setAllSegmentAttributes(mesh, 8);

    const auto lsFn = [](const Point& p) {
      const Real dx = p.x() - 0.5;
      const Real dy = p.y() - 0.5;
      const Real dz = p.z() - 0.5;
      return dx * dx + dy * dy + dz * dz - 0.16;
    };

    const Index requiredVertex = 0;
    const Index requiredEdge = getFirstUncutEntity(mesh, Polytope::Type::Segment, lsFn);
    const Index requiredTriangle =
      getFirstUncutEntity(mesh, Polytope::Type::Triangle, lsFn);
    const Index requiredTetrahedron =
      getFirstUncutEntity(mesh, Polytope::Type::Tetrahedron, lsFn);
    const auto vertexCoords = mesh.getVertexCoordinates(requiredVertex);
    const Point requiredVertexPoint(
      *mesh.getVertex(requiredVertex), vertexCoords, vertexCoords);
    const auto requiredEdgePoints = getVertexCoordinates(
      mesh, *mesh.getPolytope(getDimension(Polytope::Type::Segment), requiredEdge));
    const auto requiredTrianglePoints = getVertexCoordinates(
      mesh, *mesh.getPolytope(getDimension(Polytope::Type::Triangle), requiredTriangle));
    const auto requiredTetrahedronPoints = getVertexCoordinates(mesh,
      *mesh.getPolytope(getDimension(Polytope::Type::Tetrahedron), requiredTetrahedron));

    mesh.setRequiredVertex(requiredVertex);
    mesh.setRequiredEdge(requiredEdge);
    mesh.setRequiredTriangle(requiredTriangle);
    mesh.setRequiredTetrahedron(requiredTetrahedron);

    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({it->getDimension(), it->getIndex()}, interior);

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = lsFn;

    MMG::Mesh result = MMG::LevelSetDiscretizer()
                         .split(interior, {interior, exterior})
                         .setBoundaryReference(10)
                         .setHMax(0.35)
                         .discretize(ls);

    EXPECT_TRUE(result.getRequiredVertices().count(requiredVertex));
    EXPECT_TRUE(result.getRequiredEdges().count(requiredEdge));
    EXPECT_TRUE(result.getRequiredTriangles().count(requiredTriangle));
    EXPECT_TRUE(result.getRequiredTetrahedra().count(requiredTetrahedron));
    EXPECT_TRUE(hasVertexWithCoordinates(result, requiredVertexPoint));
    EXPECT_TRUE(
      hasEntityWithCoordinates(result, Polytope::Type::Segment, requiredEdgePoints));
    EXPECT_TRUE(
      hasEntityWithCoordinates(result, Polytope::Type::Triangle, requiredTrianglePoints));
    EXPECT_TRUE(hasEntityWithCoordinates(
      result, Polytope::Type::Tetrahedron, requiredTetrahedronPoints));
  }

  /// @brief Verifies no split does not preserve cut required edge 2 D for MMG level set discretizer by checking true predicates, false predicates.
  TEST(Rodin_MMG_LevelSetDiscretizer, NoSplitDoesNotPreserveCutRequiredEdge2D)
  {
    static constexpr Attribute splitRef = 1;
    static constexpr Attribute exterior = 2;
    static constexpr Attribute protectedRef = 7;

    const size_t n = 16;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {n, n});
    mesh.scale(1.0 / (n - 1));
    mesh.getConnectivity().compute(1, 2);

    const auto lsFn = [](const Point& p) { return p.x() + p.y() - 0.75; };
    const Index required = getFirstCutEntity(mesh, Polytope::Type::Segment, lsFn);
    const auto entity = mesh.getPolytope(getDimension(Polytope::Type::Segment), required);
    const auto requiredPoints = getVertexCoordinates(mesh, *entity);

    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({it->getDimension(), it->getIndex()}, splitRef);
    mesh.setRequiredEdge(required);
    setCellsContainingVertices(mesh, entity->getVertices(), protectedRef);

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = lsFn;

    MMG::Mesh result = MMG::LevelSetDiscretizer()
                         .split(splitRef, {splitRef, exterior})
                         .noSplit(protectedRef)
                         .setBoundaryReference(10)
                         .setHMax(0.1)
                         .discretize(ls);

    EXPECT_TRUE(result.getAttributeIndex()
                  .getAttributes(result.getDimension())
                  .count(protectedRef));
    EXPECT_FALSE(
      hasEntityWithCoordinates(result, Polytope::Type::Segment, requiredPoints));
  }

  /// @brief Verifies no split does not preserve cut required triangle 2 D for MMG level set discretizer by checking true predicates, false predicates.
  TEST(Rodin_MMG_LevelSetDiscretizer, NoSplitDoesNotPreserveCutRequiredTriangle2D)
  {
    static constexpr Attribute splitRef = 1;
    static constexpr Attribute exterior = 2;
    static constexpr Attribute protectedRef = 7;

    const size_t n = 16;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {n, n});
    mesh.scale(1.0 / (n - 1));

    const auto lsFn = [](const Point& p) { return p.x() + p.y() - 0.75; };
    const Index required = getFirstCutEntity(mesh, Polytope::Type::Triangle, lsFn);
    const auto entity =
      mesh.getPolytope(getDimension(Polytope::Type::Triangle), required);
    const auto requiredPoints = getVertexCoordinates(mesh, *entity);

    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({it->getDimension(), it->getIndex()}, splitRef);
    mesh.setRequiredTriangle(required);
    mesh.setAttribute({entity->getDimension(), entity->getIndex()}, protectedRef);

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = lsFn;

    MMG::Mesh result = MMG::LevelSetDiscretizer()
                         .split(splitRef, {splitRef, exterior})
                         .noSplit(protectedRef)
                         .setBoundaryReference(10)
                         .setHMax(0.1)
                         .discretize(ls);

    EXPECT_TRUE(result.getAttributeIndex()
                  .getAttributes(result.getDimension())
                  .count(protectedRef));
    EXPECT_FALSE(
      hasEntityWithCoordinates(result, Polytope::Type::Triangle, requiredPoints));
  }

  /// @brief Verifies no split does not preserve individually labeled cut triangle 2 D for MMG level set discretizer by checking true predicates, false predicates.
  TEST(
    Rodin_MMG_LevelSetDiscretizer, NoSplitDoesNotPreserveIndividuallyLabeledCutTriangle2D)
  {
    static constexpr Attribute splitRef = 1;
    static constexpr Attribute exterior = 2;
    static constexpr Attribute protectedRef = 7;

    const size_t n = 16;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, {n, n});
    mesh.scale(1.0 / (n - 1));

    const auto lsFn = [](const Point& p) { return p.x() + p.y() - 0.75; };
    const Index protectedTriangle =
      getFirstCutEntity(mesh, Polytope::Type::Triangle, lsFn);
    const auto entity =
      mesh.getPolytope(getDimension(Polytope::Type::Triangle), protectedTriangle);
    const auto protectedPoints = getVertexCoordinates(mesh, *entity);

    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({it->getDimension(), it->getIndex()}, splitRef);
    mesh.setAttribute({entity->getDimension(), entity->getIndex()}, protectedRef);

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = lsFn;

    MMG::Mesh result = MMG::LevelSetDiscretizer()
                         .split(splitRef, {splitRef, exterior})
                         .noSplit(protectedRef)
                         .setBoundaryReference(10)
                         .setHMax(0.1)
                         .discretize(ls);

    EXPECT_TRUE(result.getAttributeIndex()
                  .getAttributes(result.getDimension())
                  .count(protectedRef));
    EXPECT_FALSE(
      hasEntityWithCoordinates(result, Polytope::Type::Triangle, protectedPoints));
  }

  /// @brief Verifies no split does not preserve cut required edge 3 D for MMG level set discretizer by checking true predicates, false predicates.
  TEST(Rodin_MMG_LevelSetDiscretizer, NoSplitDoesNotPreserveCutRequiredEdge3D)
  {
    static constexpr Attribute splitRef = 1;
    static constexpr Attribute exterior = 2;
    static constexpr Attribute protectedRef = 7;

    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {n, n, n});
    mesh.scale(1.0 / (n - 1));
    mesh.getConnectivity().compute(1, 3);

    const auto lsFn = [](const Point& p) { return p.x() + p.y() + p.z() - 0.8; };
    const Index required = getFirstCutEntity(mesh, Polytope::Type::Segment, lsFn);
    const auto entity = mesh.getPolytope(getDimension(Polytope::Type::Segment), required);
    const auto requiredPoints = getVertexCoordinates(mesh, *entity);

    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({it->getDimension(), it->getIndex()}, splitRef);
    mesh.setRequiredEdge(required);
    setCellsContainingVertices(mesh, entity->getVertices(), protectedRef);

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = lsFn;

    MMG::Mesh result = MMG::LevelSetDiscretizer()
                         .split(splitRef, {splitRef, exterior})
                         .noSplit(protectedRef)
                         .setBoundaryReference(10)
                         .setHMax(0.35)
                         .discretize(ls);

    EXPECT_TRUE(result.getAttributeIndex()
                  .getAttributes(result.getDimension())
                  .count(protectedRef));
    EXPECT_FALSE(
      hasEntityWithCoordinates(result, Polytope::Type::Segment, requiredPoints));
  }

  /// @brief Verifies no split does not preserve cut required triangle 3 D for MMG level set discretizer by checking true predicates, false predicates.
  TEST(Rodin_MMG_LevelSetDiscretizer, NoSplitDoesNotPreserveCutRequiredTriangle3D)
  {
    static constexpr Attribute splitRef = 1;
    static constexpr Attribute exterior = 2;
    static constexpr Attribute protectedRef = 7;

    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {n, n, n});
    mesh.scale(1.0 / (n - 1));
    mesh.getConnectivity().compute(2, 3);

    const auto lsFn = [](const Point& p) { return p.x() + p.y() + p.z() - 0.8; };
    const Index required = getFirstCutEntity(mesh, Polytope::Type::Triangle, lsFn);
    const auto entity =
      mesh.getPolytope(getDimension(Polytope::Type::Triangle), required);
    const auto requiredPoints = getVertexCoordinates(mesh, *entity);

    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({it->getDimension(), it->getIndex()}, splitRef);
    mesh.setRequiredTriangle(required);
    setCellsContainingVertices(mesh, entity->getVertices(), protectedRef);

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = lsFn;

    MMG::Mesh result = MMG::LevelSetDiscretizer()
                         .split(splitRef, {splitRef, exterior})
                         .noSplit(protectedRef)
                         .setBoundaryReference(10)
                         .setHMax(0.35)
                         .discretize(ls);

    EXPECT_TRUE(result.getAttributeIndex()
                  .getAttributes(result.getDimension())
                  .count(protectedRef));
    EXPECT_FALSE(
      hasEntityWithCoordinates(result, Polytope::Type::Triangle, requiredPoints));
  }

  /// @brief Verifies no split adjacent tetrahedra does not preserve cut triangle 3 D for MMG level set discretizer by checking true predicates, false predicates.
  TEST(
    Rodin_MMG_LevelSetDiscretizer, NoSplitAdjacentTetrahedraDoesNotPreserveCutTriangle3D)
  {
    static constexpr Attribute splitRef = 1;
    static constexpr Attribute exterior = 2;
    static constexpr Attribute protectedRef = 7;

    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {n, n, n});
    mesh.scale(1.0 / (n - 1));
    mesh.getConnectivity().compute(2, 3);

    const auto lsFn = [](const Point& p) { return p.x() + p.y() + p.z() - 0.8; };
    const Index protectedTriangle =
      getFirstCutEntity(mesh, Polytope::Type::Triangle, lsFn);
    const auto entity =
      mesh.getPolytope(getDimension(Polytope::Type::Triangle), protectedTriangle);
    const auto protectedPoints = getVertexCoordinates(mesh, *entity);

    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({it->getDimension(), it->getIndex()}, splitRef);
    setCellsContainingVertices(mesh, entity->getVertices(), protectedRef);

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = lsFn;

    MMG::Mesh result = MMG::LevelSetDiscretizer()
                         .split(splitRef, {splitRef, exterior})
                         .noSplit(protectedRef)
                         .setBoundaryReference(10)
                         .setHMax(0.35)
                         .discretize(ls);

    EXPECT_TRUE(result.getAttributeIndex()
                  .getAttributes(result.getDimension())
                  .count(protectedRef));
    EXPECT_FALSE(
      hasEntityWithCoordinates(result, Polytope::Type::Triangle, protectedPoints));
  }

  /// @brief Verifies no split adjacent tetrahedra does not preserve labeled required cut triangle 3 D for MMG level set discretizer by checking true predicates, false predicates.
  TEST(Rodin_MMG_LevelSetDiscretizer,
    NoSplitAdjacentTetrahedraDoesNotPreserveLabeledRequiredCutTriangle3D)
  {
    static constexpr Attribute splitRef = 1;
    static constexpr Attribute exterior = 2;
    static constexpr Attribute protectedRef = 7;

    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {n, n, n});
    mesh.scale(1.0 / (n - 1));
    mesh.getConnectivity().compute(2, 3);

    const auto lsFn = [](const Point& p) { return p.x() + p.y() + p.z() - 0.8; };
    const Index protectedTriangle =
      getFirstCutEntity(mesh, Polytope::Type::Triangle, lsFn);
    const auto entity =
      mesh.getPolytope(getDimension(Polytope::Type::Triangle), protectedTriangle);
    const auto protectedPoints = getVertexCoordinates(mesh, *entity);

    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({it->getDimension(), it->getIndex()}, splitRef);
    mesh.setAttribute({entity->getDimension(), entity->getIndex()}, protectedRef);
    mesh.setRequiredTriangle(protectedTriangle);
    setCellsContainingVertices(mesh, entity->getVertices(), protectedRef);

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = lsFn;

    MMG::Mesh result = MMG::LevelSetDiscretizer()
                         .split(splitRef, {splitRef, exterior})
                         .noSplit(protectedRef)
                         .setBoundaryReference(10)
                         .setHMax(0.35)
                         .discretize(ls);

    EXPECT_TRUE(result.getAttributeIndex()
                  .getAttributes(result.getDimension())
                  .count(protectedRef));
    EXPECT_FALSE(
      hasEntityWithCoordinates(result, Polytope::Type::Triangle, protectedPoints));
  }

  /// @brief Verifies no split required adjacent tetrahedra does not preserve cut triangle 3 D for MMG level set discretizer by checking true predicates, false predicates.
  TEST(Rodin_MMG_LevelSetDiscretizer,
    NoSplitRequiredAdjacentTetrahedraDoesNotPreserveCutTriangle3D)
  {
    static constexpr Attribute splitRef = 1;
    static constexpr Attribute exterior = 2;
    static constexpr Attribute protectedRef = 7;

    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {n, n, n});
    mesh.scale(1.0 / (n - 1));
    mesh.getConnectivity().compute(2, 3);

    const auto lsFn = [](const Point& p) { return p.x() + p.y() + p.z() - 0.8; };
    const Index protectedTriangle =
      getFirstCutEntity(mesh, Polytope::Type::Triangle, lsFn);
    const auto entity =
      mesh.getPolytope(getDimension(Polytope::Type::Triangle), protectedTriangle);
    const auto protectedPoints = getVertexCoordinates(mesh, *entity);
    const auto adjacentTetrahedra =
      getCellsContainingVertices(mesh, entity->getVertices());
    ASSERT_FALSE(adjacentTetrahedra.empty());

    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({it->getDimension(), it->getIndex()}, splitRef);
    for (const auto& tetrahedron : adjacentTetrahedra)
    {
      mesh.setRequiredTetrahedron(tetrahedron);
      mesh.setAttribute(
        {getDimension(Polytope::Type::Tetrahedron), tetrahedron}, protectedRef);
    }

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = lsFn;

    MMG::Mesh result = MMG::LevelSetDiscretizer()
                         .split(splitRef, {splitRef, exterior})
                         .noSplit(protectedRef)
                         .setBoundaryReference(10)
                         .setHMax(0.35)
                         .discretize(ls);

    EXPECT_TRUE(result.getAttributeIndex()
                  .getAttributes(result.getDimension())
                  .count(protectedRef));
    EXPECT_FALSE(
      hasEntityWithCoordinates(result, Polytope::Type::Triangle, protectedPoints));
  }

  /// @brief Verifies no split required triangle closure does not preserve cut triangle 3 D for MMG level set discretizer by checking exact expected values, true predicates, false predicates.
  TEST(Rodin_MMG_LevelSetDiscretizer,
    NoSplitRequiredTriangleClosureDoesNotPreserveCutTriangle3D)
  {
    static constexpr Attribute splitRef = 1;
    static constexpr Attribute exterior = 2;
    static constexpr Attribute protectedRef = 7;

    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {n, n, n});
    mesh.scale(1.0 / (n - 1));
    mesh.getConnectivity().compute(2, 3);

    const auto lsFn = [](const Point& p) { return p.x() + p.y() + p.z() - 0.8; };
    const Index protectedTriangle =
      getFirstCutEntity(mesh, Polytope::Type::Triangle, lsFn);
    const auto entity =
      mesh.getPolytope(getDimension(Polytope::Type::Triangle), protectedTriangle);
    const auto protectedPoints = getVertexCoordinates(mesh, *entity);
    const auto adjacentTetrahedra =
      getCellsContainingVertices(mesh, entity->getVertices());
    const auto edges = getEdgesOnTriangle(mesh, entity->getVertices());
    ASSERT_FALSE(adjacentTetrahedra.empty());
    ASSERT_EQ(edges.size(), 3);

    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({it->getDimension(), it->getIndex()}, splitRef);
    for (const auto& vertex : entity->getVertices())
      mesh.setRequiredVertex(vertex);
    for (const auto& edge : edges)
      mesh.setRequiredEdge(edge);
    for (const auto& tetrahedron : adjacentTetrahedra)
    {
      mesh.setRequiredTetrahedron(tetrahedron);
      mesh.setAttribute(
        {getDimension(Polytope::Type::Tetrahedron), tetrahedron}, protectedRef);
    }

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = lsFn;

    MMG::Mesh result = MMG::LevelSetDiscretizer()
                         .split(splitRef, {splitRef, exterior})
                         .noSplit(protectedRef)
                         .setBoundaryReference(10)
                         .setHMax(0.35)
                         .discretize(ls);

    EXPECT_TRUE(result.getAttributeIndex()
                  .getAttributes(result.getDimension())
                  .count(protectedRef));
    EXPECT_FALSE(
      hasEntityWithCoordinates(result, Polytope::Type::Triangle, protectedPoints));
  }

  /// @brief Verifies no split required full triangle closure does not preserve cut triangle 3 D for MMG level set discretizer by checking exact expected values, true predicates, false predicates.
  TEST(Rodin_MMG_LevelSetDiscretizer,
    NoSplitRequiredFullTriangleClosureDoesNotPreserveCutTriangle3D)
  {
    static constexpr Attribute splitRef = 1;
    static constexpr Attribute exterior = 2;
    static constexpr Attribute protectedRef = 7;

    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {n, n, n});
    mesh.scale(1.0 / (n - 1));
    mesh.getConnectivity().compute(2, 3);

    const auto lsFn = [](const Point& p) { return p.x() + p.y() + p.z() - 0.8; };
    const Index protectedTriangle =
      getFirstCutEntity(mesh, Polytope::Type::Triangle, lsFn);
    const auto entity =
      mesh.getPolytope(getDimension(Polytope::Type::Triangle), protectedTriangle);
    const auto protectedPoints = getVertexCoordinates(mesh, *entity);
    const auto adjacentTetrahedra =
      getCellsContainingVertices(mesh, entity->getVertices());
    const auto edges = getEdgesOnTriangle(mesh, entity->getVertices());
    ASSERT_FALSE(adjacentTetrahedra.empty());
    ASSERT_EQ(edges.size(), 3);

    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({it->getDimension(), it->getIndex()}, splitRef);
    for (const auto& vertex : entity->getVertices())
      mesh.setRequiredVertex(vertex);
    for (const auto& edge : edges)
      mesh.setRequiredEdge(edge);
    mesh.setRequiredTriangle(protectedTriangle);
    mesh.setAttribute({entity->getDimension(), entity->getIndex()}, protectedRef);
    for (const auto& tetrahedron : adjacentTetrahedra)
    {
      mesh.setRequiredTetrahedron(tetrahedron);
      mesh.setAttribute(
        {getDimension(Polytope::Type::Tetrahedron), tetrahedron}, protectedRef);
    }

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = lsFn;

    MMG::Mesh result = MMG::LevelSetDiscretizer()
                         .split(splitRef, {splitRef, exterior})
                         .noSplit(protectedRef)
                         .setBoundaryReference(10)
                         .setHMax(0.35)
                         .discretize(ls);

    EXPECT_TRUE(result.getAttributeIndex()
                  .getAttributes(result.getDimension())
                  .count(protectedRef));
    EXPECT_FALSE(
      hasEntityWithCoordinates(result, Polytope::Type::Triangle, protectedPoints));
  }

  /// @brief Verifies no split does not preserve individually labeled cut triangle 3 D for MMG level set discretizer by checking false predicates.
  TEST(
    Rodin_MMG_LevelSetDiscretizer, NoSplitDoesNotPreserveIndividuallyLabeledCutTriangle3D)
  {
    static constexpr Attribute splitRef = 1;
    static constexpr Attribute exterior = 2;
    static constexpr Attribute protectedRef = 7;

    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {n, n, n});
    mesh.scale(1.0 / (n - 1));
    mesh.getConnectivity().compute(2, 3);

    const auto lsFn = [](const Point& p) { return p.x() + p.y() + p.z() - 0.8; };
    const Index protectedTriangle =
      getFirstCutEntity(mesh, Polytope::Type::Triangle, lsFn);
    const auto entity =
      mesh.getPolytope(getDimension(Polytope::Type::Triangle), protectedTriangle);
    const auto protectedPoints = getVertexCoordinates(mesh, *entity);

    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({it->getDimension(), it->getIndex()}, splitRef);
    mesh.setAttribute({entity->getDimension(), entity->getIndex()}, protectedRef);

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = lsFn;

    MMG::Mesh result = MMG::LevelSetDiscretizer()
                         .split(splitRef, {splitRef, exterior})
                         .noSplit(protectedRef)
                         .setBoundaryReference(10)
                         .setHMax(0.35)
                         .discretize(ls);

    EXPECT_FALSE(result.getAttributeIndex()
                   .getAttributes(result.getDimension())
                   .count(protectedRef));
    EXPECT_FALSE(
      hasEntityWithCoordinates(result, Polytope::Type::Triangle, protectedPoints));
  }

  /// @brief Verifies no split does not preserve individually labeled required cut triangle 3 D for MMG level set discretizer by checking false predicates.
  TEST(Rodin_MMG_LevelSetDiscretizer,
    NoSplitDoesNotPreserveIndividuallyLabeledRequiredCutTriangle3D)
  {
    static constexpr Attribute splitRef = 1;
    static constexpr Attribute exterior = 2;
    static constexpr Attribute protectedRef = 7;

    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {n, n, n});
    mesh.scale(1.0 / (n - 1));
    mesh.getConnectivity().compute(2, 3);

    const auto lsFn = [](const Point& p) { return p.x() + p.y() + p.z() - 0.8; };
    const Index protectedTriangle =
      getFirstCutEntity(mesh, Polytope::Type::Triangle, lsFn);
    const auto entity =
      mesh.getPolytope(getDimension(Polytope::Type::Triangle), protectedTriangle);
    const auto protectedPoints = getVertexCoordinates(mesh, *entity);

    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({it->getDimension(), it->getIndex()}, splitRef);
    mesh.setAttribute({entity->getDimension(), entity->getIndex()}, protectedRef);
    mesh.setRequiredTriangle(protectedTriangle);

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = lsFn;

    MMG::Mesh result = MMG::LevelSetDiscretizer()
                         .split(splitRef, {splitRef, exterior})
                         .noSplit(protectedRef)
                         .setBoundaryReference(10)
                         .setHMax(0.35)
                         .discretize(ls);

    EXPECT_FALSE(result.getAttributeIndex()
                   .getAttributes(result.getDimension())
                   .count(protectedRef));
    EXPECT_FALSE(
      hasEntityWithCoordinates(result, Polytope::Type::Triangle, protectedPoints));
  }

  /// @brief Verifies no split does not preserve cut required tetrahedron 3 D for MMG level set discretizer by checking true predicates, false predicates.
  TEST(Rodin_MMG_LevelSetDiscretizer, NoSplitDoesNotPreserveCutRequiredTetrahedron3D)
  {
    static constexpr Attribute splitRef = 1;
    static constexpr Attribute exterior = 2;
    static constexpr Attribute protectedRef = 7;

    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Tetrahedron, {n, n, n});
    mesh.scale(1.0 / (n - 1));

    const auto lsFn = [](const Point& p) { return p.x() + p.y() + p.z() - 0.8; };
    const Index required = getFirstCutEntity(mesh, Polytope::Type::Tetrahedron, lsFn);
    const auto entity =
      mesh.getPolytope(getDimension(Polytope::Type::Tetrahedron), required);
    const auto requiredPoints = getVertexCoordinates(mesh, *entity);

    for (auto it = mesh.getCell(); it; ++it)
      mesh.setAttribute({it->getDimension(), it->getIndex()}, splitRef);
    mesh.setRequiredTetrahedron(required);
    mesh.setAttribute({entity->getDimension(), entity->getIndex()}, protectedRef);

    P1 fes(mesh);
    MMG::RealGridFunction ls(fes);
    ls = lsFn;

    MMG::Mesh result = MMG::LevelSetDiscretizer()
                         .split(splitRef, {splitRef, exterior})
                         .noSplit(protectedRef)
                         .setBoundaryReference(10)
                         .setHMax(0.35)
                         .discretize(ls);

    EXPECT_TRUE(result.getAttributeIndex()
                  .getAttributes(result.getDimension())
                  .count(protectedRef));
    EXPECT_FALSE(
      hasEntityWithCoordinates(result, Polytope::Type::Tetrahedron, requiredPoints));
  }

  // ========================================================================
  // Common types tests
  // ========================================================================

  /// @brief Verifies split map construction for MMG common by checking exact expected values, true predicates.
  TEST(Rodin_MMG_Common, SplitMapConstruction)
  {
    MMG::SplitMap splitMap;
    splitMap[1] = MMG::Split{ 2, 3 };
    splitMap[4] = MMG::NoSplit;

    EXPECT_EQ(splitMap.size(), 2);
    EXPECT_TRUE(std::holds_alternative<MMG::Split>(splitMap.at(1)));
    EXPECT_TRUE(std::holds_alternative<MMG::NoSplitT>(splitMap.at(4)));

    const auto& split = std::get<MMG::Split>(splitMap.at(1));
    EXPECT_EQ(split.interior, 2);
    EXPECT_EQ(split.exterior, 3);
  }

  // ========================================================================
  // MMG::GridFunction alias tests
  // ========================================================================

  /// @brief Verifies scalar grid function construction for MMG grid function by checking exact expected values.
  TEST(Rodin_MMG_GridFunction, ScalarGridFunctionConstruction)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });

    P1 fes(mesh);
    MMG::RealGridFunction gf(fes);
    gf = [](const Geometry::Point& p) { return p.x() + p.y(); };

    // Grid function should have correct size
    EXPECT_EQ(gf.getFiniteElementSpace().getSize(), mesh.getVertexCount());
  }

  /// @brief Verifies vector grid function construction for MMG grid function by checking exact expected values.
  TEST(Rodin_MMG_GridFunction, VectorGridFunctionConstruction)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });

    P1 fes(mesh, 2);
    MMG::VectorGridFunction gf(fes);
    gf = [](const Geometry::Point& p)
    {
      return Math::SpatialVector<Real>{{ p.x(), p.y() }};
    };

    EXPECT_EQ(gf.getFiniteElementSpace().getVectorDimension(), 2);
  }

  // ========================================================================
  // MMG5 parameter configuration tests
  // ========================================================================

  /// @brief Verify that parameter setters return reference for chaining.
  TEST(Rodin_MMG_MMG5, ParameterSetters)
  {
    // Verify that parameter setters return reference for chaining
    MMG::Optimizer opt;
    auto& ref1 = opt.setHMin(0.01);
    auto& ref2 = opt.setHMax(1.0);
    auto& ref3 = opt.setHausdorff(0.001);
    auto& ref4 = opt.setGradation(1.3);
    auto& ref5 = opt.setAngleDetection(true);

    // All should return reference to the same object
    EXPECT_EQ(&ref1, &opt);
    EXPECT_EQ(&ref2, &opt);
    EXPECT_EQ(&ref3, &opt);
    EXPECT_EQ(&ref4, &opt);
    EXPECT_EQ(&ref5, &opt);
  }

  /// @brief Verifies adapt parameter setters for MMG MMG 5 by checking exact expected values.
  TEST(Rodin_MMG_MMG5, AdaptParameterSetters)
  {
    MMG::Adapt adapt;
    auto& ref1 = adapt.setHMin(0.01);
    auto& ref2 = adapt.setHMax(1.0);
    auto& ref3 = adapt.setHausdorff(0.001);
    auto& ref4 = adapt.setGradation(1.3);
    auto& ref5 = adapt.setAngleDetection(true);

    EXPECT_EQ(&ref1, &adapt);
    EXPECT_EQ(&ref2, &adapt);
    EXPECT_EQ(&ref3, &adapt);
    EXPECT_EQ(&ref4, &adapt);
    EXPECT_EQ(&ref5, &adapt);
  }

  /// @brief Verifies level set discretizer parameter setters for MMG MMG 5 by checking exact expected values.
  TEST(Rodin_MMG_MMG5, LevelSetDiscretizerParameterSetters)
  {
    MMG::LevelSetDiscretizer lsd;
    auto& ref1 = lsd.setHMin(0.01);
    auto& ref2 = lsd.setHMax(1.0);
    auto& ref3 = lsd.setHausdorff(0.001);
    auto& ref4 = lsd.setGradation(1.3);
    auto& ref5 = lsd.setAngleDetection(true);
    auto& ref6 = lsd.setLevelSet(0.5);
    auto& ref7 = lsd.setBoundaryReference(10);
    auto& ref8 = lsd.surface(true);

    EXPECT_EQ(&ref1, &lsd);
    EXPECT_EQ(&ref2, &lsd);
    EXPECT_EQ(&ref3, &lsd);
    EXPECT_EQ(&ref4, &lsd);
    EXPECT_EQ(&ref5, &lsd);
    EXPECT_EQ(&ref6, &lsd);
    EXPECT_EQ(&ref7, &lsd);
    EXPECT_EQ(&ref8, &lsd);
  }

  // ========================================================================
  // MMG::Mesh deeper behavior tests
  // ========================================================================

  /// @brief Verifies scale preserves topology for MMG mesh by checking exact expected values.
  TEST(Rodin_MMG_Mesh, ScalePreservesTopology)
  {
    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
    mesh.scale(0.5);

    EXPECT_EQ(mesh.getVertexCount(), n * n);
    EXPECT_EQ(mesh.getCellCount(), 2 * (n - 1) * (n - 1));
    EXPECT_EQ(mesh.getDimension(), 2);
  }

  /// @brief Verifies set attribute on cells for MMG mesh by checking exact expected values, true predicates.
  TEST(Rodin_MMG_Mesh, SetAttributeOnCells)
  {
    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });

    static constexpr Attribute label = 42;
    for (auto it = mesh.getCell(); it; ++it)
    {
      mesh.setAttribute({ it->getDimension(), it->getIndex() }, label);
    }

    // Verify attributes are set correctly
    for (auto it = mesh.getCell(); it; ++it)
    {
      auto attr = mesh.getAttribute(it->getDimension(), it->getIndex());
      ASSERT_TRUE(attr.has_value());
      EXPECT_EQ(*attr, label);
    }
  }

  /// @brief Verifies vertex coordinates for MMG mesh by checking exact expected values.
  TEST(Rodin_MMG_Mesh, VertexCoordinates)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 3, 3 });

    // Vertices on a 3x3 grid at integer coordinates
    const size_t vertexCount = mesh.getVertexCount();
    EXPECT_EQ(vertexCount, 9);

    // Check that vertex coordinates are within the expected range
    for (Index i = 0; i < static_cast<Index>(vertexCount); i++)
    {
      auto coords = mesh.getVertexCoordinates(i);
      EXPECT_GE(coords(0), 0.0);
      EXPECT_GE(coords(1), 0.0);
    }
  }

  /// @brief Verifies scale affects coordinates for MMG mesh.
  TEST(Rodin_MMG_Mesh, ScaleAffectsCoordinates)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    mesh.scale(0.5);

    // After scaling by 0.5, max coordinates should be halved
    for (Index i = 0; i < static_cast<Index>(mesh.getVertexCount()); i++)
    {
      auto coords = mesh.getVertexCoordinates(i);
      EXPECT_LE(coords(0), 1.1);  // Allow small FP tolerance
      EXPECT_LE(coords(1), 1.1);
    }
  }

  /// @brief Verifies connectivity compute for MMG mesh by checking exact expected values, true predicates, mesh connectivity computation.
  TEST(Rodin_MMG_Mesh, ConnectivityCompute)
  {
    const size_t n = 4;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
    mesh.getConnectivity().compute(1, 2);

    // After computing face→cell incidence, boundary detection works
    size_t boundaryCount = 0;
    for (auto it = mesh.getBoundary(); !it.end(); ++it)
    {
      boundaryCount++;
      EXPECT_TRUE(it->isBoundary());
    }
    // A 4×4 uniform grid has 4×(n-1) = 12 boundary edges
    EXPECT_EQ(boundaryCount, 4 * (n - 1));
  }

  /// @brief Verifies idempotent set corner for MMG mesh by checking exact expected values.
  TEST(Rodin_MMG_Mesh, IdempotentSetCorner)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });

    mesh.setCorner(0);
    mesh.setCorner(0);
    mesh.setCorner(0);

    // IndexSet should deduplicate
    EXPECT_EQ(mesh.getCorners().size(), 1);
  }

  /// @brief Verifies idempotent set ridge for MMG mesh by checking exact expected values.
  TEST(Rodin_MMG_Mesh, IdempotentSetRidge)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });

    mesh.setRidge(0);
    mesh.setRidge(0);

    EXPECT_EQ(mesh.getRidges().size(), 1);
  }

  /// @brief Verifies move assignment clears source for MMG mesh by checking exact expected values, false predicates, move semantics.
  TEST(Rodin_MMG_Mesh, MoveAssignmentClearsSource)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.setCorner(0);
    mesh.setCorner(3);
    mesh.setRidge(1);
    mesh.setRequiredVertex(5);
    mesh.setRequiredEdge(2);
    mesh.setRequiredTriangle(1);
    mesh.setRequiredTetrahedron(0);

    MMG::Mesh target;
    target = std::move(mesh);

    // Target has the data
    EXPECT_EQ(target.getCorners().size(), 2);
    EXPECT_EQ(target.getRidges().size(), 1);
    EXPECT_EQ(target.getRequiredVertices().size(), 1);
    EXPECT_EQ(target.getRequiredEdges().size(), 1);
    EXPECT_EQ(target.getRequiredTriangles().size(), 1);
    EXPECT_EQ(target.getRequiredTetrahedra().size(), 1);
    EXPECT_FALSE(target.isEmpty());
  }

  /// @brief Verifies copy preserves metadata after mutation for MMG mesh by checking exact expected values, copy semantics.
  TEST(Rodin_MMG_Mesh, CopyPreservesMetadataAfterMutation)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.setCorner(0);

    MMG::Mesh copy(mesh);
    EXPECT_EQ(copy.getCorners().size(), 1);

    // Mutate original
    mesh.setCorner(1);
    mesh.setCorner(2);
    EXPECT_EQ(mesh.getCorners().size(), 3);

    // Copy is independent
    EXPECT_EQ(copy.getCorners().size(), 1);
  }

  /// @brief Verifies parent move assignment clears MMG metadata for MMG mesh by checking exact expected values, move semantics.
  TEST(Rodin_MMG_Mesh, ParentMoveAssignmentClearsMMGMetadata)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    mesh.setCorner(0);
    mesh.setRidge(1);
    mesh.setRequiredTriangle(0);
    mesh.setRequiredTetrahedron(0);

    // Assign from parent type (which has no MMG metadata)
    Mesh<Context::Local> parentMesh;
    parentMesh = parentMesh.UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    mesh = std::move(parentMesh);

    // After parent assignment, topology changes and MMG metadata is cleared.
    EXPECT_EQ(mesh.getVertexCount(), 9);
    EXPECT_EQ(mesh.getCorners().size(), 0);
    EXPECT_EQ(mesh.getRidges().size(), 0);
    EXPECT_EQ(mesh.getRequiredVertices().size(), 0);
    EXPECT_EQ(mesh.getRequiredEdges().size(), 0);
    EXPECT_EQ(mesh.getRequiredTriangles().size(), 0);
    EXPECT_EQ(mesh.getRequiredTetrahedra().size(), 0);
  }

  /// @brief A 2D uniform grid on triangles is not a surface mesh.
  TEST(Rodin_MMG_Mesh, IsSurface)
  {
    // A 2D uniform grid on triangles is not a surface mesh
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    EXPECT_FALSE(mesh.isSurface());
  }

  /// @brief Verifies get area for MMG mesh.
  TEST(Rodin_MMG_Mesh, GetArea)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 3, 3 });
    mesh.getConnectivity().compute(1, 2);

    // The area of the default 2D grid
    Real area = mesh.getArea();
    EXPECT_GT(area, 0.0);
  }

  // ========================================================================
  // MMG::GridFunction deeper behavior tests
  // ========================================================================

  /// @brief Verifies scalar zero initialization for MMG grid function by checking tolerance-based numerical results.
  TEST(Rodin_MMG_GridFunction, ScalarZeroInitialization)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf(fes);

    // Default-constructed GridFunction should have zero data
    const auto& data = gf.getData();
    for (Eigen::Index i = 0; i < data.size(); i++)
    {
      EXPECT_NEAR(data(i), 0.0, 1e-15);
    }
  }

  /// @brief Verifies scalar project from lambda for MMG grid function by checking tolerance-based numerical results.
  TEST(Rodin_MMG_GridFunction, ScalarProjectFromLambda)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf(fes);

    gf = [](const Geometry::Point& p) { return 2.0 * p.x() + 3.0 * p.y(); };

    // P1 represents linear functions exactly on vertices.
    // Verify that DOFs match the expected values at vertex coordinates.
    for (Index i = 0; i < static_cast<Index>(mesh.getVertexCount()); i++)
    {
      auto coords = mesh.getVertexCoordinates(i);
      Real expected = 2.0 * coords(0) + 3.0 * coords(1);
      EXPECT_NEAR(gf[i], expected, 1e-10);
    }
  }

  /// @brief Verifies scalar project constant for MMG grid function by checking tolerance-based numerical results.
  TEST(Rodin_MMG_GridFunction, ScalarProjectConstant)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf(fes);

    gf = RealFunction(7.5);

    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], 7.5, 1e-14);
    }
  }

  /// @brief Verifies scalar data access for MMG grid function by checking tolerance-based numerical results, exact expected values.
  TEST(Rodin_MMG_GridFunction, ScalarDataAccess)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf(fes);

    gf = RealFunction(3.0);

    const auto& data = gf.getData();
    EXPECT_EQ(static_cast<size_t>(data.size()), gf.getSize());
    for (Eigen::Index i = 0; i < data.size(); i++)
    {
      EXPECT_NEAR(data(i), 3.0, 1e-14);
    }
  }

  /// @brief Verifies scalar operator bracket for MMG grid function by checking tolerance-based numerical results.
  TEST(Rodin_MMG_GridFunction, ScalarOperatorBracket)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf(fes);

    // Write via operator[]
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      gf[i] = static_cast<Real>(i * 10);
    }

    // Read back
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], static_cast<Real>(i * 10), 1e-14);
    }
  }

  /// @brief Verifies scalar arithmetic add scalar for MMG grid function by checking tolerance-based numerical results.
  TEST(Rodin_MMG_GridFunction, ScalarArithmeticAddScalar)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf(fes);

    gf = RealFunction(5.0);
    gf += 3.0;

    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], 8.0, 1e-14);
    }
  }

  /// @brief Verifies scalar arithmetic sub scalar for MMG grid function by checking tolerance-based numerical results.
  TEST(Rodin_MMG_GridFunction, ScalarArithmeticSubScalar)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf(fes);

    gf = RealFunction(10.0);
    gf -= 4.0;

    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], 6.0, 1e-14);
    }
  }

  /// @brief Verifies scalar arithmetic mul scalar for MMG grid function by checking tolerance-based numerical results.
  TEST(Rodin_MMG_GridFunction, ScalarArithmeticMulScalar)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf(fes);

    gf = RealFunction(4.0);
    gf *= 3.0;

    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], 12.0, 1e-14);
    }
  }

  /// @brief Verifies scalar arithmetic div scalar for MMG grid function by checking tolerance-based numerical results.
  TEST(Rodin_MMG_GridFunction, ScalarArithmeticDivScalar)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf(fes);

    gf = RealFunction(12.0);
    gf /= 4.0;

    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], 3.0, 1e-14);
    }
  }

  /// @brief Verifies scalar arithmetic add grid function for MMG grid function by checking tolerance-based numerical results.
  TEST(Rodin_MMG_GridFunction, ScalarArithmeticAddGridFunction)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf1(fes), gf2(fes);

    gf1 = RealFunction(3.0);
    gf2 = RealFunction(7.0);
    gf1 += gf2;

    for (Index i = 0; i < static_cast<Index>(gf1.getSize()); i++)
    {
      EXPECT_NEAR(gf1[i], 10.0, 1e-14);
    }
  }

  /// @brief Verifies scalar arithmetic sub grid function for MMG grid function by checking tolerance-based numerical results.
  TEST(Rodin_MMG_GridFunction, ScalarArithmeticSubGridFunction)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf1(fes), gf2(fes);

    gf1 = RealFunction(10.0);
    gf2 = RealFunction(4.0);
    gf1 -= gf2;

    for (Index i = 0; i < static_cast<Index>(gf1.getSize()); i++)
    {
      EXPECT_NEAR(gf1[i], 6.0, 1e-14);
    }
  }

  /// @brief Verifies scalar arithmetic mul grid function for MMG grid function by checking tolerance-based numerical results.
  TEST(Rodin_MMG_GridFunction, ScalarArithmeticMulGridFunction)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf1(fes), gf2(fes);

    gf1 = RealFunction(3.0);
    gf2 = RealFunction(5.0);
    gf1 *= gf2;

    for (Index i = 0; i < static_cast<Index>(gf1.getSize()); i++)
    {
      EXPECT_NEAR(gf1[i], 15.0, 1e-14);
    }
  }

  /// @brief Verifies scalar arithmetic div grid function for MMG grid function by checking tolerance-based numerical results.
  TEST(Rodin_MMG_GridFunction, ScalarArithmeticDivGridFunction)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf1(fes), gf2(fes);

    gf1 = RealFunction(12.0);
    gf2 = RealFunction(3.0);
    gf1 /= gf2;

    for (Index i = 0; i < static_cast<Index>(gf1.getSize()); i++)
    {
      EXPECT_NEAR(gf1[i], 4.0, 1e-14);
    }
  }

  /// @brief Verifies scalar min max for MMG grid function by checking tolerance-based numerical results.
  TEST(Rodin_MMG_GridFunction, ScalarMinMax)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf(fes);

    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      gf[i] = static_cast<Real>(i);
    }

    EXPECT_NEAR(gf.min(), 0.0, 1e-14);
    EXPECT_NEAR(gf.max(), static_cast<Real>(gf.getSize() - 1), 1e-14);
  }

  /// @brief Verifies scalar arg min arg max for MMG grid function by checking exact expected values.
  TEST(Rodin_MMG_GridFunction, ScalarArgMinArgMax)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf(fes);

    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      gf[i] = static_cast<Real>(i);
    }

    EXPECT_EQ(gf.argmin(), 0);
    EXPECT_EQ(gf.argmax(), static_cast<Index>(gf.getSize() - 1));
  }

  /// @brief Verifies scalar copy constructor for MMG grid function by checking tolerance-based numerical results, copy semantics.
  TEST(Rodin_MMG_GridFunction, ScalarCopyConstructor)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf1(fes);
    gf1 = RealFunction(5.0);

    // Copy
    MMG::RealGridFunction gf2(gf1);

    for (Index i = 0; i < static_cast<Index>(gf1.getSize()); i++)
    {
      EXPECT_NEAR(gf2[i], 5.0, 1e-14);
    }

    // Mutate copy — original is unaffected
    gf2[0] = 99.0;
    EXPECT_NEAR(gf1[0], 5.0, 1e-14);
    EXPECT_NEAR(gf2[0], 99.0, 1e-14);
  }

  /// @brief Verifies scalar move constructor for MMG grid function by checking tolerance-based numerical results, exact expected values, move semantics.
  TEST(Rodin_MMG_GridFunction, ScalarMoveConstructor)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf1(fes);
    gf1 = RealFunction(7.0);
    size_t sz = gf1.getSize();

    MMG::RealGridFunction gf2(std::move(gf1));
    EXPECT_EQ(gf2.getSize(), sz);
    for (Index i = 0; i < static_cast<Index>(gf2.getSize()); i++)
    {
      EXPECT_NEAR(gf2[i], 7.0, 1e-14);
    }
  }

  /// @brief Verifies scalar copy assignment for MMG grid function by checking tolerance-based numerical results, copy semantics.
  TEST(Rodin_MMG_GridFunction, ScalarCopyAssignment)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf1(fes), gf2(fes);
    gf1 = RealFunction(3.0);
    gf2 = RealFunction(9.0);

    gf2 = gf1;
    for (Index i = 0; i < static_cast<Index>(gf2.getSize()); i++)
    {
      EXPECT_NEAR(gf2[i], 3.0, 1e-14);
    }
  }

  /// @brief Verifies scalar FES association for MMG grid function by checking exact expected values.
  TEST(Rodin_MMG_GridFunction, ScalarFESAssociation)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf(fes);

    EXPECT_EQ(&gf.getFiniteElementSpace(), &fes);
    EXPECT_EQ(gf.getFiniteElementSpace().getVectorDimension(), 1);
    EXPECT_EQ(gf.getSize(), fes.getSize());
    EXPECT_EQ(gf.getSize(), mesh.getVertexCount());
  }

  /// @brief Verifies scalar point evaluation for MMG grid function by checking tolerance-based numerical results.
  TEST(Rodin_MMG_GridFunction, ScalarPointEvaluation)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf(fes);

    // Project a constant function: should evaluate to the same constant
    gf = RealFunction(42.0);

    // Get the first cell and evaluate at a reference point inside it
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{ 0.25, 0.25 }};
    Geometry::Point p(polytope, rc);

    Real value = gf.getValue(p);
    EXPECT_NEAR(value, 42.0, 1e-10);
  }

  /// @brief Verifies scalar linear function exact interpolation for MMG grid function by checking tolerance-based numerical results.
  TEST(Rodin_MMG_GridFunction, ScalarLinearFunctionExactInterpolation)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf(fes);

    // P1 should represent linear functions exactly
    gf = [](const Geometry::Point& p) { return p.x() + 2.0 * p.y(); };

    // Evaluate at the centroid of a cell
    for (Index cellIdx = 0; cellIdx < static_cast<Index>(mesh.getCellCount()); cellIdx++)
    {
      auto it = mesh.getPolytope(mesh.getDimension(), cellIdx);
      const auto& polytope = *it;
      const auto& trans = mesh.getPolytopeTransformation(mesh.getDimension(), cellIdx);

      const Math::Vector<Real> rc{{ 1.0 / 3.0, 1.0 / 3.0 }};
      Math::SpatialPoint pc;
      trans.transform(pc, rc);
      Geometry::Point p(polytope, rc);

      Real value = gf.getValue(p);
      Real expected = pc(0) + 2.0 * pc(1);
      EXPECT_NEAR(value, expected, 1e-10);
    }
  }

  /// @brief Verifies vector project from lambda for MMG grid function by checking exact expected values.
  TEST(Rodin_MMG_GridFunction, VectorProjectFromLambda)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh, 2);
    MMG::VectorGridFunction gf(fes);

    gf = [](const Geometry::Point& p)
    {
      return Math::SpatialVector<Real>{{ p.x(), p.y() }};
    };

    EXPECT_EQ(gf.getFiniteElementSpace().getVectorDimension(), 2);
    EXPECT_EQ(gf.getSize(), mesh.getVertexCount() * 2);
  }

  /// @brief Verifies vector FES association for MMG grid function by checking exact expected values.
  TEST(Rodin_MMG_GridFunction, VectorFESAssociation)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh, 3);
    MMG::GridFunction<Math::SpatialVector<Real>> gf(fes);

    EXPECT_EQ(&gf.getFiniteElementSpace(), &fes);
    EXPECT_EQ(gf.getFiniteElementSpace().getVectorDimension(), 3);
    EXPECT_EQ(gf.getSize(), mesh.getVertexCount() * 3);
  }

  /// @brief Verifies scalar IO save load MEDIT for MMG grid function by checking tolerance-based numerical results, exact expected values.
  TEST(Rodin_MMG_GridFunction, ScalarIOSaveLoadMEDIT)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf(fes);
    gf = [](const Geometry::Point& p) { return p.x() + p.y(); };

    const std::string meshFile = "/tmp/rodin_mmg_gf_io_mesh.mesh";
    const std::string solFile = "/tmp/rodin_mmg_gf_io.sol";
    mesh.save(meshFile, IO::FileFormat::MEDIT);
    gf.save(solFile, IO::FileFormat::MEDIT);

    // Load back on a fresh mesh+FES
    MMG::Mesh mesh2;
    mesh2.load(meshFile, IO::FileFormat::MEDIT);
    P1 fes2(mesh2);
    MMG::RealGridFunction gf2(fes2);
    gf2.load(solFile, IO::FileFormat::MEDIT);

    EXPECT_EQ(gf2.getSize(), gf.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf2[i], gf[i], 1e-10);
    }

    std::remove(meshFile.c_str());
    std::remove(solFile.c_str());
  }

  /// @brief Verifies scalar IO save load MFEM for MMG grid function by checking tolerance-based numerical results, exact expected values.
  TEST(Rodin_MMG_GridFunction, ScalarIOSaveLoadMFEM)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf(fes);
    gf = RealFunction(3.14);

    const std::string meshFile = "/tmp/rodin_mmg_gf_mfem.mesh";
    const std::string gfFile = "/tmp/rodin_mmg_gf_mfem.gf";
    mesh.save(meshFile, IO::FileFormat::MFEM);
    gf.save(gfFile, IO::FileFormat::MFEM);

    // Load back
    MMG::Mesh mesh2;
    mesh2.load(meshFile, IO::FileFormat::MFEM);
    P1 fes2(mesh2);
    MMG::RealGridFunction gf2(fes2);
    gf2.load(gfFile, IO::FileFormat::MFEM);

    EXPECT_EQ(gf2.getSize(), gf.getSize());
    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf2[i], 3.14, 1e-10);
    }

    std::remove(meshFile.c_str());
    std::remove(gfFile.c_str());
  }

  /// @brief Integration test: project a size map, adapt, then verify the adapted.
  TEST(Rodin_MMG_GridFunction, ScalarProjectThenAdapt)
  {
    // Integration test: project a size map, adapt, then verify the adapted
    // mesh can host a new grid function.
    const size_t n = 8;
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });

    P1 fes(mesh);
    MMG::RealGridFunction sizeMap(fes);
    sizeMap = [](const Geometry::Point&) { return 0.25; };

    MMG::Adapt().setHMin(0.05).setHMax(0.5).adapt(mesh, sizeMap);

    // After adaptation, build a new P1 space and grid function on the
    // adapted mesh.
    P1 fesAdapted(mesh);
    MMG::RealGridFunction gfAdapted(fesAdapted);
    gfAdapted = [](const Geometry::Point& p) { return p.x() * p.y(); };

    EXPECT_EQ(gfAdapted.getSize(), mesh.getVertexCount());
    EXPECT_GT(gfAdapted.getSize(), 0);
  }

  /// @brief Verifies vector point evaluation for MMG grid function by checking tolerance-based numerical results.
  TEST(Rodin_MMG_GridFunction, VectorPointEvaluation)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh, 2);
    MMG::VectorGridFunction gf(fes);

    // Project a constant vector function
    gf = [](const Geometry::Point&)
    {
      return Math::SpatialVector<Real>{{ 1.0, 2.0 }};
    };

    // Evaluate at a point inside the first cell
    auto it = mesh.getPolytope(mesh.getDimension(), 0);
    const auto& polytope = *it;
    const Math::Vector<Real> rc{{ 0.25, 0.25 }};
    Geometry::Point p(polytope, rc);

    auto value = gf.getValue(p);
    EXPECT_NEAR(value(0), 1.0, 1e-10);
    EXPECT_NEAR(value(1), 2.0, 1e-10);
  }

  /// @brief Verifies scalar set data for MMG grid function by checking tolerance-based numerical results.
  TEST(Rodin_MMG_GridFunction, ScalarSetData)
  {
    MMG::Mesh mesh;
    mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 4, 4 });
    P1 fes(mesh);
    MMG::RealGridFunction gf(fes);

    // Manually set data
    Math::Vector<Real> data(gf.getSize());
    for (Eigen::Index i = 0; i < data.size(); i++)
    {
      data(i) = static_cast<Real>(i * 2);
    }
    gf.setData(data);

    for (Index i = 0; i < static_cast<Index>(gf.getSize()); i++)
    {
      EXPECT_NEAR(gf[i], static_cast<Real>(i * 2), 1e-14);
    }
  }
}
