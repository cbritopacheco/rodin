/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <initializer_list>
#include <limits>
#include <set>
#include <vector>

#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/LevelSetInterfaceGraph.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Unit
{
  static LocalMesh makeSingleTriangle()
  {
    return LocalMesh::Builder()
      .initialize(2)
      .nodes(3)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .finalize();
  }

  static LocalMesh makeTwoTriangleSquare()
  {
    return LocalMesh::Builder()
      .initialize(2)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .polytope(Polytope::Type::Triangle, {1, 3, 2})
      .finalize();
  }

  static LocalMesh makeFourTriangleFan()
  {
    return LocalMesh::Builder()
      .initialize(2)
      .nodes(5)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({1, 1})
      .vertex({0, 1})
      .vertex({0.5, 0.5})
      .polytope(Polytope::Type::Triangle, {0, 1, 4})
      .polytope(Polytope::Type::Triangle, {1, 2, 4})
      .polytope(Polytope::Type::Triangle, {2, 3, 4})
      .polytope(Polytope::Type::Triangle, {3, 0, 4})
      .finalize();
  }

  static void computeConnectivity(LocalMesh& mesh)
  {
    auto& conn = mesh.getConnectivity();
    conn.compute(2, 1);
    conn.compute(1, 0);
  }

  static Index findEdge(const LocalMesh& mesh, Index a, Index b)
  {
    const auto& conn = mesh.getConnectivity();
    for (Index e = 0; e < conn.getCount(1); ++e)
    {
      const auto& edge = conn.getPolytope(1, e);
      if ((edge(0) == a && edge(1) == b) ||
          (edge(0) == b && edge(1) == a))
        return e;
    }
    return std::numeric_limits<Index>::max();
  }

  static P1InterfaceGraph extract(LocalMesh& mesh, std::initializer_list<Real> values)
  {
    computeConnectivity(mesh);

    P1 space(mesh);
    GridFunction phi(space);

    Index i = 0;
    for (Real value : values)
      phi[i++] = value;

    return LevelSetInterfaceGraph(phi)
      .setSignTolerance(1e-12)
      .setInterfaceAttribute(42)
      .extract();
  }

  static bool hasParentEdge(const InterfaceEdge& edge, Index parent)
  {
    return edge.parentEdges[0] == parent || edge.parentEdges[1] == parent;
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, TriangleAllNegative)
  {
    auto mesh = makeSingleTriangle();
    const auto graph = extract(mesh, {-1, -2, -3});

    EXPECT_TRUE(graph.vertices.empty());
    EXPECT_TRUE(graph.edges.empty());
    EXPECT_TRUE(graph.degenerateCells.empty());
    EXPECT_TRUE(graph.invalidCells.empty());
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, TriangleAllPositive)
  {
    auto mesh = makeSingleTriangle();
    const auto graph = extract(mesh, {1, 2, 3});

    EXPECT_TRUE(graph.vertices.empty());
    EXPECT_TRUE(graph.edges.empty());
    EXPECT_TRUE(graph.degenerateCells.empty());
    EXPECT_TRUE(graph.invalidCells.empty());
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, OneNegativeTwoPositive)
  {
    auto mesh = makeSingleTriangle();
    const auto graph = extract(mesh, {-1, 1, 1});

    ASSERT_EQ(graph.vertices.size(), 2);
    ASSERT_EQ(graph.edges.size(), 1);

    const Index e01 = findEdge(mesh, 0, 1);
    const Index e20 = findEdge(mesh, 2, 0);
    const auto& edge = graph.edges.front();

    EXPECT_EQ(edge.parentCell, 0);
    EXPECT_TRUE(hasParentEdge(edge, e01));
    EXPECT_TRUE(hasParentEdge(edge, e20));
    EXPECT_EQ(edge.interfaceAttribute, 42);
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, TwoNegativeOnePositive)
  {
    auto mesh = makeSingleTriangle();
    const auto graph = extract(mesh, {-1, -1, 1});

    ASSERT_EQ(graph.vertices.size(), 2);
    ASSERT_EQ(graph.edges.size(), 1);

    const Index e12 = findEdge(mesh, 1, 2);
    const Index e20 = findEdge(mesh, 2, 0);

    EXPECT_TRUE(hasParentEdge(graph.edges.front(), e12));
    EXPECT_TRUE(hasParentEdge(graph.edges.front(), e20));
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, OneZeroVertexOppositeSigns)
  {
    auto mesh = makeSingleTriangle();
    const auto graph = extract(mesh, {0, -1, 1});

    ASSERT_EQ(graph.vertices.size(), 2);
    ASSERT_EQ(graph.edges.size(), 1);

    const auto zero = std::find_if(
      graph.vertices.begin(), graph.vertices.end(),
      [](const InterfaceVertex& vertex)
      {
        return vertex.snappedToOriginalVertex &&
               vertex.originalVertex &&
               *vertex.originalVertex == 0;
      });

    ASSERT_NE(zero, graph.vertices.end());
    EXPECT_TRUE(hasParentEdge(graph.edges.front(), findEdge(mesh, 1, 2)));
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, ZeroEdge)
  {
    auto mesh = makeSingleTriangle();
    const auto graph = extract(mesh, {0, 0, 1});

    ASSERT_EQ(graph.vertices.size(), 2);
    ASSERT_EQ(graph.edges.size(), 1);

    const Index e01 = findEdge(mesh, 0, 1);
    EXPECT_EQ(graph.edges.front().parentEdges[0], e01);
    EXPECT_EQ(graph.edges.front().parentEdges[1], e01);

    std::set<Index> originals;
    for (const auto& vertex : graph.vertices)
    {
      ASSERT_TRUE(vertex.snappedToOriginalVertex);
      ASSERT_TRUE(vertex.originalVertex);
      originals.insert(*vertex.originalVertex);
    }

    EXPECT_EQ(originals, (std::set<Index>{0, 1}));
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, AllZeroTriangleIsDegenerate)
  {
    auto mesh = makeSingleTriangle();
    const auto graph = extract(mesh, {0, 0, 0});

    EXPECT_TRUE(graph.vertices.empty());
    EXPECT_TRUE(graph.edges.empty());
    ASSERT_EQ(graph.degenerateCells.size(), 1);
    EXPECT_EQ(graph.degenerateCells.front(), 0);
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, SharedCutEdgeProducesOneVertex)
  {
    auto mesh = makeTwoTriangleSquare();
    const auto graph = extract(mesh, {1, -1, 1, -1});

    const Index shared = findEdge(mesh, 1, 2);
    size_t sharedVertices = 0;
    for (const auto& vertex : graph.vertices)
    {
      if (!vertex.snappedToOriginalVertex && vertex.parentEdge == shared)
        sharedVertices++;
    }

    EXPECT_EQ(sharedVertices, 1);
    EXPECT_EQ(graph.vertices.size(), 3);
    EXPECT_EQ(graph.edges.size(), 2);

    size_t incidentGraphEdges = 0;
    for (const auto& edge : graph.edges)
    {
      if (hasParentEdge(edge, shared))
        incidentGraphEdges++;
    }
    EXPECT_EQ(incidentGraphEdges, 2);
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, ClosedLoopOnTriangulatedSquare)
  {
    auto mesh = makeFourTriangleFan();
    const auto graph = extract(mesh, {1, 1, 1, 1, -1});

    EXPECT_EQ(graph.vertices.size(), 4);
    EXPECT_EQ(graph.edges.size(), 4);
    ASSERT_EQ(graph.loops.size(), 1);
    EXPECT_EQ(graph.loops.front().vertices.size(), 4);
    EXPECT_EQ(graph.loops.front().edges.size(), 4);
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, PreservesCellAttributeProvenance)
  {
    auto mesh = makeSingleTriangle();
    mesh.setAttribute({2, 0}, 7);

    const auto graph = extract(mesh, {-1, 1, 1});

    ASSERT_EQ(graph.edges.size(), 1);
    ASSERT_TRUE(graph.edges.front().negativeCellAttribute);
    ASSERT_TRUE(graph.edges.front().positiveCellAttribute);
    EXPECT_EQ(*graph.edges.front().negativeCellAttribute, 7);
    EXPECT_EQ(*graph.edges.front().positiveCellAttribute, 7);
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, NoDuplicatedVerticesOnSharedMeshEdges)
  {
    auto mesh = makeTwoTriangleSquare();
    const auto graph = extract(mesh, {1, -1, 1, -1});

    std::set<Index> parentEdges;
    for (const auto& vertex : graph.vertices)
    {
      ASSERT_TRUE(parentEdges.insert(vertex.parentEdge).second)
        << "duplicated interface vertex on mesh edge " << vertex.parentEdge;
    }
  }
}
