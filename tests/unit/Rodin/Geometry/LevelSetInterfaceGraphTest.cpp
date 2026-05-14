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
#include <string>
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

  static InterfaceGraph extract(LocalMesh& mesh, std::initializer_list<Real> values)
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
    return std::any_of(
      edge.provenance.begin(), edge.provenance.end(),
      [parent](const InterfaceEdgeProvenance& provenance)
      {
        return provenance.parentEdges[0] == parent ||
               provenance.parentEdges[1] == parent;
      });
  }

  static std::set<std::pair<Index, Index>> unorderedEdgePairs(const InterfaceGraph& graph)
  {
    std::set<std::pair<Index, Index>> res;
    for (const auto& edge : graph.edges)
      res.insert(std::minmax(edge.v0, edge.v1));
    return res;
  }

  static void expectRawPSLGMatchesGraph(const InterfaceGraph& graph)
  {
    const auto pslg = InterfaceGraphPSLGBuilder().build(graph);

    EXPECT_EQ(pslg.vertices.size(), graph.vertices.size());
    ASSERT_EQ(pslg.segments.size(), graph.edges.size());

    for (Index i = 0; i < pslg.segments.size(); ++i)
    {
      const auto& segment = pslg.segments[i];
      const auto& edge = graph.edges[i];

      EXPECT_EQ(segment.v0, edge.v0);
      EXPECT_EQ(segment.v1, edge.v1);
      EXPECT_EQ(segment.attribute, edge.interfaceAttribute);
      ASSERT_EQ(segment.sourceInterfaceEdges.size(), 1);
      EXPECT_EQ(segment.sourceInterfaceEdges.front(), i);
    }
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, TriangleAllNegative)
  {
    auto mesh = makeSingleTriangle();
    const auto graph = extract(mesh, {-1, -2, -3});

    EXPECT_TRUE(graph.vertices.empty());
    EXPECT_TRUE(graph.edges.empty());
    EXPECT_TRUE(graph.chains.empty());
    EXPECT_TRUE(graph.degenerateCells.empty());
    EXPECT_TRUE(graph.invalidCells.empty());
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, TriangleAllPositive)
  {
    auto mesh = makeSingleTriangle();
    const auto graph = extract(mesh, {1, 2, 3});

    EXPECT_TRUE(graph.vertices.empty());
    EXPECT_TRUE(graph.edges.empty());
    EXPECT_TRUE(graph.chains.empty());
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

    ASSERT_EQ(edge.provenance.size(), 1);
    EXPECT_EQ(edge.provenance.front().parentCell, 0);
    EXPECT_TRUE(hasParentEdge(edge, e01));
    EXPECT_TRUE(hasParentEdge(edge, e20));
    ASSERT_TRUE(edge.interfaceAttribute);
    EXPECT_EQ(*edge.interfaceAttribute, 42);
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, DefaultInterfaceAttributeIsUnset)
  {
    auto mesh = makeSingleTriangle();
    computeConnectivity(mesh);

    P1 space(mesh);
    GridFunction phi(space);
    phi[0] = -1;
    phi[1] = 1;
    phi[2] = 1;

    const auto graph = LevelSetInterfaceGraph(phi).extract();

    ASSERT_EQ(graph.edges.size(), 1);
    EXPECT_FALSE(graph.edges.front().interfaceAttribute);
    const auto pslg = InterfaceGraphPSLGBuilder().build(graph);
    ASSERT_EQ(pslg.segments.size(), 1);
    EXPECT_FALSE(pslg.segments.front().attribute);
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, RequiresCellToEdgeConnectivity)
  {
    auto mesh = makeSingleTriangle();

    P1 space(mesh);
    GridFunction phi(space);
    phi[0] = -1;
    phi[1] = 1;
    phi[2] = 1;

    try
    {
      (void) LevelSetInterfaceGraph(phi).extract();
      FAIL() << "Expected missing cell-to-edge connectivity to throw.";
    }
    catch (const Alert::Exception& e)
    {
      EXPECT_NE(
        std::string(e.what()).find("has not been computed"),
        std::string::npos);
    }
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, RequiresEdgeToVertexConnectivity)
  {
    auto mesh = makeSingleTriangle();
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().clear(1, 0);

    P1 space(mesh);
    GridFunction phi(space);
    phi[0] = -1;
    phi[1] = 1;
    phi[2] = 1;

    try
    {
      (void) LevelSetInterfaceGraph(phi).extract();
      FAIL() << "Expected missing edge-to-vertex connectivity to throw.";
    }
    catch (const Alert::Exception& e)
    {
      EXPECT_NE(
        std::string(e.what()).find("has not been computed"),
        std::string::npos);
    }
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
        return vertex.kind == InterfaceVertexKind::OriginalVertex &&
               vertex.originalVertex &&
               *vertex.originalVertex == 0;
      });

    ASSERT_NE(zero, graph.vertices.end());
    EXPECT_FALSE(zero->parentEdge);
    EXPECT_FALSE(zero->t);
    EXPECT_TRUE(hasParentEdge(graph.edges.front(), findEdge(mesh, 1, 2)));
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, ZeroEdge)
  {
    auto mesh = makeSingleTriangle();
    const auto graph = extract(mesh, {0, 0, 1});

    ASSERT_EQ(graph.vertices.size(), 2);
    ASSERT_EQ(graph.edges.size(), 1);
    ASSERT_EQ(graph.edges.front().provenance.size(), 1);

    const Index e01 = findEdge(mesh, 0, 1);
    EXPECT_EQ(graph.edges.front().provenance.front().parentEdges[0], e01);
    EXPECT_EQ(graph.edges.front().provenance.front().parentEdges[1], e01);

    std::set<Index> originals;
    for (const auto& vertex : graph.vertices)
    {
      EXPECT_EQ(vertex.kind, InterfaceVertexKind::OriginalVertex);
      ASSERT_TRUE(vertex.originalVertex);
      EXPECT_FALSE(vertex.parentEdge);
      EXPECT_FALSE(vertex.t);
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

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, IsolatedZeroPositiveNeighborsCreatesNoInterface)
  {
    auto mesh = makeSingleTriangle();
    const auto graph = extract(mesh, {0, 1, 1});

    EXPECT_TRUE(graph.vertices.empty());
    EXPECT_TRUE(graph.edges.empty());
    EXPECT_TRUE(graph.chains.empty());
    EXPECT_TRUE(graph.degenerateCells.empty());
    EXPECT_TRUE(graph.invalidCells.empty());
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, IsolatedZeroNegativeNeighborsCreatesNoInterface)
  {
    auto mesh = makeSingleTriangle();
    const auto graph = extract(mesh, {0, -1, -1});

    EXPECT_TRUE(graph.vertices.empty());
    EXPECT_TRUE(graph.edges.empty());
    EXPECT_TRUE(graph.chains.empty());
    EXPECT_TRUE(graph.degenerateCells.empty());
    EXPECT_TRUE(graph.invalidCells.empty());
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, SharedCutEdgeProducesOneVertex)
  {
    auto mesh = makeTwoTriangleSquare();
    const auto graph = extract(mesh, {1, -1, 1, -1});

    const Index shared = findEdge(mesh, 1, 2);
    size_t sharedVertices = 0;
    for (const auto& vertex : graph.vertices)
    {
      if (vertex.kind == InterfaceVertexKind::EdgeCut &&
          vertex.parentEdge &&
          *vertex.parentEdge == shared)
      {
        sharedVertices++;
      }
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
    ASSERT_EQ(graph.chains.size(), 1);
    EXPECT_TRUE(graph.chains.front().closed);
    EXPECT_EQ(graph.chains.front().vertices.size(), 4);
    EXPECT_EQ(graph.chains.front().edges.size(), 4);
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, PreservesCellAttributeProvenance)
  {
    auto mesh = makeSingleTriangle();
    mesh.setAttribute({2, 0}, 7);

    const auto graph = extract(mesh, {-1, 1, 1});

    ASSERT_EQ(graph.edges.size(), 1);
    ASSERT_EQ(graph.edges.front().provenance.size(), 1);
    ASSERT_TRUE(graph.edges.front().provenance.front().parentCellAttribute);
    EXPECT_EQ(*graph.edges.front().provenance.front().parentCellAttribute, 7);
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, NoDuplicatedVerticesOnSharedMeshEdges)
  {
    auto mesh = makeTwoTriangleSquare();
    const auto graph = extract(mesh, {1, -1, 1, -1});

    std::set<Index> parentEdges;
    for (const auto& vertex : graph.vertices)
    {
      ASSERT_EQ(vertex.kind, InterfaceVertexKind::EdgeCut);
      ASSERT_TRUE(vertex.parentEdge);
      ASSERT_TRUE(parentEdges.insert(*vertex.parentEdge).second)
        << "duplicated interface vertex on mesh edge " << *vertex.parentEdge;
    }
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, SharedZeroEdgeCreatesOneGraphEdge)
  {
    auto mesh = makeTwoTriangleSquare();
    const auto graph = extract(mesh, {1, 0, 0, 1});

    ASSERT_EQ(graph.vertices.size(), 2);
    ASSERT_EQ(graph.edges.size(), 1);
    ASSERT_EQ(graph.edges.front().provenance.size(), 2);
    EXPECT_TRUE(hasParentEdge(graph.edges.front(), findEdge(mesh, 1, 2)));
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, ExtractsOpenChainAcrossSquare)
  {
    auto mesh = makeTwoTriangleSquare();
    const auto graph = extract(mesh, {-1, 1, -1, 1});

    ASSERT_EQ(graph.vertices.size(), 3);
    ASSERT_EQ(graph.edges.size(), 2);
    ASSERT_EQ(graph.chains.size(), 1);

    const auto& chain = graph.chains.front();
    EXPECT_FALSE(chain.closed);
    EXPECT_EQ(chain.vertices.size(), 3);
    EXPECT_EQ(chain.edges.size(), 2);
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, InvalidPhiValueMarksCell)
  {
    auto mesh = makeSingleTriangle();
    const auto graph = extract(mesh, {
      std::numeric_limits<Real>::quiet_NaN(), -1, 1});

    EXPECT_TRUE(graph.vertices.empty());
    EXPECT_TRUE(graph.edges.empty());
    ASSERT_EQ(graph.invalidCells.size(), 1);
    EXPECT_EQ(graph.invalidCells.front(), 0);
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, EdgeCutStoresExactInterpolation)
  {
    auto mesh = makeSingleTriangle();
    computeConnectivity(mesh);

    P1 space(mesh);
    GridFunction phi(space);
    phi[0] = -1;
    phi[1] = 3;
    phi[2] = 3;

    const auto graph = LevelSetInterfaceGraph(phi)
      .setSignTolerance(1e-12)
      .extract();

    ASSERT_EQ(graph.vertices.size(), 2);
    for (const auto& vertex : graph.vertices)
    {
      ASSERT_EQ(vertex.kind, InterfaceVertexKind::EdgeCut);
      ASSERT_TRUE(vertex.parentEdge);
      ASSERT_TRUE(vertex.t);

      const auto& edge = mesh.getConnectivity().getPolytope(1, *vertex.parentEdge);
      const Index a = edge(0);
      const Index b = edge(1);
      const Real expectedT = phi[a] / (phi[a] - phi[b]);
      const auto expectedX =
        (Real(1) - expectedT) * mesh.getVertexCoordinates(a) +
        expectedT * mesh.getVertexCoordinates(b);

      EXPECT_NEAR(*vertex.t, expectedT, 1e-14);
      EXPECT_NEAR(vertex.x[0], expectedX[0], 1e-14);
      EXPECT_NEAR(vertex.x[1], expectedX[1], 1e-14);
    }
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, NoDuplicatedGraphEdges)
  {
    auto mesh = makeTwoTriangleSquare();
    const auto graph = extract(mesh, {1, 0, 0, -1});

    EXPECT_EQ(unorderedEdgePairs(graph).size(), graph.edges.size());
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, BuildsRawPSLGForOpenChain)
  {
    auto mesh = makeTwoTriangleSquare();
    const auto graph = extract(mesh, {-1, 1, -1, 1});

    ASSERT_EQ(graph.chains.size(), 1);
    ASSERT_FALSE(graph.chains.front().closed);
    expectRawPSLGMatchesGraph(graph);
  }

  TEST(Rodin_Geometry_LevelSetInterfaceGraph, BuildsRawPSLGForClosedChain)
  {
    auto mesh = makeFourTriangleFan();
    const auto graph = extract(mesh, {1, 1, 1, 1, -1});

    ASSERT_EQ(graph.chains.size(), 1);
    ASSERT_TRUE(graph.chains.front().closed);
    expectRawPSLGMatchesGraph(graph);
  }
}
