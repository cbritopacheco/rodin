/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <limits>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/LevelSetDiscretizerTriangles.h>
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

  static void computeConnectivity(LocalMesh& mesh)
  {
    auto& conn = mesh.getConnectivity();
    conn.compute(2, 1);
    conn.compute(1, 0);
  }

  static LevelSetDiscretizerTrianglesResult discretize(
      LocalMesh& mesh,
      std::initializer_list<Real> values,
      Optional<Attribute> interfaceAttribute = Optional<Attribute>(42),
      Optional<Attribute> negativeAttribute = {},
      Optional<Attribute> positiveAttribute = {},
      Real diagnosticTolerance = 1e-10)
  {
    computeConnectivity(mesh);

    P1 space(mesh);
    GridFunction phi(space);

    Index i = 0;
    for (Real value : values)
      phi[i++] = value;

    return LevelSetDiscretizerTriangles(phi)
      .setSignTolerance(1e-12)
      .setDiagnosticTolerance(diagnosticTolerance)
      .setInterfaceAttribute(interfaceAttribute)
      .setNegativeCellAttribute(negativeAttribute)
      .setPositiveCellAttribute(positiveAttribute)
      .discretize();
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

  static Index sideCount(
      const LevelSetDiscretizerTrianglesReport& report,
      LevelSetSide side)
  {
    return static_cast<Index>(std::count_if(
      report.cellProvenance.begin(), report.cellProvenance.end(),
      [side](const OutputCellProvenance& provenance)
      {
        return provenance.side == side;
      }));
  }

  static Real signedTriangleArea(
      const Math::SpatialPoint& a,
      const Math::SpatialPoint& b,
      const Math::SpatialPoint& c)
  {
    return Real(0.5)
      * ((b[0] - a[0]) * (c[1] - a[1])
       - (b[1] - a[1]) * (c[0] - a[0]));
  }

  static Real triangleQuality(
      const Math::SpatialPoint& a,
      const Math::SpatialPoint& b,
      const Math::SpatialPoint& c)
  {
    const Real area = std::abs(signedTriangleArea(a, b, c));
    const Real l0 = (b - a).squaredNorm();
    const Real l1 = (c - b).squaredNorm();
    const Real l2 = (a - c).squaredNorm();
    const Real denom = l0 + l1 + l2;
    if (denom <= Real(0))
      return Real(0);
    return Real(4) * std::sqrt(Real(3)) * area / denom;
  }

  static Real minTriangleQuality(const LocalMesh& mesh)
  {
    Real qmin = std::numeric_limits<Real>::infinity();
    const auto& conn = mesh.getConnectivity();
    for (Index c = 0; c < mesh.getCellCount(); ++c)
    {
      if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
        continue;
      const auto& cell = conn.getPolytope(2, c);
      qmin = std::min(
          qmin,
          triangleQuality(
            mesh.getVertexCoordinates(cell(0)),
            mesh.getVertexCoordinates(cell(1)),
            mesh.getVertexCoordinates(cell(2))));
    }
    return std::isfinite(qmin) ? qmin : Real(0);
  }

  static void expectEveryGraphEdgeAppearsAsOutputEdge(
      const LevelSetDiscretizerTrianglesResult& result)
  {
    ASSERT_EQ(
      result.report.interfaceEdgeProvenance.size(),
      result.interfaceGraph.edges.size());

    for (Index ge = 0; ge < result.interfaceGraph.edges.size(); ++ge)
    {
      const auto& graphEdge = result.interfaceGraph.edges[ge];
      const auto out0 =
        result.report.graphVertexToOutputVertex.at(graphEdge.v0);
      const auto out1 =
        result.report.graphVertexToOutputVertex.at(graphEdge.v1);

      const Index oe = findEdge(result.mesh, out0, out1);
      ASSERT_NE(oe, std::numeric_limits<Index>::max());

      const auto provenance = result.report.interfaceEdgeProvenance.find(oe);
      ASSERT_NE(provenance, result.report.interfaceEdgeProvenance.end());
      EXPECT_EQ(provenance->second.sourceInterfaceGraphEdge, ge);
    }
  }

  static LevelSetSign classify(Real value)
  {
    if (value < -1e-12)
      return LevelSetSign::Negative;
    if (value > 1e-12)
      return LevelSetSign::Positive;
    return LevelSetSign::Zero;
  }

  static void expectNoOutputCellCrossesInterface(
      const LevelSetDiscretizerTrianglesResult& result,
      const std::vector<Real>& phi)
  {
    const auto& conn = result.mesh.getConnectivity();
    ASSERT_EQ(
      result.report.cellProvenance.size(),
      static_cast<size_t>(result.mesh.getCellCount()));

    auto vertexSign = [&](Index outputVertex)
    {
      const auto& origin = result.report.vertexOrigins[outputVertex];
      if (origin.kind == OutputVertexOriginKind::InterfaceGraphVertex)
        return LevelSetSign::Zero;
      if (origin.originalVertex)
        return classify(phi[*origin.originalVertex]);
      return LevelSetSign::Invalid;
    };

    for (Index c = 0; c < result.mesh.getCellCount(); ++c)
    {
      const auto side = result.report.cellProvenance[c].side;
      const auto& cell = conn.getPolytope(2, c);
      for (std::uint8_t i = 0; i < 3; ++i)
      {
        const auto sign = vertexSign(cell(i));
        if (side == LevelSetSide::Negative)
          EXPECT_NE(sign, LevelSetSign::Positive);
        if (side == LevelSetSide::Positive)
          EXPECT_NE(sign, LevelSetSign::Negative);
      }
    }
  }

  static Math::SpatialPoint reconstructFromParent(
      const LocalMesh& parent,
      Index parentCell,
      const std::array<Real, 3>& barycentric)
  {
    const auto& cell = parent.getConnectivity().getPolytope(2, parentCell);
    Math::SpatialPoint x{0, 0};
    for (std::uint8_t i = 0; i < 3; ++i)
    {
      const auto p = parent.getVertexCoordinates(cell(i));
      x[0] += barycentric[i] * p[0];
      x[1] += barycentric[i] * p[1];
    }
    return x;
  }

  static void expectCellReferencesReconstructVertices(
      const LocalMesh& parent,
      const LevelSetDiscretizerTrianglesResult& result)
  {
    const auto& conn = result.mesh.getConnectivity();
    ASSERT_EQ(result.report.cellReferences.size(), result.mesh.getCellCount());
    ASSERT_EQ(result.report.cellProvenance.size(), result.mesh.getCellCount());

    for (Index c = 0; c < result.mesh.getCellCount(); ++c)
    {
      const auto& reference = result.report.cellReferences[c];
      EXPECT_EQ(reference.parentCell, result.report.cellProvenance[c].parentCell);
      ASSERT_LT(reference.parentCell, parent.getCellCount());

      const auto& outCell = conn.getPolytope(2, c);
      for (std::uint8_t local = 0; local < 3; ++local)
      {
        const auto& bary = reference.vertexBarycentric[local];
        EXPECT_NEAR(bary[0] + bary[1] + bary[2], 1, 1e-12);
        const auto reconstructed =
          reconstructFromParent(parent, reference.parentCell, bary);
        const auto actual = result.mesh.getVertexCoordinates(outCell(local));
        EXPECT_NEAR((reconstructed - actual).norm(), 0, 1e-12);
      }
    }
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, AllNegativeTriangle)
  {
    auto mesh = makeSingleTriangle();
    const auto result = discretize(mesh, {-1, -2, -3});

    EXPECT_EQ(result.mesh.getCellCount(), 1);
    EXPECT_EQ(sideCount(result.report, LevelSetSide::Negative), 1);
    EXPECT_TRUE(result.interfaceGraph.edges.empty());
    EXPECT_TRUE(result.report.interfaceEdgeProvenance.empty());
    expectNoOutputCellCrossesInterface(result, {-1, -2, -3});
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, AllPositiveTriangle)
  {
    auto mesh = makeSingleTriangle();
    const auto result = discretize(mesh, {1, 2, 3});

    EXPECT_EQ(result.mesh.getCellCount(), 1);
    EXPECT_EQ(sideCount(result.report, LevelSetSide::Positive), 1);
    EXPECT_TRUE(result.interfaceGraph.edges.empty());
    EXPECT_TRUE(result.report.interfaceEdgeProvenance.empty());
    expectNoOutputCellCrossesInterface(result, {1, 2, 3});
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, OneNegativeTwoPositive)
  {
    auto mesh = makeSingleTriangle();
    const auto result = discretize(mesh, {-1, 1, 1});

    EXPECT_EQ(result.mesh.getCellCount(), 3);
    EXPECT_EQ(sideCount(result.report, LevelSetSide::Negative), 1);
    EXPECT_EQ(sideCount(result.report, LevelSetSide::Positive), 2);
    EXPECT_EQ(result.interfaceGraph.edges.size(), 1);
    expectEveryGraphEdgeAppearsAsOutputEdge(result);
    expectNoOutputCellCrossesInterface(result, {-1, 1, 1});
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, TwoNegativeOnePositive)
  {
    auto mesh = makeSingleTriangle();
    const auto result = discretize(mesh, {-1, -1, 1});

    EXPECT_EQ(result.mesh.getCellCount(), 3);
    EXPECT_EQ(sideCount(result.report, LevelSetSide::Negative), 2);
    EXPECT_EQ(sideCount(result.report, LevelSetSide::Positive), 1);
    EXPECT_EQ(result.interfaceGraph.edges.size(), 1);
    expectEveryGraphEdgeAppearsAsOutputEdge(result);
    expectNoOutputCellCrossesInterface(result, {-1, -1, 1});
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, OneZeroVertexWithOppositeSigns)
  {
    auto mesh = makeSingleTriangle();
    const auto result = discretize(mesh, {0, -1, 1});

    EXPECT_EQ(result.mesh.getCellCount(), 2);
    EXPECT_EQ(sideCount(result.report, LevelSetSide::Negative), 1);
    EXPECT_EQ(sideCount(result.report, LevelSetSide::Positive), 1);
    EXPECT_EQ(result.interfaceGraph.edges.size(), 1);
    expectEveryGraphEdgeAppearsAsOutputEdge(result);
    expectNoOutputCellCrossesInterface(result, {0, -1, 1});
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, ZeroEdgeIsPreservedAndMarked)
  {
    auto mesh = makeSingleTriangle();
    const auto result = discretize(mesh, {0, 0, 1}, 17);

    ASSERT_EQ(result.interfaceGraph.edges.size(), 1);
    ASSERT_EQ(result.report.interfaceEdgeProvenance.size(), 1);
    EXPECT_EQ(result.mesh.getCellCount(), 1);
    EXPECT_EQ(sideCount(result.report, LevelSetSide::Positive), 1);

    const Index out0 = result.report.originalVertexToOutputVertex.at(0);
    const Index out1 = result.report.originalVertexToOutputVertex.at(1);
    const Index edge = findEdge(result.mesh, out0, out1);
    ASSERT_NE(edge, std::numeric_limits<Index>::max());

    const auto attr = result.mesh.getAttribute(1, edge);
    ASSERT_TRUE(attr);
    EXPECT_EQ(*attr, 17);
    expectEveryGraphEdgeAppearsAsOutputEdge(result);
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, AllZeroTriangleIsReported)
  {
    auto mesh = makeSingleTriangle();
    const auto result = discretize(mesh, {0, 0, 0});

    EXPECT_EQ(result.mesh.getCellCount(), 1);
    EXPECT_EQ(result.report.degenerateCellCount, 1);
    ASSERT_EQ(result.report.cellProvenance.size(), 1);
    EXPECT_EQ(result.report.cellProvenance.front().side, LevelSetSide::Degenerate);
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, DegenerateCellAttributeCanBeCustomized)
  {
    auto mesh = makeSingleTriangle();
    computeConnectivity(mesh);
    mesh.setAttribute({2, 0}, Attribute(55));

    P1 space(mesh);
    GridFunction phi(space);
    phi[0] = 0;
    phi[1] = 0;
    phi[2] = 0;

    const auto result = LevelSetDiscretizerTriangles(phi)
      .setDegenerateCellAttribute(Attribute(99))
      .discretize();

    ASSERT_EQ(result.mesh.getCellCount(), 1);
    ASSERT_EQ(result.report.cellProvenance.front().side,
        LevelSetSide::Degenerate);
    const auto attr = result.mesh.getCellAttribute(0);
    ASSERT_TRUE(attr);
    EXPECT_EQ(*attr, Attribute(99));
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, SharedCutEdgeReusesOutputVertex)
  {
    auto mesh = makeTwoTriangleSquare();
    const auto result = discretize(mesh, {1, -1, 1, -1});

    EXPECT_EQ(result.interfaceGraph.vertices.size(), 3);
    EXPECT_EQ(result.mesh.getVertexCount(), 7);
    EXPECT_EQ(result.report.graphVertexToOutputVertex.size(), 3);

    std::set<Index> outputCutVertices;
    for (const auto& entry : result.report.graphVertexToOutputVertex)
      outputCutVertices.insert(entry.second);
    EXPECT_EQ(outputCutVertices.size(), 3);
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, SquareLineCutHasConnectedInterface)
  {
    auto mesh = makeTwoTriangleSquare();
    const auto result = discretize(mesh, {-1, 1, -1, 1});

    ASSERT_EQ(result.interfaceGraph.edges.size(), 2);
    ASSERT_EQ(result.interfaceGraph.chains.size(), 1);
    EXPECT_FALSE(result.interfaceGraph.chains.front().closed);
    expectEveryGraphEdgeAppearsAsOutputEdge(result);
    expectNoOutputCellCrossesInterface(result, {-1, 1, -1, 1});
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, OutputCellParentProvenance)
  {
    auto mesh = makeSingleTriangle();
    const auto result = discretize(mesh, {-1, 1, 1});

    ASSERT_EQ(result.report.cellProvenance.size(), 3);
    ASSERT_EQ(result.report.cellReferences.size(), 3);
    for (const auto& provenance : result.report.cellProvenance)
      EXPECT_EQ(provenance.parentCell, 0);
    expectCellReferencesReconstructVertices(mesh, result);
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, OutputCellReferencesTransferLinearP1Data)
  {
    auto mesh = makeTwoTriangleSquare();
    const std::vector<Real> nodalPhi = {-1, 1, -1, 1};
    const auto result = discretize(
        mesh,
        {nodalPhi[0], nodalPhi[1], nodalPhi[2], nodalPhi[3]});

    expectCellReferencesReconstructVertices(mesh, result);

    const auto& parentConn = mesh.getConnectivity();
    for (Index c = 0; c < result.mesh.getCellCount(); ++c)
    {
      const auto& reference = result.report.cellReferences[c];
      const auto& parentCell =
        parentConn.getPolytope(2, reference.parentCell);
      for (std::uint8_t local = 0; local < 3; ++local)
      {
        const auto& bary = reference.vertexBarycentric[local];
        Real transferred = 0;
        for (std::uint8_t i = 0; i < 3; ++i)
          transferred += bary[i] * nodalPhi[parentCell(i)];

        const auto x = reconstructFromParent(
            mesh, reference.parentCell, bary);
        const Real exact = Real(2) * x[0] - Real(1);
        EXPECT_NEAR(transferred, exact, 1e-12);
      }
    }
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, InterfaceEdgeProvenance)
  {
    auto mesh = makeSingleTriangle();
    const auto result = discretize(mesh, {-1, 1, 1});

    ASSERT_EQ(result.report.interfaceEdgeProvenance.size(), 1);
    const auto& provenance = result.report.interfaceEdgeProvenance.begin()->second;
    EXPECT_EQ(provenance.sourceInterfaceGraphEdge, 0);
    ASSERT_TRUE(provenance.parentCell);
    EXPECT_EQ(*provenance.parentCell, 0);
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, SideAttributesAreAssigned)
  {
    auto mesh = makeSingleTriangle();
    const auto result = discretize(mesh, {-1, 1, 1}, {}, 7, 8);

    const auto& conn = result.mesh.getConnectivity();
    for (Index c = 0; c < result.mesh.getCellCount(); ++c)
    {
      const auto attr = result.mesh.getAttribute(2, c);
      ASSERT_TRUE(attr);
      if (result.report.cellProvenance[c].side == LevelSetSide::Negative)
        EXPECT_EQ(*attr, 7);
      if (result.report.cellProvenance[c].side == LevelSetSide::Positive)
        EXPECT_EQ(*attr, 8);
      EXPECT_EQ(conn.getGeometry(2, c), Polytope::Type::Triangle);
    }
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, ParentCellAttributeIsPreserved)
  {
    auto mesh = makeSingleTriangle();
    mesh.setAttribute({2, 0}, 55);

    const auto result = discretize(mesh, {-1, 1, 1}, {});

    ASSERT_EQ(result.mesh.getCellCount(), 3);
    for (Index c = 0; c < result.mesh.getCellCount(); ++c)
    {
      const auto attr = result.mesh.getAttribute(2, c);
      ASSERT_TRUE(attr);
      EXPECT_EQ(*attr, 55);
    }
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, OriginalEdgeAttributeIsPreserved)
  {
    auto mesh = makeSingleTriangle();
    computeConnectivity(mesh);
    const Index edge01 = findEdge(mesh, 0, 1);
    mesh.setAttribute({1, edge01}, 66);

    P1 space(mesh);
    GridFunction phi(space);
    phi[0] = -1;
    phi[1] = -1;
    phi[2] = -1;

    const auto result = LevelSetDiscretizerTriangles(phi).discretize();
    const Index outEdge = findEdge(result.mesh, 0, 1);
    ASSERT_NE(outEdge, std::numeric_limits<Index>::max());

    const auto attr = result.mesh.getAttribute(1, outEdge);
    ASSERT_TRUE(attr);
    EXPECT_EQ(*attr, 66);
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, InterfaceAttributeCanBeUnset)
  {
    auto mesh = makeSingleTriangle();
    const auto result = discretize(mesh, {-1, 1, 1}, {});

    ASSERT_EQ(result.interfaceGraph.edges.size(), 1);
    const auto& graphEdge = result.interfaceGraph.edges.front();
    const auto out0 = result.report.graphVertexToOutputVertex.at(graphEdge.v0);
    const auto out1 = result.report.graphVertexToOutputVertex.at(graphEdge.v1);
    const Index edge = findEdge(result.mesh, out0, out1);
    ASSERT_NE(edge, std::numeric_limits<Index>::max());
    EXPECT_FALSE(result.mesh.getAttribute(1, edge));
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, PathologicalCutIsReported)
  {
    auto mesh = makeSingleTriangle();
    const auto result = discretize(mesh, {-0.01, 1, 1}, 42, {}, {}, 0.1);

    EXPECT_GT(result.report.pathologicalCutCount, 0);
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, CrossingSnapToleranceRemovesNearVertexSliver)
  {
    auto exactMesh = makeSingleTriangle();
    computeConnectivity(exactMesh);
    P1 exactSpace(exactMesh);
    GridFunction exactPhi(exactSpace);
    exactPhi[0] = -1e-3;
    exactPhi[1] = 1;
    exactPhi[2] = 1;

    const auto exact = LevelSetDiscretizerTriangles(exactPhi)
      .setSignTolerance(1e-12)
      .setCrossingSnapTolerance(0)
      .setInterfaceAttribute(42)
      .discretize();

    auto snappedMesh = makeSingleTriangle();
    computeConnectivity(snappedMesh);
    P1 snappedSpace(snappedMesh);
    GridFunction snappedPhi(snappedSpace);
    snappedPhi[0] = -1e-3;
    snappedPhi[1] = 1;
    snappedPhi[2] = 1;

    const Real tau = 0.01;
    const auto snapped = LevelSetDiscretizerTriangles(snappedPhi)
      .setSignTolerance(1e-12)
      .setCrossingSnapTolerance(tau)
      .setInterfaceAttribute(42)
      .discretize();

    EXPECT_EQ(exact.report.snappedCrossingCount, 0);
    EXPECT_EQ(exact.mesh.getCellCount(), 3);
    EXPECT_EQ(exact.interfaceGraph.edges.size(), 1);

    EXPECT_EQ(snapped.report.nearVertexCrossingCount, 2);
    EXPECT_EQ(snapped.report.snappedCrossingCount, 2);
    EXPECT_EQ(snapped.interfaceGraph.edges.size(), 0);
    EXPECT_EQ(snapped.mesh.getCellCount(), 1);
    EXPECT_LE(snapped.report.maxInterfaceDeviation, tau);
    EXPECT_GT(minTriangleQuality(snapped.mesh), minTriangleQuality(exact.mesh));
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, QuadCutUsesBetterInteriorDiagonal)
  {
    auto mesh = LocalMesh::Builder()
      .initialize(2)
      .nodes(3)
      .vertex({0, 0})
      .vertex({0.7613897502566094, -0.14872607731865958})
      .vertex({-0.08520978320152678, 0.8421406330301072})
      .polytope(Polytope::Type::Triangle, {0, 1, 2})
      .finalize();
    computeConnectivity(mesh);

    P1 space(mesh);
    GridFunction phi(space);
    phi[0] = -1;
    phi[1] = 7.610724711258564;
    phi[2] = 2.160684226493828;

    const auto result = LevelSetDiscretizerTriangles(phi)
      .setSignTolerance(1e-12)
      .setInterfaceAttribute(42)
      .discretize();

    EXPECT_EQ(result.report.improvedPolygonTriangulationCount, 1);
    EXPECT_EQ(result.report.snappedCrossingCount, 0);
    EXPECT_EQ(result.mesh.getCellCount(), 3);
    EXPECT_GT(minTriangleQuality(result.mesh), 0.44);
    expectEveryGraphEdgeAppearsAsOutputEdge(result);
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, RequiresPreparedConnectivity)
  {
    auto mesh = makeSingleTriangle();

    P1 space(mesh);
    GridFunction phi(space);
    phi[0] = -1;
    phi[1] = 1;
    phi[2] = 1;

    try
    {
      (void) LevelSetDiscretizerTriangles(phi).discretize();
      FAIL() << "Expected missing connectivity to throw.";
    }
    catch (const Alert::Exception& e)
    {
      EXPECT_NE(
        std::string(e.what()).find("has not been computed"),
        std::string::npos);
    }
  }

  struct TriangleTopologyCase
  {
    std::array<Real, 3> phi;
    Index cells;
    Index negative;
    Index positive;
    Index degenerate;
    Index graphEdges;
  };

  class Rodin_Geometry_LevelSetDiscretizerTriangles_StrictTopology
    : public testing::TestWithParam<TriangleTopologyCase>
  {};

  TEST_P(
      Rodin_Geometry_LevelSetDiscretizerTriangles_StrictTopology,
      StrictSignPatternIsCutCorrectly)
  {
    const auto c = GetParam();
    auto mesh = makeSingleTriangle();
    const auto result = discretize(mesh, {c.phi[0], c.phi[1], c.phi[2]});

    EXPECT_EQ(result.mesh.getCellCount(), c.cells);
    EXPECT_EQ(sideCount(result.report, LevelSetSide::Negative), c.negative);
    EXPECT_EQ(sideCount(result.report, LevelSetSide::Positive), c.positive);
    EXPECT_EQ(sideCount(result.report, LevelSetSide::Degenerate), c.degenerate);
    EXPECT_EQ(result.interfaceGraph.edges.size(), c.graphEdges);
    expectNoOutputCellCrossesInterface(
      result, {c.phi[0], c.phi[1], c.phi[2]});
    expectEveryGraphEdgeAppearsAsOutputEdge(result);
  }

  INSTANTIATE_TEST_SUITE_P(
      StrictPatterns,
      Rodin_Geometry_LevelSetDiscretizerTriangles_StrictTopology,
      testing::Values(
        TriangleTopologyCase{{{-1, -1, -1}}, 1, 1, 0, 0, 0},
        TriangleTopologyCase{{{ 1,  1,  1}}, 1, 0, 1, 0, 0},
        TriangleTopologyCase{{{-1,  1,  1}}, 3, 1, 2, 0, 1},
        TriangleTopologyCase{{{ 1, -1,  1}}, 3, 1, 2, 0, 1},
        TriangleTopologyCase{{{ 1,  1, -1}}, 3, 1, 2, 0, 1},
        TriangleTopologyCase{{{ 1, -1, -1}}, 3, 2, 1, 0, 1},
        TriangleTopologyCase{{{-1,  1, -1}}, 3, 2, 1, 0, 1},
        TriangleTopologyCase{{{-1, -1,  1}}, 3, 2, 1, 0, 1}));

  class Rodin_Geometry_LevelSetDiscretizerTriangles_ZeroTopology
    : public testing::TestWithParam<TriangleTopologyCase>
  {};

  TEST_P(
      Rodin_Geometry_LevelSetDiscretizerTriangles_ZeroTopology,
      ZeroSignPatternIsCutCorrectly)
  {
    const auto c = GetParam();
    auto mesh = makeSingleTriangle();
    const auto result = discretize(mesh, {c.phi[0], c.phi[1], c.phi[2]});

    EXPECT_EQ(result.mesh.getCellCount(), c.cells);
    EXPECT_EQ(sideCount(result.report, LevelSetSide::Negative), c.negative);
    EXPECT_EQ(sideCount(result.report, LevelSetSide::Positive), c.positive);
    EXPECT_EQ(sideCount(result.report, LevelSetSide::Degenerate), c.degenerate);
    EXPECT_EQ(result.interfaceGraph.edges.size(), c.graphEdges);
    expectNoOutputCellCrossesInterface(
      result, {c.phi[0], c.phi[1], c.phi[2]});
    expectEveryGraphEdgeAppearsAsOutputEdge(result);
  }

  INSTANTIATE_TEST_SUITE_P(
      ZeroPatterns,
      Rodin_Geometry_LevelSetDiscretizerTriangles_ZeroTopology,
      testing::Values(
        TriangleTopologyCase{{{ 0,  1,  1}}, 1, 0, 1, 0, 0},
        TriangleTopologyCase{{{ 1,  0,  1}}, 1, 0, 1, 0, 0},
        TriangleTopologyCase{{{ 1,  1,  0}}, 1, 0, 1, 0, 0},
        TriangleTopologyCase{{{ 0, -1, -1}}, 1, 1, 0, 0, 0},
        TriangleTopologyCase{{{-1,  0, -1}}, 1, 1, 0, 0, 0},
        TriangleTopologyCase{{{-1, -1,  0}}, 1, 1, 0, 0, 0},
        TriangleTopologyCase{{{ 0, -1,  1}}, 2, 1, 1, 0, 1},
        TriangleTopologyCase{{{ 0,  1, -1}}, 2, 1, 1, 0, 1},
        TriangleTopologyCase{{{-1,  0,  1}}, 2, 1, 1, 0, 1},
        TriangleTopologyCase{{{ 1,  0, -1}}, 2, 1, 1, 0, 1},
        TriangleTopologyCase{{{-1,  1,  0}}, 2, 1, 1, 0, 1},
        TriangleTopologyCase{{{ 1, -1,  0}}, 2, 1, 1, 0, 1},
        TriangleTopologyCase{{{ 0,  0,  1}}, 1, 0, 1, 0, 1},
        TriangleTopologyCase{{{ 0,  1,  0}}, 1, 0, 1, 0, 1},
        TriangleTopologyCase{{{ 1,  0,  0}}, 1, 0, 1, 0, 1},
        TriangleTopologyCase{{{ 0,  0, -1}}, 1, 1, 0, 0, 1},
        TriangleTopologyCase{{{ 0, -1,  0}}, 1, 1, 0, 0, 1},
        TriangleTopologyCase{{{-1,  0,  0}}, 1, 1, 0, 0, 1},
        TriangleTopologyCase{{{ 0,  0,  0}}, 1, 0, 0, 1, 0}));

  struct EdgeCutCase
  {
    std::array<Real, 3> phi;
    Index a;
    Index b;
  };

  class Rodin_Geometry_LevelSetDiscretizerTriangles_EdgeCuts
    : public testing::TestWithParam<EdgeCutCase>
  {};

  TEST_P(
      Rodin_Geometry_LevelSetDiscretizerTriangles_EdgeCuts,
      EdgeCutInterpolationMatchesParentEdge)
  {
    const auto c = GetParam();
    auto mesh = makeSingleTriangle();
    const auto result = discretize(mesh, {c.phi[0], c.phi[1], c.phi[2]});
    const Index parentEdge = findEdge(mesh, c.a, c.b);

    const auto it = std::find_if(
      result.interfaceGraph.vertices.begin(),
      result.interfaceGraph.vertices.end(),
      [parentEdge](const InterfaceVertex& vertex)
      {
        return vertex.kind == InterfaceVertexKind::EdgeCut &&
               vertex.parentEdge &&
               *vertex.parentEdge == parentEdge;
      });
    ASSERT_NE(it, result.interfaceGraph.vertices.end());
    ASSERT_TRUE(it->t);

    const auto& edge = mesh.getConnectivity().getPolytope(1, parentEdge);
    const Real expectedT =
      c.phi[edge(0)] / (c.phi[edge(0)] - c.phi[edge(1)]);
    const auto expectedX =
      (Real(1) - expectedT) * mesh.getVertexCoordinates(edge(0)) +
      expectedT * mesh.getVertexCoordinates(edge(1));

    EXPECT_NEAR(*it->t, expectedT, 1e-14);
    EXPECT_NEAR(it->x[0], expectedX[0], 1e-14);
    EXPECT_NEAR(it->x[1], expectedX[1], 1e-14);
  }

  INSTANTIATE_TEST_SUITE_P(
      ExactInterpolation,
      Rodin_Geometry_LevelSetDiscretizerTriangles_EdgeCuts,
      testing::Values(
        EdgeCutCase{{{-1,  3,  3}}, 0, 1},
        EdgeCutCase{{{-3,  1,  1}}, 0, 1},
        EdgeCutCase{{{ 3, -1,  3}}, 0, 1},
        EdgeCutCase{{{ 1, -3,  1}}, 0, 1},
        EdgeCutCase{{{ 3,  3, -1}}, 2, 0},
        EdgeCutCase{{{ 1,  1, -3}}, 2, 0},
        EdgeCutCase{{{-2,  5,  5}}, 0, 1},
        EdgeCutCase{{{-5,  2,  2}}, 0, 1}));

  struct SquareCutCase
  {
    std::array<Real, 4> phi;
    Index graphEdges;
    bool openChain;
  };

  class Rodin_Geometry_LevelSetDiscretizerTriangles_SquareCuts
    : public testing::TestWithParam<SquareCutCase>
  {};

  TEST_P(
      Rodin_Geometry_LevelSetDiscretizerTriangles_SquareCuts,
      SquareCutsRemainConnectedAndConstrained)
  {
    const auto c = GetParam();
    auto mesh = makeTwoTriangleSquare();
    const auto result =
      discretize(mesh, {c.phi[0], c.phi[1], c.phi[2], c.phi[3]});

    ASSERT_EQ(result.interfaceGraph.edges.size(), c.graphEdges);
    if (c.graphEdges > 0)
      expectEveryGraphEdgeAppearsAsOutputEdge(result);
    if (!result.interfaceGraph.chains.empty())
      EXPECT_EQ(result.interfaceGraph.chains.front().closed, !c.openChain);
    expectNoOutputCellCrossesInterface(
      result, {c.phi[0], c.phi[1], c.phi[2], c.phi[3]});
  }

  INSTANTIATE_TEST_SUITE_P(
      TwoTriangleSquare,
      Rodin_Geometry_LevelSetDiscretizerTriangles_SquareCuts,
      testing::Values(
        SquareCutCase{{{-1,  1, -1,  1}}, 2, true},
        SquareCutCase{{{-1, -1,  1,  1}}, 2, true},
        SquareCutCase{{{ 1, -1,  1, -1}}, 2, true},
        SquareCutCase{{{ 1,  1, -1, -1}}, 2, true},
        SquareCutCase{{{ 1,  0,  0, -1}}, 1, true},
        SquareCutCase{{{-1,  0,  0,  1}}, 1, true},
        SquareCutCase{{{ 0,  1, -1,  0}}, 2, true},
        SquareCutCase{{{ 0, -1,  1,  0}}, 2, true}));

  class Rodin_Geometry_LevelSetDiscretizerTriangles_InvalidPhi
    : public testing::TestWithParam<std::array<Real, 3>>
  {};

  TEST_P(
      Rodin_Geometry_LevelSetDiscretizerTriangles_InvalidPhi,
      InvalidPhiPreservesDegenerateParentCell)
  {
    const auto phi = GetParam();
    auto mesh = makeSingleTriangle();
    const auto result = discretize(mesh, {phi[0], phi[1], phi[2]});

    EXPECT_EQ(result.mesh.getCellCount(), 1);
    EXPECT_EQ(result.report.degenerateCellCount, 1);
    EXPECT_EQ(sideCount(result.report, LevelSetSide::Degenerate), 1);
    EXPECT_TRUE(result.interfaceGraph.invalidCells.size() == 1);
  }

  INSTANTIATE_TEST_SUITE_P(
      InvalidValues,
      Rodin_Geometry_LevelSetDiscretizerTriangles_InvalidPhi,
      testing::Values(
        std::array<Real, 3>{
          std::numeric_limits<Real>::quiet_NaN(), -1, 1},
        std::array<Real, 3>{
          -1, std::numeric_limits<Real>::quiet_NaN(), 1},
        std::array<Real, 3>{
          -1, 1, std::numeric_limits<Real>::quiet_NaN()}));

  class Rodin_Geometry_LevelSetDiscretizerTriangles_Pathology
    : public testing::TestWithParam<std::array<Real, 3>>
  {};

  TEST_P(
      Rodin_Geometry_LevelSetDiscretizerTriangles_Pathology,
      NearVertexCutIsReportedButAccepted)
  {
    const auto phi = GetParam();
    auto mesh = makeSingleTriangle();
    const auto result = discretize(mesh, {phi[0], phi[1], phi[2]}, 42, {}, {}, 0.05);

    EXPECT_GT(result.report.pathologicalCutCount, 0);
    EXPECT_GT(result.mesh.getCellCount(), 0);
    expectEveryGraphEdgeAppearsAsOutputEdge(result);
  }

  INSTANTIATE_TEST_SUITE_P(
      NearVertexCuts,
      Rodin_Geometry_LevelSetDiscretizerTriangles_Pathology,
      testing::Values(
        std::array<Real, 3>{{-0.01, 1, 1}},
        std::array<Real, 3>{{1, -0.01, 1}},
        std::array<Real, 3>{{1, 1, -0.01}},
        std::array<Real, 3>{{0.01, -1, -1}},
        std::array<Real, 3>{{-1, 0.01, -1}},
        std::array<Real, 3>{{-1, -1, 0.01}}));

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, MinCutQualityZeroPreservesExactCut)
  {
    auto mesh = makeSingleTriangle();
    computeConnectivity(mesh);
    P1 space(mesh);
    GridFunction phi(space);
    phi[0] = -1e-3; phi[1] = 1; phi[2] = 1;  // near-vertex sliver cut

    const auto r = LevelSetDiscretizerTriangles(phi)
      .setSignTolerance(1e-12)
      .setInterfaceAttribute(42)
      .setMinCutQuality(0)
      .discretize();

    EXPECT_EQ(r.report.uncutCellCount, 0);
    EXPECT_EQ(r.mesh.getCellCount(), 3);
    EXPECT_EQ(r.interfaceGraph.edges.size(), 1);
  }

  static Real totalArea(const LocalMesh& m)
  {
    const auto& conn = m.getConnectivity();
    Real s = 0;
    for (Index c = 0; c < static_cast<Index>(m.getCellCount()); ++c)
    {
      const auto& t = conn.getPolytope(2, c);
      s += std::abs(signedTriangleArea(m.getVertexCoordinates(t(0)),
                                       m.getVertexCoordinates(t(1)),
                                       m.getVertexCoordinates(t(2))));
    }
    return s;
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, MinCutQualityNeverDropsAndStaysConforming)
  {
    // A near-vertex crossing that would split into a thin sliver.
    auto mesh = makeSingleTriangle();
    computeConnectivity(mesh);
    P1 space(mesh);
    GridFunction phi(space);
    phi[0] = -1e-3; phi[1] = 1; phi[2] = 1;

    // With no snap, the crossing is a real interior vertex: keeping the cell
    // whole would create a hanging node, so the cut MUST proceed (conformity
    // wins). Crucially, the produced sliver is never dropped: the output is
    // watertight (area conserved).
    const auto cut = LevelSetDiscretizerTriangles(phi)
      .setSignTolerance(1e-12)
      .setInterfaceAttribute(42)
      .setMinCutQuality(0.1)
      .discretize();
    EXPECT_EQ(cut.report.uncutCellCount, 0u);
    EXPECT_EQ(cut.mesh.getCellCount(), 3);
    EXPECT_NEAR(totalArea(cut.mesh), 0.5, 1e-12);   // nothing dropped/holed

    // The conforming way to avoid the bad cut is the per-edge snap, which is
    // shared by both incident cells (no hanging node). It removes the sliver.
    auto mesh2 = makeSingleTriangle();
    computeConnectivity(mesh2);
    P1 s2(mesh2);
    GridFunction phi2(s2);
    phi2[0] = -1e-3; phi2[1] = 1; phi2[2] = 1;
    const auto snapped = LevelSetDiscretizerTriangles(phi2)
      .setSignTolerance(1e-12)
      .setInterfaceAttribute(42)
      .setCrossingSnapTolerance(0.01)
      .discretize();
    EXPECT_EQ(snapped.mesh.getCellCount(), 1);
    EXPECT_NEAR(totalArea(snapped.mesh), 0.5, 1e-12);
  }

  TEST(Rodin_Geometry_LevelSetDiscretizerTriangles, CustomQualityMeasureIsUsed)
  {
    auto mesh = makeSingleTriangle();
    computeConnectivity(mesh);
    P1 space(mesh);
    GridFunction phi(space);
    phi[0] = -1; phi[1] = 1; phi[2] = 1;

    bool called = false;
    const auto r = LevelSetDiscretizerTriangles(phi)
      .setSignTolerance(1e-12)
      .setInterfaceAttribute(42)
      .setQualityMeasure(
        [&](const Math::SpatialPoint&, const Math::SpatialPoint&,
            const Math::SpatialPoint&) { called = true; return Real(0.5); })
      .discretize();

    EXPECT_TRUE(called);
    EXPECT_NEAR(r.report.minOutputCellQuality, 0.5, 1e-12);
  }
}
