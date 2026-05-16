/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <cmath>
#include <limits>

#include <gtest/gtest.h>

#include <Rodin/Geometry.h>
#include <Rodin/Geometry/TriangleMeshOptimizer.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Unit
{
  static LocalMesh grid(size_t n, Real scale = 1)
  {
    auto m = LocalMesh::UniformGrid(Polytope::Type::Triangle, { n, n });
    if (scale != Real(1))
      m.scale(scale);
    m.getConnectivity().compute(2, 1);
    m.getConnectivity().compute(1, 0);
    return m;
  }

  static Real minSignedArea(const LocalMesh& m)
  {
    const auto& conn = m.getConnectivity();
    Real q = std::numeric_limits<Real>::infinity();
    for (Index c = 0; c < static_cast<Index>(m.getCellCount()); ++c)
    {
      const auto& t = conn.getPolytope(2, c);
      const auto a = m.getVertexCoordinates(t(0));
      const auto b = m.getVertexCoordinates(t(1));
      const auto cc = m.getVertexCoordinates(t(2));
      const Real ar = Real(0.5) * ((b[0]-a[0])*(cc[1]-a[1])
                                 - (b[1]-a[1])*(cc[0]-a[0]));
      q = std::min(q, ar);
    }
    return q;
  }

  static std::array<Real, 4> bbox(const LocalMesh& m)
  {
    Real x0 = 1e30, y0 = 1e30, x1 = -1e30, y1 = -1e30;
    for (Index v = 0; v < static_cast<Index>(m.getVertexCount()); ++v)
    {
      const auto p = m.getVertexCoordinates(v);
      x0 = std::min(x0, p[0]); x1 = std::max(x1, p[0]);
      y0 = std::min(y0, p[1]); y1 = std::max(y1, p[1]);
    }
    return { x0, y0, x1, y1 };
  }

  TEST(Rodin_Geometry_TriangleMeshOptimizer, IdentityWithinBounds)
  {
    auto m = grid(5);
    const auto cells = m.getCellCount();
    const auto r = TriangleMeshOptimizer()
      .setHMin(0).setHMax(std::numeric_limits<Real>::infinity())
      .optimize(m);
    EXPECT_EQ(r.splits, 0u);
    EXPECT_EQ(r.collapses, 0u);
    EXPECT_EQ(m.getCellCount(), cells);
  }

  TEST(Rodin_Geometry_TriangleMeshOptimizer, SplitRefinesLongEdges)
  {
    auto m = grid(4);                       // raw coords 0..3, edges >= 1
    const auto before = m.getCellCount();
    const auto r = TriangleMeshOptimizer()
      .setHMax(0.8).setMaxIterations(3)
      .optimize(m);
    EXPECT_GT(r.splits, 0u);
    EXPECT_GT(m.getCellCount(), before);
    EXPECT_GT(minSignedArea(m), 0.0);       // no inverted/degenerate cells
  }

  TEST(Rodin_Geometry_TriangleMeshOptimizer, CollapseCoarsensAndPreservesBoundary)
  {
    auto m = grid(9, 1.0 / 8.0);            // [0,1]^2, fine elements ~0.125
    const auto before = m.getCellCount();
    const auto bb0 = bbox(m);
    const auto r = TriangleMeshOptimizer()
      .setHMin(0.4).setHMax(std::numeric_limits<Real>::infinity())
      .setMaxIterations(5)
      .optimize(m);
    EXPECT_GT(r.collapses, 0u);
    EXPECT_LT(m.getCellCount(), before);
    EXPECT_GT(minSignedArea(m), 0.0);       // valid, no inversion
    const auto bb1 = bbox(m);
    for (int i = 0; i < 4; ++i)
      EXPECT_NEAR(bb0[i], bb1[i], 1e-12);   // domain boundary preserved
  }

  TEST(Rodin_Geometry_TriangleMeshOptimizer, ProtectedVerticesAreFrozen)
  {
    auto m = grid(7, 1.0 / 6.0);                 // fine: would collapse a lot
    const Index nv = m.getVertexCount();
    std::vector<std::array<Real, 2>> before(nv);
    for (Index v = 0; v < nv; ++v)
    {
      const auto p = m.getVertexCoordinates(v);
      before[v] = { p[0], p[1] };
    }
    // Protect every vertex: no operation may fire.
    std::vector<char> all(nv, 1);
    const auto r = TriangleMeshOptimizer()
      .setHMin(0.4).setHMax(0.05).setMaxIterations(5)
      .setProtectedVertices(all)
      .optimize(m);
    EXPECT_EQ(r.splits, 0u);
    EXPECT_EQ(r.collapses, 0u);
    EXPECT_EQ(r.swaps, 0u);
    ASSERT_EQ(m.getVertexCount(), nv);
    for (Index v = 0; v < nv; ++v)
    {
      const auto p = m.getVertexCoordinates(v);
      EXPECT_NEAR(p[0], before[v][0], 1e-15);
      EXPECT_NEAR(p[1], before[v][1], 1e-15);
    }
  }

  TEST(Rodin_Geometry_TriangleMeshOptimizer, ProtectedSubsetSurvivesCoarsening)
  {
    auto m = grid(9, 1.0 / 8.0);
    const Index nv = m.getVertexCount();
    // Freeze a diagonal band of vertices; collapse everything else.
    std::vector<char> prot(nv, 0);
    std::vector<std::array<Real, 2>> kept;
    for (Index v = 0; v < nv; ++v)
    {
      const auto p = m.getVertexCoordinates(v);
      if (std::abs(p[0] - p[1]) < 1e-9)
      { prot[v] = 1; kept.push_back({ p[0], p[1] }); }
    }
    TriangleMeshOptimizer()
      .setHMin(0.5).setHMax(std::numeric_limits<Real>::infinity())
      .setMaxIterations(5)
      .setProtectedVertices(prot)
      .optimize(m);

    // Every protected coordinate must still be present in the output.
    for (const auto& k : kept)
    {
      bool found = false;
      for (Index v = 0; v < static_cast<Index>(m.getVertexCount()); ++v)
      {
        const auto p = m.getVertexCoordinates(v);
        if (std::abs(p[0]-k[0]) < 1e-12 && std::abs(p[1]-k[1]) < 1e-12)
        { found = true; break; }
      }
      EXPECT_TRUE(found);
    }
    EXPECT_GT(minSignedArea(m), 0.0);
  }

  TEST(Rodin_Geometry_TriangleMeshOptimizer, SwapImprovesWorstQualityWithoutInversion)
  {
    // A skewed quad split along its bad diagonal; a swap should help and
    // never invert.
    auto mesh = LocalMesh::Builder()
      .initialize(2)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({1, 1})
      .vertex({0.05, 0.05})
      .polytope(Polytope::Type::Triangle, {0, 1, 3})
      .polytope(Polytope::Type::Triangle, {1, 2, 3})
      .finalize();
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    const auto r = TriangleMeshOptimizer()
      .setHMin(0).setHMax(std::numeric_limits<Real>::infinity())
      .setMaxIterations(3)
      .optimize(mesh);

    EXPECT_GT(minSignedArea(mesh), 0.0);
    EXPECT_GE(r.minQualityAfter, r.minQualityBefore);
  }
}
