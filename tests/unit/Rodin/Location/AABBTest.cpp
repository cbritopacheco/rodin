/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <gtest/gtest.h>

#include <chrono>
#include <cmath>
#include <limits>
#include <random>
#include <thread>
#include <vector>

#include "Rodin/Location.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Location;

namespace
{
  Math::SpatialPoint point(std::initializer_list<Real> coords)
  {
    Math::SpatialPoint res(coords.size());
    Eigen::Index i = 0;
    for (Real x : coords)
      res[i++] = x;
    return res;
  }

  Math::SpatialPoint physicalCentroid(
    const Mesh<Context::Local>& mesh, const Polytope& polytope)
  {
    Math::SpatialPoint x;
    Polytope::Traits traits(polytope.getGeometry());
    mesh.getPolytopeTransformation(polytope.getDimension(), polytope.getIndex())
      .transform(x, traits.getCentroid());
    return x;
  }
}

TEST(Location_AABB, Locates0DVertex)
{
  Mesh mesh = Mesh<Context::Local>::Builder()
                .initialize(2)
                .nodes(1)
                .vertex({0.25, 0.75})
                .finalize();

  AABB locator(mesh);
  const auto p = locator.locate(0, point({0.25, 0.75}));
  ASSERT_TRUE(p.has_value());
  EXPECT_EQ(p->getPolytope().getDimension(), 0);
  EXPECT_EQ(p->getPolytope().getIndex(), 0);
}

TEST(Location_AABB, Locates1DSegment)
{
  Mesh mesh = Mesh<Context::Local>::Builder()
                .initialize(1)
                .nodes(2)
                .vertex({0.0})
                .vertex({1.0})
                .polytope(Polytope::Type::Segment, {0, 1})
                .finalize();

  AABB locator(mesh);
  const auto p = locator.locate(point({0.4}));
  ASSERT_TRUE(p.has_value());
  EXPECT_EQ(p->getPolytope().getDimension(), 1);
  EXPECT_NEAR(p->getReferenceCoordinates()[0], 0.4, 1e-12);
}

TEST(Location_AABB, Locates2DTriangles)
{
  Mesh mesh = Mesh<Context::Local>::Builder()
                .initialize(2)
                .nodes(4)
                .vertex({0.0, 0.0})
                .vertex({1.0, 0.0})
                .vertex({0.0, 1.0})
                .vertex({1.0, 1.0})
                .polytope(Polytope::Type::Triangle, {0, 1, 2})
                .polytope(Polytope::Type::Triangle, {1, 3, 2})
                .finalize();

  AABB locator(mesh);
  const auto p0 = locator.locate(point({0.2, 0.2}));
  ASSERT_TRUE(p0.has_value());
  EXPECT_EQ(p0->getPolytope().getIndex(), 0);

  const auto p1 = locator.locate(point({0.8, 0.8}));
  ASSERT_TRUE(p1.has_value());
  EXPECT_EQ(p1->getPolytope().getIndex(), 1);
}

TEST(Location_AABB, Locates3DTetrahedron)
{
  Mesh mesh = Mesh<Context::Local>::Builder()
                .initialize(3)
                .nodes(4)
                .vertex({0.0, 0.0, 0.0})
                .vertex({1.0, 0.0, 0.0})
                .vertex({0.0, 1.0, 0.0})
                .vertex({0.0, 0.0, 1.0})
                .polytope(Polytope::Type::Tetrahedron, {0, 1, 2, 3})
                .finalize();

  AABB locator(mesh);
  const auto p = locator.locate(point({0.2, 0.2, 0.2}));
  ASSERT_TRUE(p.has_value());
  EXPECT_EQ(p->getPolytope().getDimension(), 3);
  EXPECT_NEAR(p->getReferenceCoordinates()[0], 0.2, 1e-12);
  EXPECT_NEAR(p->getReferenceCoordinates()[1], 0.2, 1e-12);
  EXPECT_NEAR(p->getReferenceCoordinates()[2], 0.2, 1e-12);
}

TEST(Location_AABB, LocatesMixed2DCells)
{
  Mesh mesh = Mesh<Context::Local>::Builder()
                .initialize(2)
                .nodes(6)
                .vertex({0.0, 0.0})
                .vertex({1.0, 0.0})
                .vertex({0.0, 1.0})
                .vertex({1.0, 1.0})
                .vertex({2.0, 0.0})
                .vertex({2.0, 1.0})
                .polytope(Polytope::Type::Triangle, {0, 1, 2})
                .polytope(Polytope::Type::Quadrilateral, {1, 4, 5, 3})
                .finalize();

  AABB locator(mesh);
  const auto tri = locator.locate(point({0.2, 0.2}));
  ASSERT_TRUE(tri.has_value());
  EXPECT_EQ(tri->getPolytope().getGeometry(), Polytope::Type::Triangle);

  const auto quad = locator.locate(point({1.5, 0.5}));
  ASSERT_TRUE(quad.has_value());
  EXPECT_EQ(quad->getPolytope().getGeometry(), Polytope::Type::Quadrilateral);
}

TEST(Location_AABB, LocatesModerateUniformGridCentroids)
{
  const Mesh mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {12, 12});
  AABB locator(mesh);

  size_t count = 0;
  for (auto it = mesh.getCell(); it; ++it)
  {
    const Math::SpatialPoint x = physicalCentroid(mesh, *it);
    const auto p = locator.locate(x);
    ASSERT_TRUE(p.has_value());
    EXPECT_EQ(p->getPolytope().getDimension(), mesh.getDimension());
    count++;
  }
  EXPECT_EQ(count, mesh.getCellCount());
}

TEST(Location_AABB, MissesOutsidePoint)
{
  const Mesh mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {3, 3});
  AABB locator(mesh);

  const auto p = locator.locate(point({4.0, 4.0}));
  EXPECT_FALSE(p.has_value());
}

TEST(Location_AABB, RejectsNonFiniteQuery)
{
  const Mesh mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {3, 3});
  AABB locator(mesh);

  const Real nan = std::numeric_limits<Real>::quiet_NaN();
  const Real inf = std::numeric_limits<Real>::infinity();
  EXPECT_FALSE(locator.locate(point({nan, 0.5})).has_value());
  EXPECT_FALSE(locator.locate(point({0.5, nan})).has_value());
  EXPECT_FALSE(locator.locate(point({inf, 0.5})).has_value());
  EXPECT_FALSE(locator.locate(point({-inf, -inf})).has_value());
}

TEST(Location_AABB, RejectsWrongDimensionQuery)
{
  const Mesh mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {3, 3});
  AABB locator(mesh);

  EXPECT_FALSE(locator.locate(point({0.5})).has_value());
  EXPECT_FALSE(locator.locate(point({0.5, 0.5, 0.5})).has_value());
  EXPECT_FALSE(locator.locate(5, point({0.5, 0.5})).has_value());
}

TEST(Location_AABB, RejectsOffManifoldSurfaceQuery)
{
  // A single segment embedded in 2D: locating against the segment must
  // reject points off the line even when their projection parameter is
  // inside [0, 1].
  Mesh mesh = Mesh<Context::Local>::Builder()
                .initialize(2)
                .nodes(2)
                .vertex({0.0, 0.0})
                .vertex({1.0, 0.0})
                .polytope(Polytope::Type::Segment, {0, 1})
                .finalize();

  AABB locator(mesh);
  EXPECT_FALSE(locator.locate(1, point({0.5, 0.5})).has_value());
  EXPECT_FALSE(locator.locate(1, point({0.5, 1e-3})).has_value());

  const auto on = locator.locate(1, point({0.5, 0.0}));
  ASSERT_TRUE(on.has_value());
  EXPECT_NEAR(on->getReferenceCoordinates()[0], 0.5, 1e-12);
}

TEST(Location_AABB, LocatesInsideNonParallelogramQuad)
{
  // The bilinear map of a trapezoid is not affine; the narrow phase must
  // invert it exactly (Newton), not by linearization.
  Mesh mesh = Mesh<Context::Local>::Builder()
                .initialize(2)
                .nodes(4)
                .vertex({0.0, 0.0})
                .vertex({2.0, 0.0})
                .vertex({1.5, 1.0})
                .vertex({0.5, 1.0})
                .polytope(Polytope::Type::Quadrilateral, {0, 1, 3, 2})
                .finalize();

  AABB locator(mesh);
  const auto p = locator.locate(point({1.0, 0.9}));
  ASSERT_TRUE(p.has_value());

  // Verify the inverse by mapping the reference coordinates forward.
  Math::SpatialPoint mapped;
  p->getPolytope().getTransformation().transform(mapped, p->getReferenceCoordinates());
  EXPECT_NEAR(mapped[0], 1.0, 1e-10);
  EXPECT_NEAR(mapped[1], 0.9, 1e-10);

  // Outside the slanted edge but inside the vertex bounding box.
  EXPECT_FALSE(locator.locate(point({0.1, 0.95})).has_value());
}

TEST(Location_AABB, LocatesInsideCurvedP2Bulge)
{
  // A single triangle whose bottom edge is bent downward by a quadratic
  // (P2) transformation. Points inside the bulge lie outside the straight
  // triangle: the broad phase must be conservative for the curved box and
  // the narrow phase must converge on the curved geometry.
  Mesh mesh = Mesh<Context::Local>::Builder()
                .initialize(2)
                .nodes(3)
                .vertex({0.0, 0.0})
                .vertex({1.0, 0.0})
                .vertex({0.0, 1.0})
                .polytope(Polytope::Type::Triangle, {0, 1, 2})
                .finalize();

  const auto cell = mesh.getCell(0);
  Variational::RealH1Element<2> geomFe(Polytope::Type::Triangle);
  Geometry::PointCloud pm(2, geomFe.getCount());
  Math::SpatialPoint x;
  for (size_t a = 0; a < geomFe.getCount(); ++a)
  {
    cell->getTransformation().transform(x, geomFe.getNode(a));
    // Push the node at the bottom-edge midpoint (0.5, 0) downward.
    if (std::abs(x[0] - 0.5) < 1e-12 && std::abs(x[1]) < 1e-12)
      x[1] -= 0.25;
    pm(0, a) = x[0];
    pm(1, a) = x[1];
  }
  mesh.setPolytopeTransformation({2, 0},
    new Geometry::ParametricTransformation<Variational::RealH1Element<2>>(
      std::move(pm), geomFe));

  AABB locator(mesh);

  // Inside the bulge, below the straight edge y = 0.
  const auto p = locator.locate(point({0.5, -0.1}));
  ASSERT_TRUE(p.has_value());
  Math::SpatialPoint mapped;
  p->getPolytope().getTransformation().transform(mapped, p->getReferenceCoordinates());
  EXPECT_NEAR(mapped[0], 0.5, 1e-9);
  EXPECT_NEAR(mapped[1], -0.1, 1e-9);

  // Below the bulge apex: outside the curved cell.
  EXPECT_FALSE(locator.locate(point({0.5, -0.3})).has_value());
}

TEST(Location_AABB, LocatesAtExtremeMeshScales)
{
  for (const Real scale : {Real(1e-6), Real(1e+6)})
  {
    Mesh mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {8, 8});
    mesh.scale(scale);

    AABB locator(mesh);
    size_t found = 0;
    for (auto it = mesh.getCell(); it; ++it)
      found += locator.locate(physicalCentroid(mesh, *it)).has_value();
    EXPECT_EQ(found, mesh.getCellCount()) << "scale=" << scale;

    EXPECT_FALSE(locator.locate(point({Real(20) * scale, Real(-3) * scale})).has_value())
      << "scale=" << scale;
  }
}

TEST(Location_AABB, SharedFacetPointIsDeterministicAndValid)
{
  const Mesh mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {4, 4});
  AABB locator(mesh);

  // A point on the shared diagonal of two triangles.
  const auto x = point({1.5, 1.5});
  const auto first = locator.locate(x);
  ASSERT_TRUE(first.has_value());

  Math::SpatialPoint mapped;
  first->getPolytope().getTransformation().transform(
    mapped, first->getReferenceCoordinates());
  EXPECT_NEAR(mapped[0], x[0], 1e-10);
  EXPECT_NEAR(mapped[1], x[1], 1e-10);

  for (int repeat = 0; repeat < 8; ++repeat)
  {
    const auto again = locator.locate(x);
    ASSERT_TRUE(again.has_value());
    EXPECT_EQ(again->getPolytope().getIndex(), first->getPolytope().getIndex());
  }
}

TEST(Location_AABB, MatchesBruteForceOnRandomPoints)
{
  const Mesh mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {12, 12});
  AABB locator(mesh);

  // Independent ground truth for straight triangles: barycentric containment.
  const auto inTriangle = [&](Index cell, const Math::SpatialPoint& x, Real tol) -> bool {
    const auto& v = mesh.getCell(cell)->getVertices();
    const auto a = mesh.getVertexCoordinates(v[0]);
    const auto b = mesh.getVertexCoordinates(v[1]);
    const auto c = mesh.getVertexCoordinates(v[2]);
    const Real det = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]);
    const Real l1 = ((x[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (x[1] - a[1])) / det;
    const Real l2 = ((b[0] - a[0]) * (x[1] - a[1]) - (x[0] - a[0]) * (b[1] - a[1])) / det;
    return l1 >= -tol && l2 >= -tol && l1 + l2 <= Real(1) + tol;
  };

  std::mt19937 rng(20260703);
  std::uniform_real_distribution<Real> dist(-2.0, 14.0);
  size_t hits = 0, misses = 0;
  for (size_t trial = 0; trial < 4000; ++trial)
  {
    const auto x = point({dist(rng), dist(rng)});
    const auto p = locator.locate(x);
    bool any = false;
    for (Index cell = 0; cell < mesh.getCellCount() && !any; ++cell)
      any = inTriangle(cell, x, Real(0));
    if (p.has_value())
    {
      ++hits;
      // The located cell must actually contain the point (loose tolerance
      // to absorb the locator's boundary slack).
      EXPECT_TRUE(inTriangle(p->getPolytope().getIndex(), x, Real(1e-8)))
        << "point (" << x[0] << ", " << x[1] << ")";
    }
    else
    {
      ++misses;
      // A missed point must not be strictly inside any cell.
      EXPECT_FALSE(any) << "point (" << x[0] << ", " << x[1] << ")";
    }
  }
  // Sanity: the sample must exercise both branches.
  EXPECT_GT(hits, size_t(500));
  EXPECT_GT(misses, size_t(500));
}

TEST(Location_AABB, ConcurrentQueriesAreSafe)
{
  const Mesh mesh = Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {24, 24});
  AABB locator(mesh);

  std::vector<Math::SpatialPoint> queries;
  for (auto it = mesh.getCell(); it; ++it)
    queries.push_back(physicalCentroid(mesh, *it));

  // First queries race on the lazy build; all must succeed.
  std::atomic<size_t> found{0};
  std::vector<std::thread> threads;
  for (int t = 0; t < 4; ++t)
  {
    threads.emplace_back([&] {
      size_t local = 0;
      for (const auto& x : queries)
        local += locator.locate(x).has_value();
      found += local;
    });
  }
  for (auto& thread : threads)
    thread.join();
  EXPECT_EQ(found.load(), 4 * queries.size());
}

TEST(Location_AABB, PerfSmoke)
{
  using Clock = std::chrono::steady_clock;
  const auto usPerQuery = [](Clock::time_point t0, Clock::time_point t1, size_t n) {
    return std::chrono::duration<double, std::micro>(t1 - t0).count() /
      static_cast<double>(n);
  };

  const Mesh mesh =
    Mesh<Context::Local>::UniformGrid(Polytope::Type::Triangle, {100, 100});

  const auto tBuild0 = Clock::now();
  AABB locator(mesh);
  // The index is built lazily on first query; charge it to the build phase.
  locator.locate(point({0.5, 0.5}));
  const auto tBuild1 = Clock::now();

  std::vector<Math::SpatialPoint> hits;
  for (auto it = mesh.getCell(); it; ++it)
    hits.push_back(physicalCentroid(mesh, *it));

  const auto tHit0 = Clock::now();
  size_t found = 0;
  for (const auto& x : hits)
    found += locator.locate(x).has_value();
  const auto tHit1 = Clock::now();
  EXPECT_EQ(found, hits.size());

  const size_t missCount = 200;
  const auto tMiss0 = Clock::now();
  size_t missFound = 0;
  for (size_t i = 0; i < missCount; ++i)
    missFound += locator.locate(point({200.0 + Real(i), -50.0})).has_value();
  const auto tMiss1 = Clock::now();
  EXPECT_EQ(missFound, size_t(0));

  // Near-misses: just outside the boundary but inside the root box, so the
  // tree is traversed and the narrow phase must reject.
  const size_t nearCount = 2000;
  const auto tNear0 = Clock::now();
  size_t nearFound = 0;
  for (size_t i = 0; i < nearCount; ++i)
  {
    // Inside the box padding band, so the narrow phase must do the rejection.
    const Real s = Real(100) * Real(i) / Real(nearCount);
    nearFound += locator.locate(point({s, -1e-9})).has_value();
  }
  const auto tNear1 = Clock::now();
  EXPECT_EQ(nearFound, size_t(0));

  std::cout << "[ PERF 2D  ] build="
            << std::chrono::duration<double, std::milli>(tBuild1 - tBuild0).count()
            << "ms  hit=" << usPerQuery(tHit0, tHit1, hits.size())
            << "us/q  farMiss=" << usPerQuery(tMiss0, tMiss1, missCount)
            << "us/q  nearMiss=" << usPerQuery(tNear0, tNear1, nearCount) << "us/q  ("
            << mesh.getCellCount() << " cells)\n";
}

TEST(Location_AABB, PerfSmoke3D)
{
  using Clock = std::chrono::steady_clock;
  const auto usPerQuery = [](Clock::time_point t0, Clock::time_point t1, size_t n) {
    return std::chrono::duration<double, std::micro>(t1 - t0).count() /
      static_cast<double>(n);
  };

  const Mesh mesh =
    Mesh<Context::Local>::UniformGrid(Polytope::Type::Tetrahedron, {20, 20, 20});

  const auto tBuild0 = Clock::now();
  AABB locator(mesh);
  locator.locate(point({0.5, 0.5, 0.5}));
  const auto tBuild1 = Clock::now();

  std::vector<Math::SpatialPoint> hits;
  for (auto it = mesh.getCell(); it; ++it)
    hits.push_back(physicalCentroid(mesh, *it));

  const auto tHit0 = Clock::now();
  size_t found = 0;
  for (const auto& x : hits)
    found += locator.locate(x).has_value();
  const auto tHit1 = Clock::now();
  EXPECT_EQ(found, hits.size());

  const size_t nearCount = 2000;
  const auto tNear0 = Clock::now();
  size_t nearFound = 0;
  for (size_t i = 0; i < nearCount; ++i)
  {
    // Inside the box padding band, so the narrow phase must do the rejection.
    const Real s = Real(19) * Real(i) / Real(nearCount);
    nearFound += locator.locate(point({s, s / Real(2), -1e-9})).has_value();
  }
  const auto tNear1 = Clock::now();
  EXPECT_EQ(nearFound, size_t(0));

  std::cout << "[ PERF 3D  ] build="
            << std::chrono::duration<double, std::milli>(tBuild1 - tBuild0).count()
            << "ms  hit=" << usPerQuery(tHit0, tHit1, hits.size())
            << "us/q  nearMiss=" << usPerQuery(tNear0, tNear1, nearCount) << "us/q  ("
            << mesh.getCellCount() << " cells)\n";
}
