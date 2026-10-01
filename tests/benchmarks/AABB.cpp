/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief AABB build and point-location benchmarks for every polytope type.
 *
 * Mesh construction and query generation are outside the timed loop. Hit
 * queries use non-centroid reference coordinates so inversion is measured.
 * Point has no curved geometry; the other seven types also have P2 hit cases.
 */

#include <benchmark/benchmark.h>

#include <Rodin/Location.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Benchmarks
{
  namespace
  {
    using MeshType = Mesh<Context::Local>;
    using G = Polytope::Type;

    MeshType makeMesh(G type, bool curved = false)
    {
      if (type == G::Point)
      {
        return MeshType::Builder()
          .initialize(3)
          .nodes(1)
          .vertex({0.25, 0.5, 0.75})
          .finalize();
      }
      const size_t dimension = Polytope::Traits(type).getDimension();
      Array<size_t> grid(dimension);
      size_t resolution = curved ? 8 : 8192;
      if (dimension == 2)
        resolution = curved ? 8 : 64;
      if (dimension == 3)
        resolution = curved ? 4 : 12;
      grid.setConstant(resolution);
      return MeshType::UniformGrid(type, grid);
    }

    std::vector<Math::SpatialPoint> mappedQueries(const MeshType& mesh, G type)
    {
      if (type == G::Point)
        return {mesh.getVertexCoordinates(0)};

      const Polytope::Traits traits(type);
      const Math::SpatialPoint reference =
        Real(0.75) * traits.getCentroid() + Real(0.25) * traits.getVertex(0);
      std::vector<Math::SpatialPoint> queries;
      queries.reserve(mesh.getCellCount());
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        Math::SpatialPoint x;
        cell->getTransformation().transform(x, reference);
        queries.push_back(std::move(x));
      }
      return queries;
    }

    template <size_t K = 2>
    void curveMesh(MeshType& mesh, G type)
    {
      const size_t dimension = mesh.getDimension();
      Variational::RealH1Element<K> element(type);
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        PointCloud nodes(dimension, element.getCount());
        for (size_t a = 0; a < element.getCount(); ++a)
        {
          Math::SpatialPoint x;
          cell->getTransformation().transform(x, element.getNode(a));
          for (size_t i = 0; i < dimension; ++i)
            nodes(i, a) = x[i] + (i == dimension - 1 ? 0.15 * x[0] * x[0] : 0);
        }
        mesh.setPolytopeTransformation({dimension, cell->getIndex()},
          new ParametricTransformation<Variational::RealH1Element<K>>(
            std::move(nodes), element));
      }
    }

    bool allLocated(const Location::AABB<MeshType>& locator,
      const std::vector<Math::SpatialPoint>& queries)
    {
      for (const auto& x : queries)
      {
        if (!locator.locate(x))
          return false;
      }
      return true;
    }

    void BM_AABBBuild(benchmark::State& state, G type)
    {
      const MeshType mesh = makeMesh(type);
      const auto queries = mappedQueries(mesh, type);
      for (auto _ : state)
      {
        Location::AABB locator(mesh);
        benchmark::DoNotOptimize(locator.locate(queries.front()).has_value());
      }
      state.counters["cells"] = static_cast<double>(mesh.getCellCount());
    }

    void BM_AABBBuildWithProjections(benchmark::State& state, G type)
    {
      const MeshType mesh = makeMesh(type);
      const auto queries = mappedQueries(mesh, type);
      for (auto _ : state)
      {
        Location::AABB locator(mesh);
        locator.setProjectionPruning(true);
        benchmark::DoNotOptimize(locator.locate(queries.front()).has_value());
      }
    }

    void BM_AABBHit(benchmark::State& state, G type, bool pruning = false)
    {
      const MeshType mesh = makeMesh(type);
      const auto queries = mappedQueries(mesh, type);
      Location::AABB locator(mesh);
      locator.setProjectionPruning(pruning);
      if (!allLocated(locator, queries))
      {
        state.SkipWithError("AABB missed a mapped query");
        return;
      }
      size_t i = 0;
      for (auto _ : state)
      {
        benchmark::DoNotOptimize(locator.locate(queries[i]).has_value());
        i = (i + 1) % queries.size();
      }
      state.counters["cells"] = static_cast<double>(mesh.getCellCount());
    }

    void BM_AABBOuterMiss(benchmark::State& state, G type)
    {
      const MeshType mesh = makeMesh(type);
      Location::AABB locator(mesh);
      Math::SpatialPoint outside(mesh.getSpaceDimension());
      outside.setConstant(-1);
      if (!locator.locate(mappedQueries(mesh, type).front()) || locator.locate(outside))
      {
        state.SkipWithError("AABB setup failed");
        return;
      }
      for (auto _ : state)
        benchmark::DoNotOptimize(locator.locate(outside).has_value());
      state.counters["cells"] = static_cast<double>(mesh.getCellCount());
    }

    void BM_AABBCurvedHit(benchmark::State& state, G type, bool pruning = false)
    {
      MeshType mesh = makeMesh(type, true);
      curveMesh(mesh, type);
      const auto queries = mappedQueries(mesh, type);
      Location::AABB locator(mesh);
      locator.setProjectionPruning(pruning);
      if (!allLocated(locator, queries))
      {
        state.SkipWithError("AABB missed a mapped curved query");
        return;
      }
      size_t i = 0;
      for (auto _ : state)
      {
        benchmark::DoNotOptimize(locator.locate(queries[i]).has_value());
        i = (i + 1) % queries.size();
      }
      state.counters["cells"] = static_cast<double>(mesh.getCellCount());
    }

    template <size_t K>
    void BM_AABBCurvedBuild(benchmark::State& state, G type, bool pruning = false)
    {
      MeshType mesh = makeMesh(type, true);
      curveMesh<K>(mesh, type);
      const auto queries = mappedQueries(mesh, type);
      for (auto _ : state)
      {
        Location::AABB locator(mesh);
        locator.setProjectionPruning(pruning);
        benchmark::DoNotOptimize(locator.locate(queries.front()).has_value());
      }
      // Includes per-cell sampling, conversion and BVH build. The process-wide
      // factor conversion cache is reused after its first construction.
      state.counters["cells"] = static_cast<double>(mesh.getCellCount());
      state.counters["degree"] = static_cast<double>(K);
    }

    void BM_AABBEmbeddedCurvedHit(benchmark::State& state)
    {
      MeshType mesh = MeshType::Builder()
                        .initialize(3)
                        .nodes(2)
                        .vertex({0, 0, 0})
                        .vertex({1, 0, 0})
                        .polytope(G::Segment, {0, 1})
                        .finalize();
      Variational::RealH1Element<2> element(G::Segment);
      PointCloud nodes(3, element.getCount());
      for (size_t a = 0; a < element.getCount(); ++a)
      {
        const Real t = element.getNode(a)[0];
        nodes(0, a) = t;
        nodes(1, a) = Real(0.2) * t * (Real(1) - t);
        nodes(2, a) = 0;
      }
      mesh.setPolytopeTransformation({1, 0},
        new ParametricTransformation<Variational::RealH1Element<2>>(
          std::move(nodes), element));
      Location::AABB locator(mesh);
      const Math::SpatialPoint x{0.75, 0.0375, 0};
      if (!locator.locate(x))
      {
        state.SkipWithError("AABB missed the embedded curved query");
        return;
      }
      for (auto _ : state)
        benchmark::DoNotOptimize(locator.locate(x).has_value());
    }

    void BM_AABBGradedHit(benchmark::State& state)
    {
      constexpr size_t CellCount = 1024;
      constexpr Real SmallCellLength = 1e-6;
      constexpr Real DistantVertex = 1e6;
      MeshType::Builder builder;
      builder.initialize(1).nodes(CellCount + 2);
      for (size_t i = 0; i <= CellCount; ++i)
        builder.vertex({static_cast<Real>(i) * SmallCellLength});
      builder.vertex({DistantVertex});
      for (size_t i = 0; i < CellCount; ++i)
        builder.polytope(G::Segment, {i, i + 1});
      const MeshType mesh = builder.finalize();
      Location::AABB locator(mesh);
      const Math::SpatialPoint x{(CellCount - Real(0.25)) * SmallCellLength};
      benchmark::DoNotOptimize(locator.locate(x).has_value());
      for (auto _ : state)
        benchmark::DoNotOptimize(locator.locate(x).has_value());
    }

    void BM_AABBTriangleNarrowMiss(benchmark::State& state)
    {
      const MeshType mesh = MeshType::Builder()
                              .initialize(2)
                              .nodes(3)
                              .vertex({0, 0})
                              .vertex({1, 0})
                              .vertex({0, 1})
                              .polytope(G::Triangle, {0, 1, 2})
                              .finalize();
      Location::AABB locator(mesh);
      const Math::SpatialPoint x{0.8, 0.8};
      benchmark::DoNotOptimize(locator.locate(x).has_value());
      for (auto _ : state)
        benchmark::DoNotOptimize(locator.locate(x).has_value());
    }
  }

/// @cond RODIN_TEST_INTERNAL
#define RODIN_AABB_BENCHMARKS(GEOMETRY)                                                  \
  BENCHMARK_CAPTURE(BM_AABBBuild, GEOMETRY, G::GEOMETRY);                                \
  BENCHMARK_CAPTURE(BM_AABBHit, GEOMETRY, G::GEOMETRY);                                  \
  BENCHMARK_CAPTURE(BM_AABBOuterMiss, GEOMETRY, G::GEOMETRY)

  RODIN_AABB_BENCHMARKS(Point);
  RODIN_AABB_BENCHMARKS(Segment);
  RODIN_AABB_BENCHMARKS(Triangle);
  RODIN_AABB_BENCHMARKS(Quadrilateral);
  RODIN_AABB_BENCHMARKS(Tetrahedron);
  RODIN_AABB_BENCHMARKS(Hexahedron);
  RODIN_AABB_BENCHMARKS(Pyramid);
  RODIN_AABB_BENCHMARKS(Wedge);

#undef RODIN_AABB_BENCHMARKS

  BENCHMARK_CAPTURE(BM_AABBCurvedHit, Segment, G::Segment);
  BENCHMARK_CAPTURE(BM_AABBCurvedHit, Triangle, G::Triangle);
  BENCHMARK_CAPTURE(BM_AABBCurvedHit, Quadrilateral, G::Quadrilateral);
  BENCHMARK_CAPTURE(BM_AABBCurvedHit, Tetrahedron, G::Tetrahedron);
  BENCHMARK_CAPTURE(BM_AABBCurvedHit, Hexahedron, G::Hexahedron);
  BENCHMARK_CAPTURE(BM_AABBCurvedHit, Pyramid, G::Pyramid);
  BENCHMARK_CAPTURE(BM_AABBCurvedHit, Wedge, G::Wedge);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<2>, P2Segment, G::Segment);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<4>, P4Segment, G::Segment);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<2>, P2Triangle, G::Triangle);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<4>, P4Triangle, G::Triangle);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<2>, P2Quadrilateral, G::Quadrilateral);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<4>, P4Quadrilateral, G::Quadrilateral);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<2>, P2Tetrahedron, G::Tetrahedron);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<4>, P4Tetrahedron, G::Tetrahedron);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<2>, P2Hexahedron, G::Hexahedron);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<4>, P4Hexahedron, G::Hexahedron);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<2>, P2Pyramid, G::Pyramid);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<4>, P4Pyramid, G::Pyramid);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<2>, P2Wedge, G::Wedge);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<4>, P4Wedge, G::Wedge);
  // Run only FirstP4 cases in a fresh process to include the conversion-cache
  // construction. Running the other build cases first warms that cache.
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<4>, FirstP4Segment, G::Segment)
    ->Iterations(1)
    ->Repetitions(1);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<4>, FirstP4Triangle, G::Triangle)
    ->Iterations(1)
    ->Repetitions(1);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<4>, FirstP4Quadrilateral, G::Quadrilateral)
    ->Iterations(1)
    ->Repetitions(1);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<4>, FirstP4Tetrahedron, G::Tetrahedron)
    ->Iterations(1)
    ->Repetitions(1);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<4>, FirstP4Hexahedron, G::Hexahedron)
    ->Iterations(1)
    ->Repetitions(1);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<4>, FirstP4Pyramid, G::Pyramid)
    ->Iterations(1)
    ->Repetitions(1);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<4>, FirstP4Wedge, G::Wedge)
    ->Iterations(1)
    ->Repetitions(1);
  BENCHMARK_CAPTURE(BM_AABBBuildWithProjections, Point, G::Point);
  BENCHMARK_CAPTURE(BM_AABBBuildWithProjections, Segment, G::Segment);
  BENCHMARK_CAPTURE(BM_AABBBuildWithProjections, Triangle, G::Triangle);
  BENCHMARK_CAPTURE(BM_AABBBuildWithProjections, Quadrilateral, G::Quadrilateral);
  BENCHMARK_CAPTURE(BM_AABBBuildWithProjections, Tetrahedron, G::Tetrahedron);
  BENCHMARK_CAPTURE(BM_AABBBuildWithProjections, Hexahedron, G::Hexahedron);
  BENCHMARK_CAPTURE(BM_AABBBuildWithProjections, Pyramid, G::Pyramid);
  BENCHMARK_CAPTURE(BM_AABBBuildWithProjections, Wedge, G::Wedge);
  BENCHMARK_CAPTURE(BM_AABBHit, PrunedPoint, G::Point, true);
  BENCHMARK_CAPTURE(BM_AABBHit, PrunedSegment, G::Segment, true);
  BENCHMARK_CAPTURE(BM_AABBCurvedHit, PrunedSegment, G::Segment, true);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<2>, PrunedP2Segment, G::Segment, true);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<4>, PrunedP4Segment, G::Segment, true);
  BENCHMARK_CAPTURE(BM_AABBHit, PrunedTriangle, G::Triangle, true);
  BENCHMARK_CAPTURE(BM_AABBCurvedHit, PrunedTriangle, G::Triangle, true);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<2>, PrunedP2Triangle, G::Triangle, true);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<4>, PrunedP4Triangle, G::Triangle, true);
  BENCHMARK_CAPTURE(BM_AABBHit, PrunedQuadrilateral, G::Quadrilateral, true);
  BENCHMARK_CAPTURE(BM_AABBCurvedHit, PrunedQuadrilateral, G::Quadrilateral, true);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<2>, PrunedP2Quadrilateral, G::Quadrilateral, true);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<4>, PrunedP4Quadrilateral, G::Quadrilateral, true);
  BENCHMARK_CAPTURE(BM_AABBHit, PrunedTetrahedron, G::Tetrahedron, true);
  BENCHMARK_CAPTURE(BM_AABBCurvedHit, PrunedTetrahedron, G::Tetrahedron, true);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<2>, PrunedP2Tetrahedron, G::Tetrahedron, true);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<4>, PrunedP4Tetrahedron, G::Tetrahedron, true);
  BENCHMARK_CAPTURE(BM_AABBHit, PrunedHexahedron, G::Hexahedron, true);
  BENCHMARK_CAPTURE(BM_AABBCurvedHit, PrunedHexahedron, G::Hexahedron, true);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<2>, PrunedP2Hexahedron, G::Hexahedron, true);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<4>, PrunedP4Hexahedron, G::Hexahedron, true);
  BENCHMARK_CAPTURE(BM_AABBHit, PrunedPyramid, G::Pyramid, true);
  BENCHMARK_CAPTURE(BM_AABBCurvedHit, PrunedPyramid, G::Pyramid, true);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<2>, PrunedP2Pyramid, G::Pyramid, true);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<4>, PrunedP4Pyramid, G::Pyramid, true);
  BENCHMARK_CAPTURE(BM_AABBHit, PrunedWedge, G::Wedge, true);
  BENCHMARK_CAPTURE(BM_AABBCurvedHit, PrunedWedge, G::Wedge, true);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<2>, PrunedP2Wedge, G::Wedge, true);
  BENCHMARK_CAPTURE(BM_AABBCurvedBuild<4>, PrunedP4Wedge, G::Wedge, true);
  BENCHMARK(BM_AABBEmbeddedCurvedHit);
  BENCHMARK(BM_AABBGradedHit);
  BENCHMARK(BM_AABBTriangleNarrowMiss);
  /// @endcond
}
