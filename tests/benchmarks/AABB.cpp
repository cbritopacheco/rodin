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

    void curveMesh(MeshType& mesh, G type)
    {
      const size_t dimension = mesh.getDimension();
      Variational::RealH1Element<2> element(type);
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
          new ParametricTransformation<Variational::RealH1Element<2>>(
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

    void BM_AABBHit(benchmark::State& state, G type)
    {
      const MeshType mesh = makeMesh(type);
      const auto queries = mappedQueries(mesh, type);
      Location::AABB locator(mesh);
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

    void BM_AABBCurvedHit(benchmark::State& state, G type)
    {
      MeshType mesh = makeMesh(type, true);
      curveMesh(mesh, type);
      const auto queries = mappedQueries(mesh, type);
      Location::AABB locator(mesh);
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
  BENCHMARK(BM_AABBTriangleNarrowMiss);
}
