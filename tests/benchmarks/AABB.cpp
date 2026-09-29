/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief AABB build and point-location benchmarks.
 *
 * Mesh construction and query generation are outside the timed loop. Hit
 * queries use non-centroid reference coordinates so inversion is measured.
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

    std::vector<Math::SpatialPoint> mappedQueries(
      const MeshType& mesh, const Math::SpatialPoint& reference)
    {
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

    void BM_AABB2DBuild(benchmark::State& state)
    {
      const MeshType mesh = MeshType::UniformGrid(Polytope::Type::Triangle, {64, 64});
      const auto queries = mappedQueries(mesh, Math::SpatialPoint{0.2, 0.3});
      for (auto _ : state)
      {
        Location::AABB locator(mesh);
        benchmark::DoNotOptimize(locator.locate(queries.front()).has_value());
      }
      state.counters["cells"] = static_cast<double>(mesh.getCellCount());
    }

    void BM_AABB2DHit(benchmark::State& state)
    {
      const MeshType mesh = MeshType::UniformGrid(Polytope::Type::Triangle, {64, 64});
      const auto queries = mappedQueries(mesh, Math::SpatialPoint{0.2, 0.3});
      Location::AABB locator(mesh);
      for (const auto& x : queries)
      {
        if (!locator.locate(x))
        {
          state.SkipWithError("AABB missed a mapped 2D query");
          return;
        }
      }
      size_t i = 0;
      for (auto _ : state)
      {
        benchmark::DoNotOptimize(locator.locate(queries[i]).has_value());
        i = (i + 1) % queries.size();
      }
      state.counters["cells"] = static_cast<double>(mesh.getCellCount());
    }

    void BM_AABB3DHit(benchmark::State& state)
    {
      const MeshType mesh = MeshType::UniformGrid(Polytope::Type::Tetrahedron, {12, 12, 12});
      const auto queries = mappedQueries(mesh, Math::SpatialPoint{0.2, 0.3, 0.1});
      Location::AABB locator(mesh);
      for (const auto& x : queries)
      {
        if (!locator.locate(x))
        {
          state.SkipWithError("AABB missed a mapped 3D query");
          return;
        }
      }
      size_t i = 0;
      for (auto _ : state)
      {
        benchmark::DoNotOptimize(locator.locate(queries[i]).has_value());
        i = (i + 1) % queries.size();
      }
      state.counters["cells"] = static_cast<double>(mesh.getCellCount());
    }

    void BM_AABB2DNarrowMiss(benchmark::State& state)
    {
      const MeshType mesh = MeshType::Builder()
                              .initialize(2)
                              .nodes(3)
                              .vertex({0, 0})
                              .vertex({1, 0})
                              .vertex({0, 1})
                              .polytope(Polytope::Type::Triangle, {0, 1, 2})
                              .finalize();
      Location::AABB locator(mesh);
      const Math::SpatialPoint x{0.8, 0.8};
      benchmark::DoNotOptimize(locator.locate(x).has_value());
      for (auto _ : state)
        benchmark::DoNotOptimize(locator.locate(x).has_value());
    }

    void BM_AABB2DCurvedHit(benchmark::State& state)
    {
      MeshType mesh = MeshType::UniformGrid(Polytope::Type::Triangle, {16, 16});
      Variational::RealH1Element<2> element(Polytope::Type::Triangle);
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        PointCloud nodes(2, element.getCount());
        for (size_t a = 0; a < element.getCount(); ++a)
        {
          Math::SpatialPoint x;
          cell->getTransformation().transform(x, element.getNode(a));
          nodes(0, a) = x[0];
          nodes(1, a) = x[1] + 0.15 * x[0] * x[0];
        }
        mesh.setPolytopeTransformation({2, cell->getIndex()},
          new ParametricTransformation<Variational::RealH1Element<2>>(
            std::move(nodes), element));
      }
      const auto queries = mappedQueries(mesh, Math::SpatialPoint{0.2, 0.3});
      Location::AABB locator(mesh);
      for (const auto& x : queries)
      {
        if (!locator.locate(x))
        {
          state.SkipWithError("AABB missed a mapped curved query");
          return;
        }
      }
      size_t i = 0;
      for (auto _ : state)
      {
        benchmark::DoNotOptimize(locator.locate(queries[i]).has_value());
        i = (i + 1) % queries.size();
      }
      state.counters["cells"] = static_cast<double>(mesh.getCellCount());
    }
  }

  BENCHMARK(BM_AABB2DBuild);
  BENCHMARK(BM_AABB2DHit);
  BENCHMARK(BM_AABB3DHit);
  BENCHMARK(BM_AABB2DNarrowMiss);
  BENCHMARK(BM_AABB2DCurvedHit);
}
