/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2023.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Uniform-grid benchmark coverage
 *
 * Benchmarks construction of triangular uniform grids at increasing resolutions. The cases measure mesh generation scaling for the geometry layer without adding finite-element or solver work.
 */

#include <benchmark/benchmark.h>

#include <Rodin/Geometry.h>
#include <Rodin/Configure.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Benchmarks
{
  struct UniformGrid : public benchmark::Fixture
  {
    public:
      void SetUp(const benchmark::State&)
      {}

      void TearDown(const benchmark::State&)
      {}
  };

  /// @brief Measures construction of a 16x16 triangular uniform grid.
  BENCHMARK_F(UniformGrid, Triangular_16x16)(benchmark::State& st)
  {
    for (auto _ : st)
      LocalMesh::UniformGrid(Polytope::Type::Triangle, { 16, 16 });
  }

  /// @brief Measures construction of a 64x64 triangular uniform grid.
  BENCHMARK_F(UniformGrid, Triangular_64x64)(benchmark::State& st)
  {
    for (auto _ : st)
      LocalMesh::UniformGrid(Polytope::Type::Triangle, { 64, 64 });
  }

  /// @brief Measures construction of a 128x128 triangular uniform grid.
  BENCHMARK_F(UniformGrid, Triangular_128x128)(benchmark::State& st)
  {
    for (auto _ : st)
      LocalMesh::UniformGrid(Polytope::Type::Triangle, { 128, 128 });
  }

  /// @brief Measures construction of a 256x256 triangular uniform grid.
  BENCHMARK_F(UniformGrid, Triangular_256x256)(benchmark::State& st)
  {
    for (auto _ : st)
      LocalMesh::UniformGrid(Polytope::Type::Triangle, { 256, 256 });
  }

  /// @brief Measures construction of a 512x512 triangular uniform grid.
  BENCHMARK_F(UniformGrid, Triangular_512x512)(benchmark::State& st)
  {
    for (auto _ : st)
      LocalMesh::UniformGrid(Polytope::Type::Triangle, { 512, 512 });
  }
}
