/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Mesh I/O benchmark coverage
 *
 * Benchmarks loading representative MEDIT meshes from disk. The cases protect the cost model for mesh parser setup and repeated mesh-file ingestion used by examples and regression workflows.
 */

#include <benchmark/benchmark.h>

#include <Rodin/Geometry.h>
#include <Rodin/Configure.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Benchmarks
{
  struct MeshIO : public benchmark::Fixture
  {
    public:
      void SetUp(const benchmark::State&)
      {}

      void TearDown(const benchmark::State&)
      {}
  };

  /// @brief Measures repeated MEDIT loading of the square 2D resource mesh.
  BENCHMARK_F(MeshIO, Load_MEDIT_2D_Square)(benchmark::State& st)
  {
    static constexpr const char* filename = "medit/Square.medit.mesh";
    boost::filesystem::path meshfile;
    meshfile = boost::filesystem::path(RODIN_RESOURCES_DIR);
    meshfile.append(filename);
    Mesh mesh;
    for (auto _ : st)
      mesh.load(meshfile, IO::FileFormat::MEDIT);
  }

  /// @brief Measures repeated MEDIT loading of the 64x64 triangular resource mesh.
  BENCHMARK_F(MeshIO, Load_MEDIT_2D_UniformTriangular64)(benchmark::State& st)
  {
    static constexpr const char* filename = "medit/UniformTriangular64.medit.mesh";
    boost::filesystem::path meshfile;
    meshfile = boost::filesystem::path(RODIN_RESOURCES_DIR);
    meshfile.append(filename);
    Mesh mesh;
    for (auto _ : st)
      mesh.load(meshfile, IO::FileFormat::MEDIT);
  }
}
