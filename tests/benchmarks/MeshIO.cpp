/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <benchmark/benchmark.h>

#include <Rodin/Geometry.h>
#include <Rodin/Configure.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace RodinBenchmark
{
  struct MeshIO : public benchmark::Fixture
  {
    public:
      void SetUp(const benchmark::State&)
      {}

      void TearDown(const benchmark::State&)
      {}
  };

  BENCHMARK_F(MeshIO, Load_MEDIT_2D_Square)(benchmark::State& st)
  {
    static constexpr const char* filename = "mmg/Square.medit.mesh";
    boost::filesystem::path meshfile;
    meshfile = boost::filesystem::path(RODIN_RESOURCES_DIR);
    meshfile.append(filename);
    Mesh mesh;
    for (auto _ : st)
      mesh.load(meshfile, IO::FileFormat::MEDIT);
  }

  BENCHMARK_F(MeshIO, Load_MEDIT_2D_UniformTriangular64)(benchmark::State& st)
  {
    static constexpr const char* filename = "mmg/UniformTriangular64.medit.mesh";
    boost::filesystem::path meshfile;
    meshfile = boost::filesystem::path(RODIN_RESOURCES_DIR);
    meshfile.append(filename);
    Mesh mesh;
    for (auto _ : st)
      mesh.load(meshfile, IO::FileFormat::MEDIT);
  }
}
