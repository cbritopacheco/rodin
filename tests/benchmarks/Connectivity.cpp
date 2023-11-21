/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2023.
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
  struct Connectivity : public benchmark::Fixture
  {
    public:
      template <class Mesh>
      void test(Mesh& mesh, benchmark::State& st)
      {
        const size_t D = mesh.getDimension();
        for (auto _ : st)
        {
          for (size_t d = 0; d <= D; d++)
            for (size_t dp = 0; dp <= D; dp++)
              mesh.getConnectivity().compute(d, dp);
        }
      }

      void SetUp(const benchmark::State&)
      {}

      void TearDown(const benchmark::State&)
      {}
  };


  BENCHMARK_F(Connectivity, Triangular_16x16)(benchmark::State& st)
  {
    auto mesh = SerialMesh::UniformGrid(Polytope::Type::Triangle, 16, 16);
    test(mesh, st);
  }

  BENCHMARK_F(Connectivity, Triangular_32x32)(benchmark::State& st)
  {
    auto mesh = SerialMesh::UniformGrid(Polytope::Type::Triangle, 32, 32);
    test(mesh, st);
  }

  BENCHMARK_F(Connectivity, Triangular_64x64)(benchmark::State& st)
  {
    auto mesh = SerialMesh::UniformGrid(Polytope::Type::Triangle, 64, 64);
    test(mesh, st);
  }
}

