/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <benchmark/benchmark.h>

#include <boost/filesystem.hpp>

#include <Rodin/Geometry.h>
#include <Rodin/Configure.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace RodinBenchmark
{
  struct P1Benchmark : public benchmark::Fixture
  {
    public:
      void SetUp(const benchmark::State&)
      {
        square2d = Mesh<Rodin::Context::Serial>::Builder()
        .initialize(2)
        .nodes(4)
        .vertex({0, 0})
        .vertex({1, 0})
        .vertex({0, 1})
        .vertex({1, 1})
        .polytope(Polytope::Type::Triangle, {0, 1, 2})
        .polytope(Polytope::Type::Triangle, {1, 3, 2})
        .finalize();

        uniformTriangularMesh16 =
          SerialMesh::UniformGrid(Polytope::Type::Triangle, 16, 16);

        uniformTriangularMesh32 =
          SerialMesh::UniformGrid(Polytope::Type::Triangle, 32, 32);

        uniformTriangularMesh32 =
          SerialMesh::UniformGrid(Polytope::Type::Triangle, 32, 32);

        uniformTriangularMesh64 =
          SerialMesh::UniformGrid(Polytope::Type::Triangle, 64, 64);

        uniformTriangularMesh128 =
          SerialMesh::UniformGrid(Polytope::Type::Triangle, 128, 128);
      }

      void TearDown(const benchmark::State&)
      {}

      SerialMesh square2d;
      SerialMesh uniformTriangularMesh16;
      SerialMesh uniformTriangularMesh32;
      SerialMesh uniformTriangularMesh64;
      SerialMesh uniformTriangularMesh128;
  };

  BENCHMARK_F(P1Benchmark, UniformTriangular16_Build)
  (benchmark::State& st)
  {
    for (auto _ : st)
    {
      P1 fes(uniformTriangularMesh16);
      benchmark::DoNotOptimize(fes);
    }
  }

  BENCHMARK_F(P1Benchmark, UniformTriangular32_Build)
  (benchmark::State& st)
  {
    for (auto _ : st)
    {
      P1 fes(uniformTriangularMesh32);
      benchmark::DoNotOptimize(fes);
    }
  }

  BENCHMARK_F(P1Benchmark, UniformTriangular64_Build)
  (benchmark::State& st)
  {
    for (auto _ : st)
    {
      P1 fes(uniformTriangularMesh64);
      benchmark::DoNotOptimize(fes);
    }
  }

  BENCHMARK_F(P1Benchmark, UniformTriangular128_Build)
  (benchmark::State& st)
  {
    for (auto _ : st)
    {
      P1 fes(uniformTriangularMesh128);
      benchmark::DoNotOptimize(fes);
    }
  }

  BENCHMARK_F(P1Benchmark, 2D_Square_GridFunction_Projection_Scalar_SumOfComponents)
  (benchmark::State& st)
  {
    P1 fes(square2d);
    GridFunction gf(fes);
    ScalarFunction c([](const Geometry::Point& p) { return p.x() + p.y(); } );
    for (auto _ : st)
      gf.project(c);
  }

  BENCHMARK_F(P1Benchmark, UniformTriangular16_GridFunction_Projection_Scalar_SumOfComponents)
  (benchmark::State& st)
  {
    P1 fes(uniformTriangularMesh16);
    GridFunction gf(fes);
    ScalarFunction c([](const Geometry::Point& p) { return p.x() + p.y(); } );
    for (auto _ : st)
      gf.project(c);
  }

  BENCHMARK_F(P1Benchmark, UniformTriangular32_GridFunction_Projection_Scalar_SumOfComponents)
  (benchmark::State& st)
  {
    P1 fes(uniformTriangularMesh32);
    GridFunction gf(fes);
    ScalarFunction c([](const Geometry::Point& p) { return p.x() + p.y(); } );
    for (auto _ : st)
      gf.project(c);
  }

  BENCHMARK_F(P1Benchmark, 2D_Square_GridFunction_Projection_Vector_Components)
  (benchmark::State& st)
  {
    constexpr size_t vdim = 2;
    P1 fes(square2d, vdim);
    GridFunction gf(fes);
    VectorFunction c{
      [](const Geometry::Point& p) { return p.x(); },
      [](const Geometry::Point& p) { return p.y(); },
    };
    for (auto _ : st)
      gf.project(c);
  }

  BENCHMARK_F(P1Benchmark, UniformTriangular16_GridFunction_Projection_Vector_Components)
  (benchmark::State& st)
  {
    constexpr size_t vdim = 2;
    P1 fes(uniformTriangularMesh16, vdim);
    GridFunction gf(fes);
    VectorFunction c{
      [](const Geometry::Point& p) { return p.x(); },
      [](const Geometry::Point& p) { return p.y(); },
    };
    for (auto _ : st)
      gf.project(c);
  }

  BENCHMARK_F(P1Benchmark, UniformTriangular32_GridFunction_Projection_Vector_Components)
  (benchmark::State& st)
  {
    constexpr size_t vdim = 2;
    P1 fes(uniformTriangularMesh32, vdim);
    GridFunction gf(fes);
    VectorFunction c{
      [](const Geometry::Point& p) { return p.x(); },
      [](const Geometry::Point& p) { return p.y(); },
    };
    for (auto _ : st)
      gf.project(c);
  }
}