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
  struct P1_GridFunction_Projection : public benchmark::Fixture
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
        .polytope(Polytope::Geometry::Triangle, {0, 1, 2})
        .polytope(Polytope::Geometry::Triangle, {1, 3, 2})
        .finalize();

        uniformTriangularMesh16 =
          SerialMesh::UniformGrid(Polytope::Geometry::Triangle, 16, 16);

        uniformTriangularMesh32 =
          SerialMesh::UniformGrid(Polytope::Geometry::Triangle, 32, 32);
      }

      void TearDown(const benchmark::State&)
      {}

      SerialMesh square2d;
      SerialMesh uniformTriangularMesh16;
      SerialMesh uniformTriangularMesh32;
  };

  BENCHMARK_F(P1_GridFunction_Projection, Scalar_2D_Square_SumOfComponents)
  (benchmark::State& st)
  {
    P1 fes(square2d);
    GridFunction gf(fes);
    ScalarFunction c([](const Geometry::Point& p) { return p.x() + p.y(); } );
    for (auto _ : st)
      gf.project(c);
  }

  BENCHMARK_F(P1_GridFunction_Projection, UniformTriangular16_SumOfComponents)
  (benchmark::State& st)
  {
    P1 fes(uniformTriangularMesh16);
    GridFunction gf(fes);
    ScalarFunction c([](const Geometry::Point& p) { return p.x() + p.y(); } );
    for (auto _ : st)
      gf.project(c);
  }

  BENCHMARK_F(P1_GridFunction_Projection, UniformTriangular32_SumOfComponents)
  (benchmark::State& st)
  {
    P1 fes(uniformTriangularMesh32);
    GridFunction gf(fes);
    ScalarFunction c([](const Geometry::Point& p) { return p.x() + p.y(); } );
    for (auto _ : st)
      gf.project(c);
  }

  BENCHMARK_F(P1_GridFunction_Projection, Vector_2D_Square_Components)
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

  BENCHMARK_F(P1_GridFunction_Projection, Vector_2D_UniformTriangular16_Components)
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

  BENCHMARK_F(P1_GridFunction_Projection, Vector_2D_UniformTriangular32_Components)
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
