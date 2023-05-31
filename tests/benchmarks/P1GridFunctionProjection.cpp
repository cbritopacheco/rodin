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
  struct ScalarP1GridFunctionProjection : public benchmark::Fixture
  {
    public:
      void SetUp(const benchmark::State&)
      {}

      void TearDown(const benchmark::State&)
      {}
  };

  BENCHMARK_F(ScalarP1GridFunctionProjection, 2D_Square_SumOfComponents)
  (benchmark::State& st)
  {
    constexpr size_t sdim = 2;
    Mesh mesh =
      Mesh<Rodin::Context::Serial>::Builder()
      .initialize(sdim)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Geometry::Triangle, {0, 1, 2})
      .polytope(Polytope::Geometry::Triangle, {1, 3, 2})
      .finalize();
    P1 fes(mesh);
    GridFunction gf(fes);
    ScalarFunction c([](const Geometry::Point& p) { return p.x() + p.y(); } );
    for (auto _ : st)
      gf.project(c);
  }
}

namespace RodinBenchmark
{
  struct VectorP1GridFunctionProjection : public benchmark::Fixture
  {
    public:
      void SetUp(const benchmark::State&)
      {}

      void TearDown(const benchmark::State&)
      {}
  };

  BENCHMARK_F(VectorP1GridFunctionProjection, 2D_Square_Components)
  (benchmark::State& st)
  {
    constexpr size_t sdim = 2;
    constexpr size_t vdim = sdim;
    Mesh mesh =
      Mesh<Rodin::Context::Serial>::Builder()
      .initialize(sdim)
      .nodes(4)
      .vertex({0, 0})
      .vertex({1, 0})
      .vertex({0, 1})
      .vertex({1, 1})
      .polytope(Polytope::Geometry::Triangle, {0, 1, 2})
      .polytope(Polytope::Geometry::Triangle, {1, 3, 2})
      .finalize();
    P1 fes(mesh, vdim);
    GridFunction gf(fes);
    VectorFunction c{
      [](const Geometry::Point& p) { return p.x(); },
      [](const Geometry::Point& p) { return p.y(); },
    };
    for (auto _ : st)
      gf.project(c);
  }
}

