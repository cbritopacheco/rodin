/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Grid-function point-evaluation benchmarks.
 *
 * Compares finite-element expansion evaluation with the mapped-basis
 * contraction previously used by GridFunction at arbitrary geometric points.
 */

#include <benchmark/benchmark.h>

#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Variational/H1.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Benchmarks
{
  struct GridFunctionEvaluationBenchmark : public benchmark::Fixture
  {
    void SetUp(const benchmark::State&)
    {
      mesh = LocalMesh::UniformGrid(Polytope::Type::Triangle, { 4, 4 });
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
    }

    LocalMesh mesh;
  };

  BENCHMARK_F(GridFunctionEvaluationBenchmark, P1VectorExpansion)
  (benchmark::State& state)
  {
    P1 fes(mesh, 2);
    GridFunction gf(fes);
    gf = VectorFunction{
      [](const Point& p) { return 1.0 + p.x() - p.y(); },
      [](const Point& p) { return -2.0 + 2.0 * p.x() + p.y(); }
    };
    const auto cell = mesh.getCell(0);
    const Point p(*cell, Math::SpatialPoint{ 0.2, 0.3 });

    for (auto _ : state)
    {
      Math::SpatialVector<Real> value;
      gf.interpolate(value, p);
      benchmark::DoNotOptimize(value);
    }
  }

  BENCHMARK_F(GridFunctionEvaluationBenchmark, P1ScalarExpansion)
  (benchmark::State& state)
  {
    P1 fes(mesh);
    GridFunction gf(fes);
    gf = RealFunction(
        [](const Point& p) { return 1.0 + p.x() - p.y(); });
    const auto cell = mesh.getCell(0);
    const Point p(*cell, Math::SpatialPoint{ 0.2, 0.3 });

    for (auto _ : state)
    {
      Real value;
      gf.interpolate(value, p);
      benchmark::DoNotOptimize(value);
    }
  }

  BENCHMARK_F(GridFunctionEvaluationBenchmark, P1ScalarReferenceExpansion)
  (benchmark::State& state)
  {
    P1 fes(mesh);
    GridFunction gf(fes);
    gf = RealFunction([](const Point& p) { return 1.0 + p.x() - p.y(); });
    const auto cell = mesh.getCell(0);
    const Point p(*cell, Math::SpatialPoint{0.2, 0.3});
    const auto& fe = fes.getFiniteElement(2, 0);
    const auto& dofs = fes.getDOFs(2, 0);

    for (auto _ : state)
    {
      Real value;
      fe.evaluate(
        value, [&](size_t local) -> decltype(auto) { return gf[dofs[local]]; },
        p.getReferenceCoordinates());
      benchmark::DoNotOptimize(value);
    }
  }

  BENCHMARK_F(GridFunctionEvaluationBenchmark, P1ScalarSpaceExpansion)
  (benchmark::State& state)
  {
    P1 fes(mesh);
    GridFunction gf(fes);
    gf = RealFunction([](const Point& p) { return 1.0 + p.x() - p.y(); });
    const auto cell = mesh.getCell(0);
    const Point p(*cell, Math::SpatialPoint{0.2, 0.3});
    const auto& dofs = fes.getDOFs(2, 0);

    for (auto _ : state)
    {
      Real value;
      fes.evaluate(
        value, {2, 0}, [&](size_t local) -> decltype(auto) { return gf[dofs[local]]; },
        p);
      benchmark::DoNotOptimize(value);
    }
  }

  BENCHMARK_F(GridFunctionEvaluationBenchmark, P1ScalarMappedBasis)
  (benchmark::State& state)
  {
    P1 fes(mesh);
    GridFunction gf(fes);
    gf = RealFunction(
        [](const Point& p) { return 1.0 + p.x() - p.y(); });
    const auto cell = mesh.getCell(0);
    const Point p(*cell, Math::SpatialPoint{ 0.2, 0.3 });
    const auto& fe = fes.getFiniteElement(2, 0);
    const auto& dofs = fes.getDOFs(2, 0);

    for (auto _ : state)
    {
      Real value = 0;
      for (size_t local = 0; local < fe.getCount(); ++local)
      {
        const auto mapping = fes.getPushforward({ 2, 0 }, fe.getBasis(local));
        value += gf[dofs[local]] * mapping(p);
      }
      benchmark::DoNotOptimize(value);
    }
  }

  BENCHMARK_F(GridFunctionEvaluationBenchmark, P1VectorMappedBasis)
  (benchmark::State& state)
  {
    P1 fes(mesh, 2);
    GridFunction gf(fes);
    gf = VectorFunction{
      [](const Point& p) { return 1.0 + p.x() - p.y(); },
      [](const Point& p) { return -2.0 + 2.0 * p.x() + p.y(); }
    };
    const auto cell = mesh.getCell(0);
    const Point p(*cell, Math::SpatialPoint{ 0.2, 0.3 });
    const auto& fe = fes.getFiniteElement(2, 0);
    const auto& dofs = fes.getDOFs(2, 0);

    for (auto _ : state)
    {
      Math::SpatialVector<Real> value;
      for (size_t local = 0; local < fe.getCount(); ++local)
      {
        const auto mapping = fes.getPushforward({ 2, 0 }, fe.getBasis(local));
        const auto term = gf[dofs[local]] * mapping(p);
        if (local == 0)
          value = term;
        else
          value += term;
      }
      benchmark::DoNotOptimize(value);
    }
  }

  BENCHMARK_F(GridFunctionEvaluationBenchmark, H1P2VectorExpansion)
  (benchmark::State& state)
  {
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, size_t(2));
    GridFunction gf(fes);
    gf = VectorFunction{
      [](const Point& p) { return 1.0 + p.x() - p.y(); },
      [](const Point& p) { return -2.0 + 2.0 * p.x() + p.y(); }
    };
    const auto cell = mesh.getCell(0);
    const Point p(*cell, Math::SpatialPoint{ 0.2, 0.3 });

    for (auto _ : state)
    {
      Math::SpatialVector<Real> value;
      gf.interpolate(value, p);
      benchmark::DoNotOptimize(value);
    }
  }

  BENCHMARK_F(GridFunctionEvaluationBenchmark, H1P2ScalarExpansion)
  (benchmark::State& state)
  {
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);
    gf = RealFunction(
        [](const Point& p) { return 1.0 + p.x() - p.y(); });
    const auto cell = mesh.getCell(0);
    const Point p(*cell, Math::SpatialPoint{ 0.2, 0.3 });

    for (auto _ : state)
    {
      Real value;
      gf.interpolate(value, p);
      benchmark::DoNotOptimize(value);
    }
  }

  BENCHMARK_F(GridFunctionEvaluationBenchmark, H1P2ScalarMappedBasis)
  (benchmark::State& state)
  {
    H1 fes(std::integral_constant<size_t, 2>{}, mesh);
    GridFunction gf(fes);
    gf = RealFunction(
        [](const Point& p) { return 1.0 + p.x() - p.y(); });
    const auto cell = mesh.getCell(0);
    const Point p(*cell, Math::SpatialPoint{ 0.2, 0.3 });
    const auto& fe = fes.getFiniteElement(2, 0);
    const auto& dofs = fes.getDOFs(2, 0);

    for (auto _ : state)
    {
      Real value = 0;
      for (size_t local = 0; local < fe.getCount(); ++local)
      {
        const auto mapping = fes.getPushforward({ 2, 0 }, fe.getBasis(local));
        value += gf[dofs[local]] * mapping(p);
      }
      benchmark::DoNotOptimize(value);
    }
  }

  BENCHMARK_F(GridFunctionEvaluationBenchmark, H1P2VectorMappedBasis)
  (benchmark::State& state)
  {
    H1 fes(std::integral_constant<size_t, 2>{}, mesh, size_t(2));
    GridFunction gf(fes);
    gf = VectorFunction{
      [](const Point& p) { return 1.0 + p.x() - p.y(); },
      [](const Point& p) { return -2.0 + 2.0 * p.x() + p.y(); }
    };
    const auto cell = mesh.getCell(0);
    const Point p(*cell, Math::SpatialPoint{ 0.2, 0.3 });
    const auto& fe = fes.getFiniteElement(2, 0);
    const auto& dofs = fes.getDOFs(2, 0);

    for (auto _ : state)
    {
      Math::SpatialVector<Real> value;
      for (size_t local = 0; local < fe.getCount(); ++local)
      {
        const auto mapping = fes.getPushforward({ 2, 0 }, fe.getBasis(local));
        const auto term = gf[dofs[local]] * mapping(p);
        if (local == 0)
          value = term;
        else
          value += term;
      }
      benchmark::DoNotOptimize(value);
    }
  }
}
