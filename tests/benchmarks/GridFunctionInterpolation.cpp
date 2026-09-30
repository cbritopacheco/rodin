/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Interpolation and re-interpolation benchmarks with numerical oracles.
 *
 * Meshes, spaces, analytic functions and quadrature points are prepared before
 * timing. Each operation checks the final value against the analytic field.
 * Cache hits repeat one sample; misses alternate two samples in the same cell.
 * Cross-space assignment uses distinct spaces on the same mesh and therefore
 * interpolates the source rather than copying its coefficients.
 */
#include <benchmark/benchmark.h>
#include <cmath>
#include <string>

#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Benchmarks
{
  /// @cond RODIN_TEST_INTERNAL
  template <size_t Order, bool Vector>
  void interpolation(benchmark::State& state)
  {
    const size_t dimension = state.range(0);
    const auto geometry = dimension == 2 ? Polytope::Type::Triangle : Polytope::Type::Tetrahedron;
    LocalMesh mesh = dimension == 2 ? LocalMesh::UniformGrid(geometry, {4, 4}) :
      LocalMesh::UniformGrid(geometry, {4, 4, 4});
    for (size_t d = dimension; d > 0; --d)
      mesh.getConnectivity().compute(d, d - 1);
    const auto makeSpace = [&] {
      if constexpr (Order == 1 && !Vector)
        return P1(mesh);
      else if constexpr (Order == 1 && Vector)
        return P1(mesh, dimension);
      else if constexpr (!Vector)
        return H1(std::integral_constant<size_t, Order>{}, mesh);
      else
        return H1(std::integral_constant<size_t, Order>{}, mesh, dimension);
    };
    auto fes = makeSpace();
    auto sourceSpace = makeSpace();
    GridFunction gf(fes);
    GridFunction source0(sourceSpace);
    GridFunction source1(sourceSpace);
    Real phase = 0;
    const auto scalar = [&](const Point& p) {
      Real value = 1 + phase;
      for (size_t d = 0; d < dimension; ++d)
        value += (d + 1) * p.getCoordinates()(d);
      return value;
    };
    const auto fn = [&] {
      if constexpr (!Vector)
        return RealFunction(scalar);
      else
        return VectorFunction(dimension, [&](const Point& p) {
          Math::SpatialVector<Real> value(dimension);
          for (size_t d = 0; d < dimension; ++d)
            value(d) = (d + 1) * scalar(p);
          return value;
        });
    }();
    gf.project(fn);
    source0.project(fn);
    phase = 1;
    source1.project(fn);
    phase = 0;
    const Polytope cell(dimension, 0, mesh);
    const auto& qf = QF::PolytopeQuadratureFormula::get(2 * Order, geometry);
    const Point point0(cell, qf.getPoint(0));
    const Point point1(cell, qf.getPoint(1));
    const IntegrationPoint ip0(point0, &qf, 0);
    const IntegrationPoint ip1(point1, &qf, 1);
    const Point* lastPoint = &point0;
    typename decltype(gf)::RangeType value{};
    gf.interpolate(value, ip0); // Warm the cache and lazy geometry outside timing.
    gf.interpolate(value, ip1);
    gf.interpolate(value, ip0);
    const auto evaluate = [&](const IntegrationPoint& ip) {
      gf.interpolate(value, ip);
      benchmark::DoNotOptimize(value);
    };
    switch (state.range(1))
    {
      case 0: // Quadrature cache hit.
        for (auto _ : state)
          evaluate(ip0);
        break;
      case 1: // Alternating samples force basis-cache refresh.
        for (auto _ : state)
        {
          phase = 1 - phase;
          const auto& ip = phase == 0 ? ip0 : ip1;
          lastPoint = &ip.getPoint();
          evaluate(ip);
        }
        phase = 0; // The field itself was not changed.
        break;
      case 2: // Pointwise evaluation bypasses quadrature tabulation.
        for (auto _ : state)
        {
          gf.interpolate(value, point0);
          benchmark::DoNotOptimize(value);
        }
        break;
      case 3: // Re-interpolate changing analytic data into an existing field.
        for (auto _ : state)
        {
          phase = 1 - phase;
          gf.project(fn);
          evaluate(ip0);
          benchmark::ClobberMemory();
        }
        break;
      case 4: // Same-order, distinct-space assignment must re-interpolate.
        for (auto _ : state)
        {
          phase = 1 - phase;
          gf = phase == 0 ? source0 : source1;
          evaluate(ip0);
          benchmark::ClobberMemory();
        }
        break;
      case 5: // First interpolation includes field construction/destruction.
        for (auto _ : state)
        {
          GridFunction fresh(fes);
          fresh.project(fn);
          fresh.interpolate(value, ip0);
          benchmark::DoNotOptimize(value);
          benchmark::ClobberMemory();
        }
        break;
    }
    const auto expected = fn(*lastPoint);
    const Real error = [&] {
      if constexpr (Vector)
        return (value - expected).norm();
      else
        return std::abs(value - expected);
    }();
    state.counters["absolute_error"] = error;
    state.counters["dofs"] = fes.getSize();
    if (!std::isfinite(error) || error > 1e-10)
      state.SkipWithError("Interpolation disagrees with the analytic field");
    state.SetItemsProcessed(state.iterations());
  }

  const bool registeredInterpolation = [] {
    const char* operations[] = {"QuadratureHit", "QuadratureMiss", "Pointwise",
      "Reproject", "CrossSpaceAssignment", "FirstInterpolation"};
    for (size_t operation = 0; operation < 6; ++operation)
    {
      const auto args = [&](auto* registration) {
        registration->Args({2, static_cast<int64_t>(operation)})
          ->Args({3, static_cast<int64_t>(operation)})
          ->ArgNames({"dimension", "operation"});
      };
      args(benchmark::RegisterBenchmark((std::string("Interpolation/P1Scalar/") +
        operations[operation]).c_str(), &interpolation<1, false>));
      args(benchmark::RegisterBenchmark((std::string("Interpolation/P1Vector/") +
        operations[operation]).c_str(), &interpolation<1, true>));
      args(benchmark::RegisterBenchmark((std::string("Interpolation/P2Scalar/") +
        operations[operation]).c_str(), &interpolation<2, false>));
      args(benchmark::RegisterBenchmark((std::string("Interpolation/P2Vector/") +
        operations[operation]).c_str(), &interpolation<2, true>));
    }
    return true;
  }();
  /// @endcond
}
