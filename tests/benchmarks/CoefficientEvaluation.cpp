/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief CPU benchmarks of generic coefficient-weighted assembly. */
#include <benchmark/benchmark.h>
#include <cmath>
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  template <size_t Order, class Range>
  void coefficientAssembly(benchmark::State& state)
  {
    const size_t dimension = state.range(0);
    const bool expensive = state.range(1);
    auto mesh = dimension == 2
      ? LocalMesh::UniformGrid(Polytope::Type::Triangle, {3, 3})
      : LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {3, 3, 3});
    for (size_t d = 1; d <= dimension; ++d)
    {
      for (size_t lower = 0; lower < d; ++lower)
        mesh.getConnectivity().compute(d, lower);
    }
    auto fes = [&] {
      if constexpr (FormLanguage::IsMatrixRange<Range>::Value)
        return H1<Order, Range>(std::integral_constant<size_t, Order>{}, mesh, 2, 3);
      else if constexpr (FormLanguage::IsVectorRange<Range>::Value)
        return H1<Order, Range>(std::integral_constant<size_t, Order>{}, mesh, dimension);
      else
        return H1<Order, Range>(std::integral_constant<size_t, Order>{}, mesh);
    }();
    TrialFunction u(fes);
    TestFunction v(fes);
    auto f = RealFunction([&](const Point& point) {
      Real value = 2;
      if (expensive)
        for (size_t i = 1; i <= 64; ++i)
          value += std::sin(point.x() + Real(i) / 10) / Real(i + 1);
      return value;
    });
    // Test-side weighting deliberately exercises the generic rule. The
    // specialized trial-side weighted mass rule provides an independent path.
    auto rule = Integral(u, f * v);
    auto reference = Integral(f * u, v);
    static_assert(!decltype(rule)::Specialized);
    rule.setOrder(2 * Order + 2);
    reference.setOrder(2 * Order + 2);
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      rule.setPolytope(*cell);
      reference.setPolytope(*cell);
      const auto n =
        fes.getFiniteElement(cell->getDimension(), cell->getIndex()).getCount();
      for (size_t i = 0; i < n; ++i)
      {
        for (size_t j = 0; j < n; ++j)
        {
          const Real expected = reference.integrate(i, j);
          if (std::abs(rule.integrate(i, j) - expected) >
            1e-12 * (1 + std::abs(expected)))
          {
            state.SkipWithError("Generic and specialized weighted operators disagree");
            return;
          }
        }
      }
    }
    for (auto _ : state)
    {
      for (auto cell = mesh.getCell(); cell; ++cell)
      {
        rule.setPolytope(*cell);
        benchmark::DoNotOptimize(rule.integrate(0, 0));
      }
      benchmark::ClobberMemory();
    }
    state.SetItemsProcessed(state.iterations() * mesh.getCellCount());
  }

  const bool registered = [] {
    const auto add = []<size_t Order, class Range>(const char* range) {
      const std::string name =
        "GenericCoefficient/" + std::string(range) + "/P" + std::to_string(Order);
      benchmark::RegisterBenchmark(name.c_str(), &coefficientAssembly<Order, Range>)
        ->Args({2, 0})
        ->Args({2, 1})
        ->Args({3, 0})
        ->Args({3, 1})
        ->ArgNames({"dimension", "expensive"})
        ->MeasureProcessCPUTime();
    };
    const auto order = [&]<size_t Order>() {
      add.template operator()<Order, Real>("Scalar");
      add.template operator()<Order, Math::SpatialVector<Real>>("Vector");
      add.template operator()<Order, Math::SpatialMatrix<Real>>("Matrix");
    };
    order.template operator()<1>();
    order.template operator()<2>();
    order.template operator()<3>();
    return true;
  }();
}
