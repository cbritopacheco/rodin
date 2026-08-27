/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Quadrature rule benchmark coverage.
 *
 * Establishes the baseline against which a change of quadrature family is
 * judged. Point count alone is not the figure of merit: a rule with fewer
 * points is only faster if the per-point work dominates rule access, so both
 * are measured, and the per-point workload is representative of what an
 * element integrator actually does rather than an empty loop.
 */

#include <benchmark/benchmark.h>

#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include <Rodin/Geometry/Polytope.h>
#include <Rodin/Math.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Benchmarks
{
  namespace
  {
    Polytope::Type geometryOf(int64_t i)
    {
      switch (i)
      {
        case 0:
          return Polytope::Type::Triangle;
        case 1:
          return Polytope::Type::Tetrahedron;
        case 2:
          return Polytope::Type::Wedge;
        default:
          return Polytope::Type::Pyramid;
      }
    }
  }

  /// @brief Cost of obtaining a rule.
  ///
  /// This measures the pooled lookup, not construction. Read it with care: it
  /// is not a per-element-type cost. PolytopeQuadratureFormula holds a
  /// thread-local cache of eight entries, inserted round-robin and scanned
  /// linearly from slot zero, so the cost of a hit is proportional to the slot
  /// the key landed in, hence to the order in which distinct keys were first
  /// requested in the process. Running this family in full gives 3.5 ns for
  /// the first case and 32 ns for the eighth, resetting every eight; running
  /// any single case alone gives 3.5 ns whatever the element. The ramp is the
  /// cache, not the geometry.
  ///
  /// It is kept because the envelope is the useful quantity: 3.5 to 32 ns
  /// against a 375 ns degree-4 wedge sweep bounds the lookup at under 9% of an
  /// integrator's time in the worst slot, so no plausible change of quadrature
  /// family is limited by it.
  static void BM_QuadratureAccess(benchmark::State& state)
  {
    const auto g = geometryOf(state.range(0));
    const auto order = static_cast<size_t>(state.range(1));
    for (auto _ : state)
    {
      const auto& qf = QF::PolytopeQuadratureFormula::get(order, g);
      benchmark::DoNotOptimize(qf.getSize());
    }
    state.counters["points"] =
      static_cast<double>(QF::PolytopeQuadratureFormula::get(order, g).getSize());
  }

  /// @brief Cost of one sweep over the rule doing per-point work of the kind
  /// an element integrator does: build a deformation gradient from the point,
  /// take its determinant and its cofactor contraction. This is the quantity
  /// a reduction in point count actually reduces.
  static void BM_QuadratureSweep(benchmark::State& state)
  {
    const auto g = geometryOf(state.range(0));
    const auto order = static_cast<size_t>(state.range(1));
    const auto& qf = QF::PolytopeQuadratureFormula::get(order, g);
    const size_t d = Polytope::Traits(g).getDimension();

    for (auto _ : state)
    {
      Real acc = 0;
      for (size_t q = 0; q < qf.getSize(); ++q)
      {
        const auto& x = qf.getPoint(q);
        Math::SpatialMatrix<Real> F(d, d);
        for (size_t i = 0; i < d; ++i)
          for (size_t j = 0; j < d; ++j)
            F(i, j) = (i == j ? Real(1) : Real(0)) +
              Real(0.1) * x[static_cast<Eigen::Index>((i + j) % d)];
        const Real j0 = F.determinant();
        acc += qf.getWeight(q) * (j0 + F.squaredNorm() / j0);
      }
      benchmark::DoNotOptimize(acc);
    }
    state.counters["points"] = static_cast<double>(qf.getSize());
    state.counters["ns_per_point"] = benchmark::Counter(static_cast<double>(qf.getSize()),
      benchmark::Counter::kIsIterationInvariantRate | benchmark::Counter::kInvert);
  }

  static void Cases(benchmark::internal::Benchmark* b)
  {
    for (int64_t g = 0; g < 4; ++g)
      for (int64_t order : {2, 4, 6, 8})
        b->Args({g, order});
  }

  /// @brief Registers point-access benchmarks for the quadrature cases.
  BENCHMARK(BM_QuadratureAccess)->Apply(Cases);

  /// @brief Registers integrand-sweep benchmarks for the quadrature cases.
  BENCHMARK(BM_QuadratureSweep)->Apply(Cases);
}
