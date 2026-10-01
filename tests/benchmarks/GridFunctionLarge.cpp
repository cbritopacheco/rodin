/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <benchmark/benchmark.h>
#include <array>
#include <cmath>
#include <memory>
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
  /**
   * @brief A large, immutable mesh traversal and its analytic interpolation oracle.
   *
   * Mesh, space, field and mapped quadrature points are constructed once per
   * field/dimension. Points and expected sums are prepared outside timing.
   * Traversal workloads visit every cell, so the slow path refreshes DOF and basis
   * caches. Consecutive repetitions of each sample produce controlled hit ratios.
   * PureHit repeats one prewarmed sample of the large field without traversal.
   * The pointwise alternative uses the space's expansion evaluator and retains
   * the DOF cache, but bypasses the quadrature basis cache entirely.
   */
  template <size_t Order, bool Vector>
  class LargeInterpolationData
  {
    public:
      using Range = std::conditional_t<Vector, Math::SpatialVector<Real>, Real>;
      using Space = std::conditional_t<Order == 1, P1<Range>, H1<Order, Range>>;
      using Field = GridFunction<Space, Math::Vector<Real>>;

      explicit LargeInterpolationData(size_t dimension)
        : m_dimension(dimension),
          m_mesh(makeMesh(dimension)),
          m_space(makeSpace()),
          m_field(m_space),
          m_qf(QF::PolytopeQuadratureFormula::get(2 * Order,
            dimension == 2 ? Polytope::Type::Triangle : Polytope::Type::Tetrahedron))
      {
        const auto fn = [&] {
          if constexpr (!Vector)
            return RealFunction([this](const Point& p) { return analytic(p); });
          else
            return VectorFunction(dimension, [this](const Point& p) {
              Range value(m_dimension);
              for (size_t d = 0; d < m_dimension; ++d)
                value(d) = (d + 1) * analytic(p);
              return value;
            });
        }();
        m_field.project(fn);
        m_points.reserve(m_mesh.getCellCount() * m_qf.getSize());
        for (Index cell = 0; cell < m_mesh.getCellCount(); ++cell)
        {
          const Polytope polytope(dimension, cell, m_mesh);
          for (size_t sample = 0; sample < m_qf.getSize(); ++sample)
          {
            m_points.emplace_back(polytope, m_qf.getPoint(sample));
            // Populate lazy physical coordinates and the numerical oracle once.
            m_expectedSum += reducedExpected(m_points.back());
          }
        }
        m_firstExpected = reducedExpected(m_points.front());
        // Verify every sample and every vector component outside timing. A
        // checksum alone could hide compensating errors on a uniform mesh.
        for (size_t sample = 0; sample < m_points.size(); ++sample)
        {
          const Point& p = m_points[sample];
          const IntegrationPoint ip(p, &m_qf, sample % m_qf.getSize());
          for (const auto& value :
            {m_field.getValue(ip), m_field.getValue(ip), m_field.getValue(p)})
          {
            if constexpr (!Vector)
              m_setupError = std::max(m_setupError, std::abs(value - analytic(p)));
            else
              for (size_t d = 0; d < m_dimension; ++d)
                m_setupError =
                  std::max(m_setupError, std::abs(value(d) - (d + 1) * analytic(p)));
          }
        }
      }

      void run(benchmark::State& state)
      {
        // 0: all misses, 1: 75% hits, 2: 31/32 hits, 3: pure hits.
        const size_t workload = state.range(1);
        const bool pointwise = state.range(2) != 0;
        const size_t repeats = workload == 0 ? 1 : workload == 1 ? 4 : 32;
        const size_t evaluations = m_points.size() * repeats;
        const Real expected =
          workload == 3 ? evaluations * m_firstExpected : repeats * m_expectedSum;
        Range value{};
        // Ensures pure-hit timing starts with its single sample already cached.
        const IntegrationPoint first(m_points.front(), &m_qf, 0);
        m_field.interpolate(value, first);
        Real total = 0;
        for (auto _ : state)
        {
          total = 0;
          if (workload == 3)
          {
            for (size_t evaluation = 0; evaluation < evaluations; ++evaluation)
            {
              if (pointwise)
                m_field.interpolate(value, m_points.front());
              else
                m_field.interpolate(value, first);
              benchmark::DoNotOptimize(value);
              total += reduce(value);
            }
          }
          else
          {
            for (size_t sample = 0; sample < m_points.size(); ++sample)
            {
              const Point& point = m_points[sample];
              const IntegrationPoint ip(point, &m_qf, sample % m_qf.getSize());
              for (size_t repeat = 0; repeat < repeats; ++repeat)
              {
                if (pointwise)
                  m_field.interpolate(value, point);
                else
                  m_field.interpolate(value, ip);
                benchmark::DoNotOptimize(value);
                total += reduce(value);
              }
            }
          }
          benchmark::DoNotOptimize(total);
        }
        const Real relativeError =
          std::abs(total - expected) / std::max(Real(1), std::abs(expected));
        state.counters["relative_error"] = relativeError;
        state.counters["cells"] = m_mesh.getCellCount();
        state.counters["dofs"] = m_space.getSize();
        state.counters["samples"] = m_points.size();
        state.counters["evaluations_per_sweep"] = evaluations;
        state.counters["basis_hit_fraction"] = pointwise ? 0
          : workload == 3                                ? 1
                                                         : Real(repeats - 1) / repeats;
        if (!std::isfinite(relativeError) || relativeError > 1e-9 || m_setupError > 1e-10)
          state.SkipWithError("Large interpolation sweep disagrees with analytic field");
        state.SetItemsProcessed(state.iterations() * evaluations);
      }

    private:
      static LocalMesh makeMesh(size_t dimension)
      {
        LocalMesh mesh = dimension == 2
          ? LocalMesh::UniformGrid(Polytope::Type::Triangle, {129, 129})
          : LocalMesh::UniformGrid(Polytope::Type::Tetrahedron, {17, 17, 17});
        for (size_t d = dimension; d > 0; --d)
          mesh.getConnectivity().compute(d, d - 1);
        return mesh;
      }

      Space makeSpace()
      {
        if constexpr (Order == 1 && !Vector)
          return Space(m_mesh);
        else if constexpr (Order == 1 && Vector)
          return Space(m_mesh, m_dimension);
        else if constexpr (!Vector)
          return Space(std::integral_constant<size_t, Order>{}, m_mesh);
        else
          return Space(std::integral_constant<size_t, Order>{}, m_mesh, m_dimension);
      }

      Real analytic(const Point& p) const
      {
        Real value = 1;
        for (size_t d = 0; d < m_dimension; ++d)
          value += (d + 1) * p.getCoordinates()(d);
        return value;
      }

      Real reducedExpected(const Point& p) const
      {
        if constexpr (Vector)
          return analytic(p) * m_dimension * (m_dimension + 1) / 2;
        else
          return analytic(p);
      }

      static Real reduce(const Range& value)
      {
        if constexpr (Vector)
        {
          Real sum = 0;
          for (size_t d = 0; d < value.size(); ++d)
            sum += value(d);
          return sum;
        }
        else
          return value;
      }

      const size_t m_dimension;
      LocalMesh m_mesh;
      Space m_space;
      Field m_field;
      const QF::QuadratureFormulaBase& m_qf;
      std::vector<Point> m_points;
      Real m_expectedSum = 0;
      Real m_firstExpected = 0;
      Real m_setupError = 0;
  };

  template <size_t Order, bool Vector>
  void largeInterpolation(benchmark::State& state)
  {
    // Persist large fixtures across repetitions; setup is outside timed loops.
    static std::array<std::unique_ptr<LargeInterpolationData<Order, Vector>>, 2> data;
    const size_t index = state.range(0) - 2;
    if (!data[index])
      data[index] =
        std::make_unique<LargeInterpolationData<Order, Vector>>(state.range(0));
    data[index]->run(state);
  }

  const bool registeredLargeInterpolation = [] {
    const char* workloads[] = {"SlowSweep", "Mixed75", "FastBlocks", "PureHit"};
    for (size_t workload = 0; workload < 4; ++workload)
      for (size_t alternative = 0; alternative < 2; ++alternative)
      {
        const std::string suffix = std::string(workloads[workload]) +
          (alternative == 0 ? "/Quadrature" : "/Pointwise");
        const auto args = [&](auto* registration) {
          registration
            ->Args({2, static_cast<int64_t>(workload), static_cast<int64_t>(alternative)})
            ->Args({3, static_cast<int64_t>(workload), static_cast<int64_t>(alternative)})
            ->ArgNames({"dimension", "workload", "pointwise"});
        };
        args(
          benchmark::RegisterBenchmark(("LargeInterpolation/P1Scalar/" + suffix).c_str(),
            &largeInterpolation<1, false>));
        args(
          benchmark::RegisterBenchmark(("LargeInterpolation/P1Vector/" + suffix).c_str(),
            &largeInterpolation<1, true>));
        args(
          benchmark::RegisterBenchmark(("LargeInterpolation/P2Scalar/" + suffix).c_str(),
            &largeInterpolation<2, false>));
        args(
          benchmark::RegisterBenchmark(("LargeInterpolation/P2Vector/" + suffix).c_str(),
            &largeInterpolation<2, true>));
      }
    return true;
  }();
  /// @endcond
}
