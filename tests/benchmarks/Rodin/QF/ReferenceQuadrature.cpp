/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Profiles owned rules, traversal and controlled canonical working sets.
 */
#include <benchmark/benchmark.h>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/Geometry/Polytope.h"
#include "unit/Rodin/QF/QuadratureInvariants.h"

namespace Rodin::Tests::Benchmarks
{
  namespace Quadrature = ::Rodin::QF;
  namespace Verification = ::Rodin::Tests::QF;

  bool referenceQuadratureFailed = false;

  enum class ReferenceQuadratureOperation
  {
    Build,
    Copy,
    Traverse,
    HotSweep,
    PoolSweep,
    MoveDiscardedSource
  };

  /**
   * @brief Separates rule lifetimes from traversal and canonical lookup.
   *
   * Owned build/copy scopes include destruction. They measure repeated owned
   * construction after library initialization, not a process-first allocation.
   * MoveDiscardedSource alternates owned wrappers between two optional slots.
   * It includes move construction, moved-from destruction, source-identity
   * reading and slot selection, but no rule construction or sample traversal.
   * The moved-from wrapper is discarded without querying its payload. Exact
   * destination samples and independent moments are checked outside timing.
   * Traversal integrates @f$1+\sum_j(x_j+x_j^2)@f$. Before timing, all moments
   * through total degree two are checked against closed-form integrals, with
   * positive finite weights and coordinates in the closed reference cell. This is a
   * numerical sanity guard, not certification of every claimed order; the QF
   * unit suite supplies the full exactness tests.
   *
   * ## Canonical cache protocol
   *
   * The current implementation has eight thread-local entries. Sweeping eight
   * prewarmed distinct keys visits every cache position once, independent of
   * insertion rotation. Sweeping nine prewarmed keys in round-robin order
   * instead misses that hot cache on every access, but reuses the global pool.
   * These scopes include lookup and size access, not rule construction. A
   * change to the source cache capacity requires updating this explicit
   * experimental protocol. No cache mutation API or source hook is introduced.
   * The global pool is prewarmed with the union of all registered request
   * ranges before any measurement, independently of the CLI selection.
   * Its seeded logical-key count is reported in the benchmark context.
   */
  class ReferenceQuadrature
  {
    public:
      /// @brief Dimensionless absolute/relative budget for low-degree moments.
      static constexpr Real MomentTolerance = 1e-12;
      /// @brief Current source hot-cache capacity, not a public QF API contract.
      static constexpr size_t HotWorkingSet = 8;
      /// @brief One more distinct key forces cyclic round-robin eviction.
      static constexpr size_t PoolWorkingSet = HotWorkingSet + 1;
      /// @brief Disjoint point-rule keys used only to evict prior hot entries.
      static constexpr size_t EvictionOrder = 1024;

      /**
       * @brief Checks independent moments and reference-domain membership.
       * @param formula Rule under test.
       * @param geometry Reference family.
       * @param perturb Whether to perturb the first weight to check sensitivity.
       * @returns Whether the low-degree numerical guard accepts the rule.
       */
      static bool valid(const Quadrature::QuadratureFormulaBase& formula,
        Geometry::Polytope::Type geometry, bool perturb = false)
      {
        const size_t dim = Geometry::Polytope::Traits(geometry).getDimension();
        std::array<std::array<size_t, 3>, 10> powers{};
        size_t count = 1;
        for (size_t j = 0; j < dim; ++j)
        {
          powers[count++][j] = 1;
          powers[count++][j] = 2;
          for (size_t k = j + 1; k < dim; ++k)
          {
            powers[count][j] = powers[count][k] = 1;
            ++count;
          }
        }
        std::array<Real, 10> sums{};
        if (formula.getSize() == 0 ||
          !Verification::allPointsInside(formula, geometry, 0))
          return false;
        for (size_t qp = 0; qp < formula.getSize(); ++qp)
        {
          const auto& x = formula.getPoint(qp);
          if (size_t(x.size()) != dim)
            return false;
          for (size_t j = 0; j < dim; ++j)
          {
            if (!std::isfinite(x(j)))
              return false;
          }
          const Real weight = formula.getWeight(qp) + (perturb && qp == 0 ? Real(1) : 0);
          if (!std::isfinite(weight) || weight <= 0)
            return false;
          for (size_t m = 0; m < count; ++m)
          {
            Real value = weight;
            for (size_t j = 0; j < dim; ++j)
              value *= std::pow(x(j), powers[m][j]);
            sums[m] += value;
          }
        }
        for (size_t m = 0; m < count; ++m)
        {
          const Real expected =
            Verification::exactMoment(geometry, powers[m][0], powers[m][1], powers[m][2]);
          if (std::abs(sums[m] - expected) >
            MomentTolerance * std::max(Real(1), std::abs(expected)))
            return false;
        }
        return true;
      }

      /**
       * @brief Measures one rule or a controlled set of canonical requests.
       * @param state Requested order and timing iterations.
       * @param geometry Reference cell family.
       * @param operation Selected operation scope.
       */
      static void run(benchmark::State& state, Geometry::Polytope::Type geometry,
        ReferenceQuadratureOperation operation)
      {
        const size_t order = state.range(0);
        const size_t dim = Geometry::Polytope::Traits(geometry).getDimension();
        const auto fail = [&](const char* message) {
          referenceQuadratureFailed = true;
          state.SkipWithError(message);
        };
        const auto source = Quadrature::PolytopeQuadratureFormula::build(order, geometry);
        if (!valid(*source, geometry) || valid(*source, geometry, true))
        {
          fail("Reference quadrature failed independent moment or sensitivity guard");
          return;
        }
        state.counters["order"] = order;
        state.counters["dimension"] = dim;
        state.counters["quadrature_points"] = source->getSize();
        if (operation == ReferenceQuadratureOperation::HotSweep ||
          operation == ReferenceQuadratureOperation::PoolSweep)
        {
          const size_t count = operation == ReferenceQuadratureOperation::HotSweep
            ? HotWorkingSet
            : PoolWorkingSet;
          std::array<const Quadrature::QuadratureFormulaBase*, PoolWorkingSet> formulas{};
          size_t points = 0;
          // Point rules are always one point, even for these disjoint order keys.
          // Every preceding cache case exits with only measured keys (orders
          // below 64), so none of these eight filler keys can already be hot.
          for (size_t j = 0; j < HotWorkingSet; ++j)
          {
            benchmark::DoNotOptimize(Quadrature::PolytopeQuadratureFormula::get(
              EvictionOrder + j, Geometry::Polytope::Type::Point));
          }
          for (size_t j = 0; j < count; ++j)
          {
            formulas[j] =
              &Quadrature::PolytopeQuadratureFormula::get(order + j, geometry);
            points += formulas[j]->getSize();
            if (!valid(*formulas[j], geometry))
            {
              fail("Canonical quadrature failed independent moment guard");
              return;
            }
          }
          // Establish the exact cyclic working set after any prior registrations.
          for (size_t j = 0; j < count; ++j)
          {
            benchmark::DoNotOptimize(
              Quadrature::PolytopeQuadratureFormula::get(order + j, geometry));
          }
          for (auto _ : state)
          {
            for (size_t j = 0; j < count; ++j)
            {
              const auto& formula =
                Quadrature::PolytopeQuadratureFormula::get(order + j, geometry);
              benchmark::DoNotOptimize(&formula);
              benchmark::DoNotOptimize(formula.getSize());
            }
          }
          for (size_t j = 0; j < count; ++j)
          {
            if (&Quadrature::PolytopeQuadratureFormula::get(order + j, geometry) !=
              formulas[j])
              fail("Canonical identity changed after hot-cache eviction");
          }
          state.counters["requests_per_sweep"] = count;
          state.counters["points_in_working_set"] = points;
        }
        else if (operation == ReferenceQuadratureOperation::Traverse)
        {
          const auto integrate = [&] {
            Real result = 0;
            for (size_t qp = 0; qp < source->getSize(); ++qp)
            {
              const auto& x = source->getPoint(qp);
              Real value = 1;
              for (size_t j = 0; j < dim; ++j)
                value += x(j) + x(j) * x(j);
              result += source->getWeight(qp) * value;
            }
            return result;
          };
          Real result = 0;
          for (auto _ : state)
          {
            benchmark::DoNotOptimize(source.get());
            benchmark::ClobberMemory();
            result = integrate();
            benchmark::DoNotOptimize(result);
          }
          std::array<size_t, 3> powers{};
          Real expected =
            Verification::exactMoment(geometry, powers[0], powers[1], powers[2]);
          for (size_t j = 0; j < dim; ++j)
          {
            powers[j] = 1;
            expected +=
              Verification::exactMoment(geometry, powers[0], powers[1], powers[2]);
            powers[j] = 2;
            expected +=
              Verification::exactMoment(geometry, powers[0], powers[1], powers[2]);
            powers[j] = 0;
          }
          state.counters["moment_error"] = std::abs(result - expected);
          if (!std::isfinite(result) ||
            std::abs(result - expected) >
              MomentTolerance * std::max(Real(1), std::abs(expected)))
            fail("Timed traversal disagrees with analytic integral");
        }
        else if (operation == ReferenceQuadratureOperation::MoveDiscardedSource)
        {
          const Quadrature::PolytopeQuadratureFormula expected(order, geometry);
          std::array<Optional<Quadrature::PolytopeQuadratureFormula>, 2> rules;
          rules[0].emplace(expected);
          size_t active = 0;
          size_t sourceIdentity = rules[active]->getCacheIdentity();
          const auto move = [&] {
            sourceIdentity = rules[active]->getCacheIdentity();
            const size_t destination = 1 - active;
            rules[destination].emplace(std::move(*rules[active]));
            rules[active].reset();
            active = destination;
          };
          const auto unchanged = [&] {
            if (!rules[active] || rules[1 - active] ||
              rules[active]->getCacheIdentity() == sourceIdentity ||
              rules[active]->getGeometry() != geometry ||
              rules[active]->getSize() != expected.getSize())
              return false;
            for (size_t qp = 0; qp < expected.getSize(); ++qp)
            {
              if (rules[active]->getWeight(qp) != expected.getWeight(qp) ||
                rules[active]->getPoint(qp).size() != expected.getPoint(qp).size())
                return false;
              for (size_t j = 0; j < dim; ++j)
              {
                if (rules[active]->getPoint(qp)(j) != expected.getPoint(qp)(j))
                  return false;
              }
            }
            return valid(*rules[active], geometry);
          };
          move();
          if (!unchanged())
          {
            fail("Rule move changed destination samples, identity or ownership");
            return;
          }
          for (auto _ : state)
          {
            move();
            auto* current = &*rules[active];
            benchmark::DoNotOptimize(current);
            benchmark::ClobberMemory();
          }
          if (!unchanged())
            fail("Timed rule moves changed destination samples or ownership");
          state.counters["moves_per_iteration"] = 1;
          state.counters["source_discarded"] = 1;
        }
        else
        {
          const auto construct = [&] {
            return operation == ReferenceQuadratureOperation::Build
              ? Quadrature::PolytopeQuadratureFormula::build(order, geometry)
              : std::unique_ptr<const Quadrature::QuadratureFormulaBase>(source->copy());
          };
          const auto control = construct();
          if (!valid(*control, geometry) || control->getSize() != source->getSize())
          {
            fail("Owned rule construction or copy failed numerical guard");
            return;
          }
          for (auto _ : state)
          {
            const auto formula = construct();
            benchmark::DoNotOptimize(formula.get());
            benchmark::DoNotOptimize(formula->getSize());
          }
        }
      }

      /**
       * @brief Registers boundary cases and seeds a fixed canonical working set.
       *
       * Orders fifteen through seventeen compare the tetrahedral table limit
       * with adjacent Gaussian product sizes on every family. Pyramid orders
       * nine through eleven and triangle/wedge orders forty-nine through
       * fifty-one expose their own table limits. The requested exactness is
       * retained; these registrations do not change any integration workload.
       * The union of measured key ranges and hot-cache eviction keys is
       * prewarmed, preventing CLI selection from changing pool occupancy.
       *
       * @returns Number of distinct seeded canonical keys.
       */
      static size_t registerCases()
      {
        using G = Geometry::Polytope::Type;
        FlatSet<std::pair<G, size_t>> poolKeys;
        for (const auto& [geometry, name] :
          {std::pair{G::Point, "Point"}, {G::Segment, "Segment"},
            {G::Triangle, "Triangle"}, {G::Quadrilateral, "Quadrilateral"},
            {G::Tetrahedron, "Tetrahedron"}, {G::Pyramid, "Pyramid"},
            {G::Hexahedron, "Hexahedron"}, {G::Wedge, "Wedge"}})
        {
          std::vector<int> orders{2, 8, 15, 16, 17};
          if (geometry == G::Triangle || geometry == G::Wedge)
            orders.insert(orders.end(), {49, 50, 51});
          else if (geometry == G::Pyramid)
            orders.insert(orders.end(), {9, 10, 11});
          else if (geometry == G::Tetrahedron)
            orders.push_back(14);
          for (const int order : orders)
          {
            for (size_t j = 0; j < PoolWorkingSet; ++j)
              poolKeys.emplace(geometry, static_cast<size_t>(order) + j);
          }
          for (const auto& [operation, scope] :
            {std::pair{ReferenceQuadratureOperation::Build, "OwnedBuild"},
              {ReferenceQuadratureOperation::Copy, "OwnedCopy"},
              {ReferenceQuadratureOperation::Traverse, "Traversal"},
              {ReferenceQuadratureOperation::HotSweep, "CanonicalHotSweep8"},
              {ReferenceQuadratureOperation::PoolSweep, "CanonicalPoolSweep9"},
              {ReferenceQuadratureOperation::MoveDiscardedSource,
                "OwnedMoveDiscardedSource"}})
          {
            const std::string label =
              std::string("ReferenceQuadrature/") + scope + "/" + name;
            auto* registration = benchmark::RegisterBenchmark(label.c_str(),
              [geometry, operation](auto& state) { run(state, geometry, operation); });
            for (const int order : orders)
              registration->Arg(order);
            registration->UseRealTime();
          }
        }
        for (size_t j = 0; j < HotWorkingSet; ++j)
          poolKeys.emplace(G::Point, EvictionOrder + j);
        for (const auto& [geometry, order] : poolKeys)
        {
          benchmark::DoNotOptimize(
            Quadrature::PolytopeQuadratureFormula::get(order, geometry));
        }
        return poolKeys.size();
      }
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  const size_t poolKeys = ReferenceQuadrature::registerCases();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext("canonical_pool_seed_keys", std::to_string(poolKeys));
  const size_t count = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return count == 0 || referenceQuadratureFailed ? 1 : 0;
}
