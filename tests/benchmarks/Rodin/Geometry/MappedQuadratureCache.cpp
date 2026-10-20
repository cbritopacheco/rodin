/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Mesh-owned mapped-quadrature lookup along the assembly call tree.
 */
#include <benchmark/benchmark.h>
#include <memory>
#include <string>
#include <vector>
#include "Rodin/Geometry.h"
#include "Rodin/QF/GaussLegendre.h"
#include "convergence/Convergence.h"
#include "../QF/RepeatedFormula.h"

namespace Rodin::Tests::Benchmarks
{
  bool mappedQuadratureCacheFailed = false;

  /**
   * @brief Measures cold index storage separately from owned mapped quadratures.
   *
   * FirstCell constructs one mapped quadrature per live formula; CellSweep
   * constructs one per cell. Cache scopes include index initialization, dense
   * empty-slot allocation, publication and destruction. Owned scopes construct
   * and destroy quadratures directly, without retaining cell-wide storage.
   * These are complete lifecycles, not isolated allocation latencies. All
   * transformations are warmed outside timing, with fixed live formulas.
   * One-formula scopes retain warm reference tables. Twenty-formula scopes
   * exceed the eight-entry reference store and include its retabulation policy;
   * they must not be interpreted as pure mapped-index allocation measurements.
   *
   * Comparing FirstCell across mesh sizes exposes the cost of empty slots at
   * fixed mapped-point count. Comparing full sweeps instead includes retention
   * versus immediate destruction, and must not be interpreted as the same
   * memory-lifetime contract. Slot counts describe the source representation,
   * not measured resident bytes. No collective or field solve is performed.
   *
   * @param state Grid points per axis and number of live formula versions.
   * @param geometry Reference-cell family.
   * @param cache Whether to use a fresh index rather than owned quadratures.
   * @param sweep Whether to visit every cell rather than the first cell only.
   */
  void mappedQuadratureLifecycle(
    benchmark::State& state, Geometry::Polytope::Type geometry, bool cache, bool sweep)
  {
    auto mesh = Convergence::UniformGrid(geometry).makeMesh(state.range(0));
    const size_t dimension = mesh.getDimension();
    const size_t held = mesh.getCellCount();
    const size_t visited = sweep ? held : 1;
    const size_t formulaCount = state.range(1);
    std::vector<const QF::QuadratureFormula*> formulas;
    const auto& source =
      QF::QuadratureFormula::get<QF::GaussLegendre>(geometry, size_t{2});
    for (size_t i = 0; i < formulaCount; ++i)
      formulas.push_back(&QF::QuadratureFormula::get<RepeatedFormula>(source, i));
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      for (const auto& formula : formulas)
      {
        // Warm every transformation once; duplicate rules need one first-cell
        // oracle each, not twenty identical numerical checks on every cell.
        if (cell->getIndex() != 0 && formula != formulas.front())
          continue;
        const Geometry::PolytopeQuadrature warm(*cell, *formula);
        if (warm.getSize() != formula->getSize() ||
          &warm.getQuadratureFormula() != formula)
        {
          mappedQuadratureCacheFailed = true;
          state.SkipWithError("Cold quadrature oracle has incorrect formula association");
          return;
        }
        for (size_t q = 0; q < warm.getSize(); ++q)
        {
          const Geometry::Point direct(*cell, formula->getPoint(q));
          const auto& actual = warm.getPoint(q).getJacobian();
          const auto& expected = direct.getJacobian();
          if (actual.rows() != expected.rows() || actual.cols() != expected.cols())
          {
            mappedQuadratureCacheFailed = true;
            state.SkipWithError("Mapped Jacobian has incorrect dimensions");
            return;
          }
          for (size_t row = 0; row < actual.rows(); ++row)
          {
            for (size_t col = 0; col < actual.cols(); ++col)
            {
              if (actual(row, col) != expected(row, col))
              {
                mappedQuadratureCacheFailed = true;
                state.SkipWithError("Bulk mapped Jacobian differs from pointwise oracle");
                return;
              }
            }
          }
        }
      }
    }
    for (auto _ : state)
    {
      Geometry::PolytopeQuadratureIndex index;
      if (cache)
        index.initialize(dimension);
      for (size_t cell = 0; cell < visited; ++cell)
      {
        const auto polytope = *mesh.getCell(cell);
        for (const auto& formula : formulas)
        {
          if (cache)
          {
            const auto& quadrature = index.get({dimension, cell}, held, *formula, [&]() {
              return std::make_unique<Geometry::PolytopeQuadrature>(polytope, *formula);
            });
            benchmark::DoNotOptimize(&quadrature);
          }
          else
          {
            const Geometry::PolytopeQuadrature quadrature(polytope, *formula);
            benchmark::DoNotOptimize(&quadrature);
          }
        }
      }
    }
    state.counters["dimension"] = dimension;
    state.counters["held_cells"] = held;
    state.counters["visited_cells"] = visited;
    state.counters["formulas"] = formulaCount;
    state.counters["dense_slots"] = cache ? held * formulaCount : 0;
    state.counters["empty_slots"] = cache ? (held - visited) * formulaCount : 0;
    state.counters["mapped_quadratures"] = visited * formulaCount;
    state.counters["mapped_points"] =
      visited * formulaCount * formulas.front()->getSize();
    state.SetItemsProcessed(state.iterations() * visited * formulaCount);
  }

  /**
   * @brief Measures warm mapped-quadrature lookup without mapping or assembly.
   *
   * Hot repeatedly requests the first cell and one formula. CellSweep visits
   * every cell at fixed formula identity. FormulaRotation requests twenty
   * distinct live formulas on one cell, exceeding the sixteen-entry fast
   * formula lookup. The latter isolates lookup order and full-index traversal,
   * not geometrical complexity: all formulas have the same two-point Gaussian
   * construction parameter. Creation and mapped-coordinate guards are untimed.
   * Outputs are escaped and work is normalized by actual public lookups.
   * No reference-table miss, mapping, solver or collective is timed. This
   * standalone executable is opt-in and is not a CTest or CI timing gate.
   *
   * @param state Benchmark iterations.
   * @param geometry Reference-cell family.
   * @param operation Hot, cell-sweep or formula-rotation scope.
   */
  void mappedQuadratureCache(
    benchmark::State& state, Geometry::Polytope::Type geometry, size_t operation)
  {
    auto mesh = Convergence::UniformGrid(geometry).makeMesh(3);
    std::vector<const QF::QuadratureFormula*> formulas;
    const auto& source =
      QF::QuadratureFormula::get<QF::GaussLegendre>(geometry, size_t{2});
    const size_t formulaCount = operation == 2 ? 20 : 1;
    for (size_t i = 0; i < formulaCount; ++i)
      formulas.push_back(&QF::QuadratureFormula::get<RepeatedFormula>(source, i));
    const size_t cellCount = operation == 1 ? mesh.getCellCount() : 1;
    for (size_t cell = 0; cell < cellCount; ++cell)
    {
      for (const auto& formula : formulas)
      {
        const auto& quadrature = mesh.getQuadrature(mesh.getDimension(), cell, *formula);
        if (quadrature.getSize() != formula->getSize() ||
          &quadrature.getQuadratureFormula() != formula)
        {
          mappedQuadratureCacheFailed = true;
          state.SkipWithError("Mapped quadrature has incorrect formula association");
          return;
        }
        for (size_t q = 0; q < quadrature.getSize(); ++q)
        {
          const auto& actual = quadrature.getPoint(q).getReferenceCoordinates();
          const auto& expected = formula->getPoint(q);
          for (size_t j = 0; j < actual.size(); ++j)
          {
            if (actual(j) != expected(j))
            {
              mappedQuadratureCacheFailed = true;
              state.SkipWithError("Mapped quadrature substituted reference coordinates");
              return;
            }
          }
        }
      }
    }
    for (auto _ : state)
    {
      for (size_t cell = 0; cell < cellCount; ++cell)
      {
        for (const auto& formula : formulas)
        {
          const auto& quadrature =
            mesh.getQuadrature(mesh.getDimension(), cell, *formula);
          benchmark::DoNotOptimize(&quadrature);
        }
      }
    }
    state.counters["dimension"] = mesh.getDimension();
    state.counters["held_cells"] = mesh.getCellCount();
    state.counters["visited_cells"] = cellCount;
    state.counters["formulas"] = formulaCount;
    state.counters["points_per_formula"] = formulas.front()->getSize();
    state.counters["lookups_per_iteration"] = cellCount * formulaCount;
    state.SetItemsProcessed(state.iterations() * cellCount * formulaCount);
  }
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Tests;
  for (const auto geometry : {Geometry::Polytope::Type::Segment,
         Geometry::Polytope::Type::Triangle, Geometry::Polytope::Type::Quadrilateral,
         Geometry::Polytope::Type::Tetrahedron, Geometry::Polytope::Type::Pyramid,
         Geometry::Polytope::Type::Hexahedron, Geometry::Polytope::Type::Wedge})
  {
    for (size_t operation = 0; operation < 3; ++operation)
    {
      const std::string scope = operation == 0 ? "Hot"
        : operation == 1                       ? "CellSweep"
                                               : "FormulaRotation";
      const std::string name = "MappedQuadratureCache/" + scope + "/" +
        std::string(Convergence::UniformGrid::getGeometryName(geometry));
      benchmark::RegisterBenchmark(
        name.c_str(), [geometry, operation](benchmark::State& state) {
          Benchmarks::mappedQuadratureCache(state, geometry, operation);
        });
    }
    for (const bool cache : {false, true})
    {
      for (const bool sweep : {false, true})
      {
        const std::string scope =
          std::string(cache ? "Cache" : "Owned") + (sweep ? "CellSweep" : "FirstCell");
        const std::string name = "MappedQuadratureLifecycle/" + scope + "/" +
          std::string(Convergence::UniformGrid::getGeometryName(geometry));
        auto* registration = benchmark::RegisterBenchmark(
          name.c_str(), [geometry, cache, sweep](benchmark::State& state) {
            Benchmarks::mappedQuadratureLifecycle(state, geometry, cache, sweep);
          });
        for (const size_t n : {3, 9, 17})
        {
          registration->Args({static_cast<int64_t>(n), 1});
          if (!sweep)
            registration->Args({static_cast<int64_t>(n), 20});
        }
      }
    }
  }
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return Benchmarks::mappedQuadratureCacheFailed ? 1 : 0;
}
