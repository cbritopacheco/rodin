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

namespace Rodin::Tests::Benchmarks
{
  bool mappedQuadratureCacheFailed = false;

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
    std::vector<std::unique_ptr<QF::GaussLegendre>> formulas;
    const size_t formulaCount = operation == 2 ? 20 : 1;
    for (size_t i = 0; i < formulaCount; ++i)
      formulas.push_back(std::make_unique<QF::GaussLegendre>(geometry, 2));
    const size_t cellCount = operation == 1 ? mesh.getCellCount() : 1;
    for (size_t cell = 0; cell < cellCount; ++cell)
    {
      for (const auto& formula : formulas)
      {
        const auto& quadrature = mesh.getQuadrature(mesh.getDimension(), cell, *formula);
        if (quadrature.getSize() != formula->getSize() ||
          quadrature.getQuadratureFormulaIdentity() != formula->getCacheIdentity())
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
  }
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return Benchmarks::mappedQuadratureCacheFailed ? 1 : 0;
}
