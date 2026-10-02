/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Overlapping sequential assembly-stage measurements. */
#include <benchmark/benchmark.h>

#include "Rodin/Assembly.h"
#include "Rodin/QF/PolytopeQuadratureFormula.h"
#include "Rodin/Variational/H1.h"
#include "../convergence/Convergence.h"
#include "PhysicsForm.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Benchmarks
{
  enum class Stage
  {
    Binding,
    Kernel,
    Triplets,
    Finalize
  };
  bool stageCheckFailed = false;

  /**
   * @brief Separates local evaluation from canonical sparse assembly costs.
   * Binding visits all cells and binds each integrator. Kernel additionally
   * evaluates every local entry. Triplets uses the library's sequential
   * assembler, including binding, evaluation, and scatter. Finalize replaces
   * an Eigen sparse matrix using previously generated triplets. These scopes
   * overlap and must not be added or interpreted as disjoint percentages.
   */
  template <size_t Degree, bool Elasticity>
  class PhysicsStages
  {
    public:
      static void run(
        benchmark::State& state, Polytope::Type geometry, Stage stage, bool conductivity)
      {
        auto mesh = Convergence::UniformGrid(geometry).makeMesh(state.range(0));
        const size_t dim = mesh.getDimension();
        using Range = std::conditional_t<Elasticity, Math::SpatialVector<Real>, Real>;
        auto space = [&]() {
          if constexpr (Elasticity)
            return H1<Degree, Range, LocalMesh>(
              std::integral_constant<size_t, Degree>{}, mesh, dim);
          else
            return H1<Degree, Range, LocalMesh>(
              std::integral_constant<size_t, Degree>{}, mesh);
        }();
        TrialFunction u(space);
        TestFunction v(space);
        BilinearForm form(u, v);
        const Real expected =
          PhysicsForm<Elasticity>::configure(u, v, form, dim, conductivity);
        form.assemble();
        BilinearForm original(u, v);
        PhysicsForm<Elasticity>::template configure<true>(
          u, v, original, dim, conductivity);
        original.assemble();
        const auto baseline = original.getOperator();
        const auto& coefficients = u.getSolution().getData();
        using FES = std::decay_t<decltype(space)>;
        using Triplets = std::vector<Math::SparseTriplet<Real>>;
        using TripletForm =
          BilinearForm<typename decltype(u)::SolutionType, FES, FES, Triplets>;
        Assembly::Sequential<Triplets, TripletForm> assembler;
        typename decltype(assembler)::InputType input(
          space, space, form.getLocalIntegrators(), form.getGlobalIntegrators());
        Triplets triplets;
        assembler.execute(triplets, input);
        Math::SparseMatrix<Real> matrix(space.getSize(), space.getSize());
        matrix.setFromTriplets(triplets.begin(), triplets.end());

        const size_t cells = mesh.getPolytopeCount(dim);
        size_t entries = 0;
        for (size_t cell = 0; cell < cells; ++cell)
        {
          const size_t count = space.getDOFs(dim, cell).size();
          entries += count * count * form.getLocalIntegrators().size();
        }
        const auto evaluate = [&](bool bindOnly, bool energy) {
          Real result = 0;
          for (auto& integral : form.getLocalIntegrators())
            for (size_t cell = 0; cell < cells; ++cell)
            {
              const auto polytope = *mesh.getPolytope(dim, cell);
              integral.setPolytope(polytope);
              if (bindOnly)
                continue;
              const auto& dofs = space.getDOFs(dim, cell);
              for (size_t l = 0; l < static_cast<size_t>(dofs.size()); ++l)
                for (size_t m = 0; m < static_cast<size_t>(dofs.size()); ++m)
                {
                  const Real value = integral.integrate(m, l);
                  result += energy ? coefficients(dofs(l)) * value * coefficients(dofs(m))
                                   : value;
                }
            }
          return result;
        };
        const auto valid = [&]() {
          const Real localEnergy = evaluate(false, true);
          const Real matrixEnergy = coefficients.dot(matrix * coefficients);
          return std::isfinite(localEnergy) && std::isfinite(matrixEnergy) &&
            std::abs(localEnergy / expected - 1) < 1e-9 &&
            std::abs(matrixEnergy / expected - 1) < 1e-9 &&
            (matrix - baseline).norm() < 1e-12 * std::max(Real(1), baseline.norm());
        };
        if (!valid())
        {
          stageCheckFailed = true;
          state.SkipWithError("Stage failed local energy or canonical operator check");
          return;
        }
        for (auto _ : state)
        {
          switch (stage)
          {
            case Stage::Binding:
              benchmark::DoNotOptimize(evaluate(true, false));
              break;
            case Stage::Kernel:
              benchmark::DoNotOptimize(evaluate(false, false));
              break;
            case Stage::Triplets:
              assembler.execute(triplets, input);
              benchmark::DoNotOptimize(triplets.data());
              break;
            case Stage::Finalize:
              matrix.setFromTriplets(triplets.begin(), triplets.end());
              benchmark::DoNotOptimize(matrix.nonZeros());
              break;
          }
          benchmark::ClobberMemory();
        }
        // Finalization of the last triplet output is deliberately untimed.
        if (stage == Stage::Triplets)
          matrix.setFromTriplets(triplets.begin(), triplets.end());
        if (!valid())
        {
          stageCheckFailed = true;
          state.SkipWithError("Timed stage changed local energy or operator");
        }
        state.counters["cells"] = cells;
        state.counters["dofs"] = space.getSize();
        state.counters["nnz"] = matrix.nonZeros();
        state.counters["local_entries"] = entries;
        state.counters["triplets"] = triplets.size();
        state.counters["degree"] = Degree;
        state.counters["components"] = Elasticity ? dim : 1;
        state.counters["quadrature_order"] = 6;
        state.counters["quadrature_points_per_cell"] =
          QF::PolytopeQuadratureFormula::get(6, geometry).getSize();
      }
  };

  template <size_t Degree>
  void registerStages()
  {
    for (auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
      for (auto stage : {Stage::Binding, Stage::Kernel, Stage::Triplets, Stage::Finalize})
        for (size_t physics = 0; physics < 3; ++physics)
        {
          const char* stageName = stage == Stage::Binding ? "Binding"
            : stage == Stage::Kernel                      ? "Kernel"
            : stage == Stage::Triplets                    ? "Triplets"
                                                          : "Finalize";
          const std::string name = std::string(stageName) + "/" +
            (physics == 2      ? "LinearElasticity"
                : physics == 1 ? "Conductivity"
                               : "Poisson") +
            "/P" + std::to_string(Degree) + "/" +
            std::string(Convergence::UniformGrid::getGeometryName(geometry));
          benchmark::RegisterBenchmark(name.c_str(),
            [geometry, stage, physics](auto& state) {
              if (physics == 2)
                PhysicsStages<Degree, true>::run(state, geometry, stage, false);
              else
                PhysicsStages<Degree, false>::run(state, geometry, stage, physics == 1);
            })
            ->Arg(3)
            ->Arg(5)
            ->Arg(9)
            ->UseRealTime();
        }
  }
}

int main(int argc, char** argv)
{
  Rodin::Tests::Benchmarks::registerStages<1>();
  Rodin::Tests::Benchmarks::registerStages<2>();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("assembly_backend", "Eigen/Sequential stages");
  benchmark::AddCustomContext("compiler", __VERSION__);
  const auto matched = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return !matched || Rodin::Tests::Benchmarks::stageCheckFailed ? 1 : 0;
}
