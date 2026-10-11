/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief Warmed real P1/P2 physics operator assembly on all geometries. */

#include <benchmark/benchmark.h>

#include "Rodin/Assembly.h"
#include "Rodin/Variational/H1.h"
#include "../convergence/Convergence.h"
#include "PhysicsForm.h"

#ifdef RODIN_USE_OPENMP
#include <omp.h>
#endif

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Benchmarks
{
  bool assemblyCheckFailed = false;

  /**
   * @brief Measures repeated assembly of a fixed physical bilinear form.
   *
   * Architecture: mesh, space, and form construction precede timing. An
   * interpolated affine field supplies an independent analytic energy oracle.
   * Two untimed assemblies verify replacement rather than accumulation. The
   * timed loop includes the complete public assembly operation; the final
   * operator is checked outside timing against the initial operator and energy.
   * No solve, constraint elimination, or mesh construction is timed.
   */
  template <size_t Degree, bool Elasticity>
  class PhysicsAssembly
  {
    public:
      static void run(
        benchmark::State& state, Polytope::Type geometry, bool conductivity = false)
      {
        const size_t n = state.range(0);
        auto mesh = Convergence::UniformGrid(geometry).makeMesh(n);
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
        const auto valid = [&]() {
          const auto& matrix = form.getOperator();
          const Real energy = coefficients.dot(matrix * coefficients);
          return std::isfinite(energy) && std::abs(energy / expected - 1) < 1e-9 &&
            (matrix - baseline).norm() < 1e-12 * std::max(Real(1), baseline.norm());
        };
        form.assemble();
        if (!valid())
        {
          assemblyCheckFailed = true;
          state.SkipWithError("Assembly failed analytic-energy or replacement check");
          return;
        }

        for (auto _ : state)
        {
          form.assemble();
          benchmark::DoNotOptimize(form.getOperator().nonZeros());
          benchmark::ClobberMemory();
        }
        if (!valid())
        {
          assemblyCheckFailed = true;
          state.SkipWithError("Timed assembly changed operator or analytic energy");
        }
        const size_t cells = mesh.getPolytopeCount(dim);
        state.counters["cells"] = cells;
        state.counters["dofs"] = space.getSize();
        state.counters["nnz"] = form.getOperator().nonZeros();
        state.counters["degree"] = Degree;
        state.counters["quadrature_order"] = 6;
        state.counters["components"] = Elasticity ? dim : 1;
        state.counters["cells_per_second"] =
          benchmark::Counter(cells, benchmark::Counter::kIsIterationInvariantRate);
      }
  };

  /// @brief Registers three mesh sizes for each physical context and geometry.
  template <size_t Degree>
  void registerPhysics()
  {
    for (auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      const std::string suffix = "/P" + std::to_string(Degree) + "/" +
        std::string(Convergence::UniformGrid::getGeometryName(geometry));
      for (bool conductivity : {false, true})
      {
        const std::string name =
          std::string(conductivity ? "Conductivity" : "Poisson") + suffix;
        benchmark::RegisterBenchmark(name.c_str(),
          [geometry, conductivity](auto& state) {
            PhysicsAssembly<Degree, false>::run(state, geometry, conductivity);
          })
          ->Arg(3)
          ->Arg(5)
          ->Arg(9)
          ->UseRealTime();
      }
      const std::string name = "LinearElasticity" + suffix;
      benchmark::RegisterBenchmark(name.c_str(),
        [geometry](auto& state) { PhysicsAssembly<Degree, true>::run(state, geometry); })
        ->Arg(3)
        ->Arg(5)
        ->Arg(9)
        ->UseRealTime();
    }
  }
}

int main(int argc, char** argv)
{
  Rodin::Tests::Benchmarks::registerPhysics<1>();
  Rodin::Tests::Benchmarks::registerPhysics<2>();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
#ifdef RODIN_USE_OPENMP
  benchmark::AddCustomContext("assembly_backend", "Eigen/OpenMP");
  benchmark::AddCustomContext(
    "openmp_max_threads", std::to_string(omp_get_max_threads()));
#else
  benchmark::AddCustomContext("assembly_backend", "Eigen/Sequential");
#endif
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return Rodin::Tests::Benchmarks::assemblyCheckFailed ? 1 : 0;
}
