/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Warmed constrained Eigen problem assembly with manufactured oracles. */
#include "Rodin/Assembly.h"
#include <iostream>
#include "ConstrainedWorkload.h"
#ifdef RODIN_USE_OPENMP
#include <omp.h>
#endif

namespace Rodin::Tests::Benchmarks
{
  bool constrainedCheckFailed = false;

  /** @brief Times complete problem assembly; validates operators, loads and patches. */
  template <size_t K, class Scalar>
  class ConstrainedAssembly
  {
    public:
      static void run(
        benchmark::State& state, Geometry::Polytope::Type geometry, unsigned physics)
      {
        using namespace Variational;
        using Workload = ConstrainedWorkload<K, Scalar>;
        auto mesh = Convergence::UniformGrid(geometry).makeMesh(state.range(0));
        const size_t dim = mesh.getDimension();
        H1<K, Scalar> space(std::integral_constant<size_t, K>{}, mesh);
        TrialFunction u(space);
        TestFunction v(space);
        Problem problem(u, v), baseline(u, v);
        Real operatorScale = 1, boundaryScale = 1;
        Workload::configure(u, v, problem, dim, physics, operatorScale, boundaryScale);
        Workload::template configure<true>(
          u, v, baseline, dim, physics, operatorScale, boundaryScale);
        baseline.assemble();
        problem.assemble();
        const auto valid = [&]() {
          const auto& system = problem.getLinearSystem();
          const auto& reference = baseline.getLinearSystem();
          const auto& matrix = system.getOperator();
          const auto& rhs = system.getVector();
          const auto& exact = u.getSolution().getData();
          const Real denominator = std::max(Real(1), rhs.norm());
          const Real residual = (matrix * exact - rhs).norm() / denominator;
          const Math::Vector<Scalar> wrong = exact.array() + Scalar(1);
          const Real sensitivity = (matrix * wrong - rhs).norm() / denominator;
          if (!std::isfinite(residual) || residual >= 1e-10)
            std::cerr << "patch=" << residual << " sensitivity=" << sensitivity
                      << " matrix_difference="
                      << (matrix - reference.getOperator()).norm()
                      << " load_difference=" << (rhs - reference.getVector()).norm()
                      << '\n';
          return std::isfinite(residual) && residual < 1e-10 && sensitivity > 1e-4 &&
            (matrix - reference.getOperator()).norm() <
            1e-11 * std::max(Real(1), reference.getOperator().norm()) &&
            (rhs - reference.getVector()).norm() < 1e-11 * denominator;
        };
        problem.assemble();
        if (!valid())
        {
          constrainedCheckFailed = true;
          state.SkipWithError(
            "Constrained matrix/load/patch or sensitivity check failed");
          return;
        }
        for (auto _ : state)
        {
          problem.assemble();
          benchmark::DoNotOptimize(problem.getLinearSystem().getOperator().nonZeros());
          benchmark::ClobberMemory();
        }
        bool ok = valid();
        operatorScale = 1.5;
        boundaryScale = 2;
        u.getSolution() = Workload::field(dim, boundaryScale);
        baseline.assemble();
        problem.assemble();
        ok = valid() && ok;
        if (!ok)
        {
          constrainedCheckFailed = true;
          state.SkipWithError("Timed assembly or coefficient/boundary update failed");
        }
        const size_t cells = mesh.getPolytopeCount(dim);
        state.counters["cells"] = cells;
        state.counters["dofs"] = space.getSize();
        state.counters["nnz"] = problem.getLinearSystem().getOperator().nonZeros();
        state.counters["degree"] = K;
        state.counters["quadrature_order"] = Workload::quadratureOrder;
        state.counters["cells_per_second"] =
          benchmark::Counter(cells, benchmark::Counter::kIsIterationInvariantRate);
      }
  };

  template <size_t K, class Scalar>
  void registerConstrained()
  {
    ConstrainedWorkload<K, Scalar>::registerCases(
      [](auto& state, auto geometry, auto physics) {
        ConstrainedAssembly<K, Scalar>::run(state, geometry, physics);
      },
      false);
  }
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Tests::Benchmarks;
  registerConstrained<1, Real>();
  registerConstrained<2, Real>();
  registerConstrained<3, Real>();
  registerConstrained<1, Complex>();
  registerConstrained<2, Complex>();
  registerConstrained<3, Complex>();
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
  const size_t cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || constrainedCheckFailed ? 1 : 0;
}
