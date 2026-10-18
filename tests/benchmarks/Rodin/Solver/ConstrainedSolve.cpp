/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Zero-initial-guess solves paired with constrained assembly workloads.
 */
#include <benchmark/benchmark.h>
#include "Rodin/Assembly.h"
#include "Rodin/Solver/CG.h"
#include "Workloads/ConstrainedWorkload.h"

#ifdef RODIN_USE_OPENMP
#include <omp.h>
#endif

namespace Rodin::Tests::Benchmarks
{
  bool constrainedSolveFailed = false;

  /**
   * @brief Measures CG on the existing manufactured constrained systems.
   *
   * The same @f$A x=b@f$ and exact coefficient vector @f$x_\star@f$ used by
   * the assembly benchmark are retained. The relative residual is
   * @f$\|Ax-b\|_2/\max(1,\|b\|_2)@f$ and the coefficient error is
   * @f$\|x-x_\star\|_2/\max(1,\|x_\star\|_2)@f$.
   * These are algebraic diagnostics, not continuous error norms or rates.
   *
   * ## Architecture
   *
   * Setup, assembly and an original-kernel reference precede timing. A warmup
   * solve supplies the expected iteration count. Every measured solve starts
   * from zero, reset outside timing; otherwise repeated solves would measure
   * an already converged initial guess. The timed public linear-system solve
   * includes CG preconditioner setup and iterations, but not assembly, field
   * writeback, resetting, or validation. Operators, loads, residuals, exact
   * coefficients and iteration counts are checked outside timing. No CTest
   * test, CI timing gate, or additional convergence certification is added.
   * Eigen's reported CG iteration count is retained without reinterpretation;
   * some versions report zero when the first update reaches convergence.
   * Poisson, conductivity and scalar reaction--diffusion are covered here;
   * the coupled reaction--diffusion convergence context is not represented.
   *
   * Complex cases use the same real positive diffusion and reaction
   * coefficients, yielding Hermitian positive-definite systems. They are not
   * indefinite Helmholtz solves. Both real and imaginary manufactured
   * components are checked through the complex residual and coefficient norm.
   *
   * @tparam K Degree of the H1 space.
   * @tparam Scalar Real or complex coefficient field.
   */
  template <size_t K, class Scalar>
  class ConstrainedSolve
  {
    public:
      /// @brief Matches the tolerance in the native Poisson convergence solve.
      static constexpr Real SolverTolerance = 1e-13;
      /// @brief Matches the iteration cap in the native Poisson convergence solve.
      static constexpr size_t MaximumIterations = 20000;
      /// @brief Allows roundoff amplification beyond the requested CG tolerance.
      static constexpr Real ResidualTolerance = 1e-11;
      /// @brief Allows conditioning to amplify algebraic residual into coefficient error.
      static constexpr Real CoefficientTolerance = 1e-9;
      /// @brief Matches the constrained-assembly original-kernel comparison budget.
      static constexpr Real AssemblyTolerance = 1e-11;
      /// @brief Matches the existing wrong-field sensitivity floor.
      static constexpr Real SensitivityFloor = 1e-4;

      /**
       * @brief Profiles solves of a fixed, assembled constrained problem.
       * @param state Grid points per axis and benchmark iterations.
       * @param geometry Reference-cell geometry on the unit box.
       * @param physics Existing workload selector: Poisson, conductivity or reaction.
       */
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
        Problem problem(u, v), reference(u, v);
        const Real operatorScale = 1, boundaryScale = 1;
        Workload::configure(u, v, problem, dim, physics, operatorScale, boundaryScale);
        Workload::template configure<true>(
          u, v, reference, dim, physics, operatorScale, boundaryScale);
        const Math::Vector<Scalar> exact = u.getSolution().getData();
        u.getSolution().getData().setZero();
        reference.assemble();
        problem.assemble();
        auto& system = problem.getLinearSystem();
        const auto& baseline = reference.getLinearSystem();
        Solver::CG solver(problem);
        solver.setTolerance(SolverTolerance).setMaxIterations(MaximumIterations);
        system.getSolution().setZero();
        solver.solve(system);
        const size_t iterations = solver.getIterationNumber();
        const Real rhsScale = std::max(Real(1), baseline.getVector().norm());
        const Real exactScale = std::max(Real(1), exact.norm());
        const auto residual = [&]() {
          return (system.getOperator() * system.getSolution() - system.getVector())
                   .norm() /
            rhsScale;
        };
        const auto coefficientError = [&]() {
          return (system.getSolution() - exact).norm() / exactScale;
        };
        const auto valid = [&]() {
          return solver.success() && solver.getIterationNumber() == iterations &&
            std::isfinite(residual()) && residual() < ResidualTolerance &&
            std::isfinite(coefficientError()) &&
            coefficientError() < CoefficientTolerance &&
            (system.getOperator() - baseline.getOperator()).norm() <
            AssemblyTolerance * std::max(Real(1), baseline.getOperator().norm()) &&
            (system.getVector() - baseline.getVector()).norm() <
            AssemblyTolerance * rhsScale;
        };
        const Math::Vector<Scalar> wrong = exact.array() + Scalar(1);
        const Real sensitivity =
          (system.getOperator() * wrong - system.getVector()).norm() / rhsScale;
        if (!valid() || !(sensitivity > SensitivityFloor))
        {
          constrainedSolveFailed = true;
          state.SkipWithError("Warmup failed manufactured-solution or reference check");
          return;
        }
        bool repeatedSolveValid = true;
        for (auto _ : state)
        {
          state.PauseTiming();
          system.getSolution().setZero();
          state.ResumeTiming();
          solver.solve(system);
          benchmark::DoNotOptimize(system.getSolution().data());
          benchmark::ClobberMemory();
          state.PauseTiming();
          repeatedSolveValid = repeatedSolveValid && solver.success() &&
            solver.getIterationNumber() == iterations;
          state.ResumeTiming();
        }
        if (!repeatedSolveValid || !valid())
        {
          constrainedSolveFailed = true;
          state.SkipWithError("Repeated solve changed solution, iterations or system");
        }
        const size_t cells = mesh.getPolytopeCount(dim);
        size_t localEntries = 0;
        for (size_t cell = 0; cell < cells; ++cell)
        {
          const size_t count = space.getDOFs(dim, cell).size();
          localEntries += count * count;
        }
        state.counters["cells"] = cells;
        state.counters["dofs"] = space.getSize();
        state.counters["system_rows"] = system.getOperator().rows();
        state.counters["nnz"] = system.getOperator().nonZeros();
        state.counters["local_entries_per_operator_term"] = localEntries;
        state.counters["degree"] = K;
        state.counters["quadrature_order"] = Workload::quadratureOrder;
        state.counters["quadrature_points_per_cell"] =
          QF::PolytopeQuadratureFormula::get(Workload::quadratureOrder, geometry)
            .getSize();
        state.counters["solver_reported_iterations"] = iterations;
        state.counters["initial_relative_residual"] =
          system.getVector().norm() / rhsScale;
        state.counters["relative_residual"] = residual();
        state.counters["relative_coefficient_error"] = coefficientError();
      }
  };

  /**
   * @brief Reuses the assembly matrix, names and three grid sizes for solves.
   * @tparam K Degree of the H1 space.
   * @tparam Scalar Real or complex coefficient field.
   */
  template <size_t K, class Scalar>
  void registerConstrainedSolves()
  {
    ConstrainedWorkload<K, Scalar>::registerCases(
      [](auto& state, auto geometry, auto physics) {
        ConstrainedSolve<K, Scalar>::run(state, geometry, physics);
      },
      false);
  }
}

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  registerConstrainedSolves<1, Rodin::Real>();
  registerConstrainedSolves<2, Rodin::Real>();
  registerConstrainedSolves<3, Rodin::Real>();
  registerConstrainedSolves<1, Rodin::Complex>();
  registerConstrainedSolves<2, Rodin::Complex>();
  registerConstrainedSolves<3, Rodin::Complex>();
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
  benchmark::AddCustomContext(
    "timed_operation", "CG setup and solve; zero initial guess");
  benchmark::AddCustomContext("compiler", __VERSION__);
  const size_t cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || constrainedSolveFailed ? 1 : 0;
}
