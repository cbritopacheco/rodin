/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Fixed-system CG on manufactured two-field reaction--diffusion operators.
 */
#include <cmath>
#include <iostream>
#include "Rodin/Assembly.h"
#include "Rodin/Solver/CG.h"
#include "Workloads/CoupledReactionDiffusionWorkload.h"

#ifdef RODIN_USE_OPENMP
#include <omp.h>
#endif

namespace Rodin::Tests::Benchmarks
{
  bool coupledSolveFailed = false;

  /**
   * @brief Measures public CG on the same block operators as coupled assembly.
   *
   * ## Architecture
   *
   * The manufactured fields, original-kernel reference and omitted-coupling
   * control match the assembly benchmark. Public block assembly gathers exact
   * coefficients and exposes logical test offsets; no scalar-space interleaving
   * is assumed. The omitted reaction is normalized by an independent scalar
   * mass action, with essential output rows masked by boundary DOF indices.
   *
   * Setup, interpolation, all assembly and numerical guards precede timing.
   * Each measured public linear-system solve includes preconditioner setup and
   * CG iterations. The initial guess is reset to zero outside timing. Repeated
   * iteration counts, true residuals, component coefficients and reference
   * matrix/load comparisons qualify the result. Assembly, block-field scatter
   * and error integration are not included. No CTest timing gate is added.
   *
   * Real positive diffusion and the positive-definite real reaction matrix
   * yield SPD or Hermitian positive-definite systems. Complex component
   * multipliers are distinct, so one common real phase is not assumed.
   *
   * @tparam K Degree of both H1 component spaces.
   * @tparam Scalar Real or complex coefficient field.
   */
  template <size_t K, class Scalar>
  class CoupledReactionDiffusionSolve
  {
    public:
      /// @brief Matches the native coupled convergence solver tolerance.
      static constexpr Real SolverTolerance = 1e-13;

      /// @brief Existing generous native scalar/mixed iterative safety cap.
      static constexpr size_t MaximumIterations = 20000;

      /// @brief Allows conditioning amplification of the requested residual budget.
      static constexpr Real CoefficientTolerance = 1e-9;

      /// @brief Existing full-operator and load reference-comparison budget.
      static constexpr Real AssemblyTolerance = 1e-11;

      /// @brief Existing independent true-residual budget for native solves.
      static constexpr Real ResidualTolerance = 1e-11;

      /// @brief Wrong-field and omitted-term rejection floor, outside timing.
      static constexpr Real SensitivityFloor = 1e-4;

      /**
       * @brief Measures zero-initial-guess solves on one fixed coupled discretization.
       * @param state Points per axis and benchmark repetitions.
       * @param geometry Unit-box reference-cell family.
       */
      static void run(benchmark::State& state, Geometry::Polytope::Type geometry)
      {
        using namespace Variational;
        using Workload = CoupledReactionDiffusionWorkload<K, Scalar>;
        using Convergence::ReactionDiffusionData;
        auto mesh = Convergence::UniformGrid(geometry).makeMesh(state.range(0));
        const size_t dim = mesh.getDimension();
        const ReactionDiffusionData data(dim,
          K == 1 ? ReactionDiffusionData::Field::Affine
                 : ReactionDiffusionData::Field::Quadratic);
        H1<K, Scalar> space(std::integral_constant<size_t, K>{}, mesh);
        TrialFunction u(space), w(space);
        TestFunction v(space), z(space);
        Problem problem(u, w, v, z), reference(u, w, v, z), omitted(u, w, v, z);
        Workload::configure(u, w, v, z, problem, data);
        Workload::template configure<true>(u, w, v, z, reference, data);
        Workload::configure(u, w, v, z, omitted, data, true);
        reference.assemble();
        omitted.assemble();
        problem.assemble();
        auto& system = problem.getLinearSystem();
        const auto& baseline = reference.getLinearSystem();
        const Math::Vector<Scalar> exact = system.getSolution();
        const Real rhsScale = std::max(Real(1), baseline.getVector().norm());
        const Real matrixScale = std::max(Real(1), baseline.getOperator().norm());
        BilinearForm mass(u, v);
        auto massIntegral = Integral(u, v);
        massIntegral.setOrder(Workload::QuadratureOrder);
        mass = massIntegral;
        mass.assemble();
        Math::Vector<Scalar> couplingAction = Math::Vector<Scalar>::Zero(exact.size());
        const auto& offsets = problem.getTestOffsets();
        couplingAction.segment(offsets[0], space.getSize()) =
          ReactionDiffusionData::Coupling *
          (mass.getOperator() * w.getSolution().getData());
        couplingAction.segment(offsets[1], space.getSize()) =
          ReactionDiffusionData::Coupling *
          (mass.getOperator() * u.getSolution().getData());
        auto boundaryU = DirichletBC(u, Workload::getSolution(data, 0));
        auto boundaryW = DirichletBC(w, Workload::getSolution(data, 1));
        boundaryU.assemble();
        boundaryW.assemble();
        for (const auto& [index, value] : std::get<IndexMap<Scalar>>(boundaryU.getDOFs()))
          couplingAction(offsets[0] + index) = 0;
        for (const auto& [index, value] : std::get<IndexMap<Scalar>>(boundaryW.getDOFs()))
          couplingAction(offsets[1] + index) = 0;
        const Real couplingNorm = couplingAction.norm();
        if (!std::isfinite(couplingNorm) || !(couplingNorm > 0))
        {
          coupledSolveFailed = true;
          state.SkipWithError(
            "Coupled manufactured field has no finite nonzero reaction action");
          return;
        }
        const auto& wrongSystem = omitted.getLinearSystem();
        const Real omittedResidual =
          (wrongSystem.getOperator() * exact - wrongSystem.getVector()).norm() /
          couplingNorm;
        const Math::Vector<Scalar> wrong = exact.array() + Scalar(1);
        Solver::CG solver(problem);
        solver.setTolerance(SolverTolerance).setMaxIterations(MaximumIterations);
        system.getSolution().setZero();
        solver.solve(system);
        const size_t iterations = solver.getIterationNumber();
        const Real exactScale = std::max(Real(1), exact.norm());
        const auto coefficientError = [&]() {
          return (system.getSolution() - exact).norm() / exactScale;
        };
        const auto residual = [&]() {
          return (system.getOperator() * system.getSolution() - system.getVector())
                   .norm() /
            rhsScale;
        };
        const auto matrixError = [&]() {
          return (system.getOperator() - baseline.getOperator()).norm() / matrixScale;
        };
        const auto loadError = [&]() {
          return (system.getVector() - baseline.getVector()).norm() / rhsScale;
        };
        const auto valid = [&]() {
          const Real sensitivity =
            (system.getOperator() * wrong - system.getVector()).norm() / rhsScale;
          return solver.success() && solver.getIterationNumber() == iterations &&
            std::isfinite(solver.getError()) && std::isfinite(residual()) &&
            residual() < ResidualTolerance && std::isfinite(coefficientError()) &&
            coefficientError() < CoefficientTolerance && std::isfinite(matrixError()) &&
            matrixError() < AssemblyTolerance && std::isfinite(loadError()) &&
            loadError() < AssemblyTolerance && std::isfinite(omittedResidual) &&
            omittedResidual > SensitivityFloor && sensitivity > SensitivityFloor;
        };
        if (!valid())
        {
          coupledSolveFailed = true;
          std::cerr << "Coupled solve qualification: residual=" << residual()
                    << " coefficient_error=" << coefficientError()
                    << " solver_success=" << solver.success()
                    << " matrix_error=" << matrixError() << " load_error=" << loadError()
                    << " omitted_coupling_residual=" << omittedResidual << '\n';
          state.SkipWithError(
            "Initial coupled solve, reference or sensitivity check failed");
          return;
        }
        bool repeatedValid = true;
        for (auto _ : state)
        {
          state.PauseTiming();
          system.getSolution().setZero();
          state.ResumeTiming();
          solver.solve(system);
          benchmark::DoNotOptimize(system.getSolution().data());
          benchmark::ClobberMemory();
          state.PauseTiming();
          repeatedValid = repeatedValid && solver.success() &&
            solver.getIterationNumber() == iterations;
          state.ResumeTiming();
        }
        if (!repeatedValid || !valid())
        {
          coupledSolveFailed = true;
          state.SkipWithError(
            "Repeated coupled solve changed the solution, iterations or system");
        }
        const size_t cells = mesh.getPolytopeCount(dim);
        size_t localEntries = 0;
        for (size_t cell = 0; cell < cells; ++cell)
        {
          const size_t count = space.getDOFs(dim, cell).size();
          localEntries += count * count;
        }
        state.counters["dimension"] = dim;
        state.counters["degree"] = K;
        state.counters["fields"] = 2;
        state.counters["diffusion_terms"] = 2;
        state.counters["reaction_terms"] = 4;
        state.counters["cells"] = cells;
        state.counters["rows"] = system.getOperator().rows();
        state.counters["nnz"] = system.getOperator().nonZeros();
        state.counters["local_basis_per_field"] = space.getDOFs(dim, 0).size();
        state.counters["local_entries_per_operator_term"] = localEntries;
        state.counters["quadrature_order"] = Workload::QuadratureOrder;
        state.counters["quadrature_points_per_cell"] =
          QF::PolytopeQuadratureFormula::get(Workload::QuadratureOrder, geometry)
            .getSize();
        state.counters["scalar_bytes"] = sizeof(Scalar);
        state.counters["solver_reported_iterations"] = iterations;
        state.counters["solver_reported_error"] = solver.getError();
        state.counters["relative_coefficient_error"] = coefficientError();
        state.counters["relative_residual"] = residual();
        state.counters["relative_matrix_error"] = matrixError();
        state.counters["relative_load_error"] = loadError();
        state.counters["omitted_coupling_normalized_residual"] = omittedResidual;
      }

      /// @brief Registers the common physical workload for this operation.
      static void registerCases()
      {
        CoupledReactionDiffusionWorkload<K, Scalar>::registerCases(
          "CG", [](auto& state, auto geometry) { run(state, geometry); });
      }
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  CoupledReactionDiffusionSolve<1, Rodin::Real>::registerCases();
  CoupledReactionDiffusionSolve<2, Rodin::Real>::registerCases();
  CoupledReactionDiffusionSolve<3, Rodin::Real>::registerCases();
  CoupledReactionDiffusionSolve<1, Rodin::Complex>::registerCases();
  CoupledReactionDiffusionSolve<2, Rodin::Complex>::registerCases();
  CoupledReactionDiffusionSolve<3, Rodin::Complex>::registerCases();
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
    "timed_operation", "Coupled CG setup and solve; zero initial guess");
  benchmark::AddCustomContext("compiler", __VERSION__);
  const size_t cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || coupledSolveFailed ? 1 : 0;
}
