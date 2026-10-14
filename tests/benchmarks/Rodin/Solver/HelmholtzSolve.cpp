/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Zero-initial-guess GMRES on manufactured complex Helmholtz operators.
 */
#include <benchmark/benchmark.h>
#include <iostream>
#include "Rodin/Assembly.h"
#include "Rodin/Solver/GMRES.h"
#include "Workloads/HelmholtzWorkload.h"

#ifdef RODIN_USE_OPENMP
#include <omp.h>
#endif

namespace Rodin::Tests::Benchmarks
{
  bool helmholtzSolveFailed = false;

  /**
   * @brief Profiles GMRES without assembly or an already converged initial guess.
   *
   * On the unit box, @f$-\Delta u+\alpha u=f@f$ has full Dirichlet data from
   * @ref Convergence::HelmholtzData. Affine fields are used at degree one and
   * quadratic fields otherwise. The low-wave-number case uses
   * @f$\alpha=-1/4@f$, as in convergence certification. The shifted case uses
   * @f$\alpha=-64+i@f$: its continuous real part is indefinite in dimensions
   * one through three and its imaginary mass term is nonzero.
   * This is not a claim that every coarse discrete interior block has both
   * positive and negative eigenvalues. A smooth zero-boundary probe records
   * and verifies a negative real Rayleigh quotient for the shifted operator.
   *
   * ## Architecture
   *
   * Setup, interpolation, assembly and the original local-kernel comparison
   * are untimed. Each timed public solve includes diagonal preconditioner
   * setup and restarted GMRES iterations. The solution is reset to zero
   * outside timing. Independent true residuals, manufactured coefficients,
   * backend diagnostics and repeated iteration counts qualify the measured
   * output. A wrong field and an omitted mass action must be rejected.
   * The omitted free-free mass block is normalized by its own nonzero action,
   * not by a boundary-dominated load whose scaling changes with refinement.
   * Essential rows and columns are masked through logical boundary DOF indices.
   * There is no CTest registration, timing gate or new convergence claim.
   *
   * @tparam K Degree of the complex H1 space.
   */
  template <size_t K>
  class HelmholtzSolve
  {
    public:
    /// @brief Matches the native scalar convergence solver tolerance.
    static constexpr Real SolverTolerance = 1e-13;
    /// @brief Algebraic validation budget allowing conditioning amplification.
    static constexpr Real ResidualTolerance = 1e-11;
    /// @brief Manufactured coefficient budget, distinct from a continuous norm.
    static constexpr Real CoefficientTolerance = 1e-9;
    /// @brief Original-kernel assembly comparison budget.
    static constexpr Real AssemblyTolerance = 1e-11;
    /// @brief Existing benchmark wrong-field sensitivity floor.
    static constexpr Real SensitivityFloor = 1e-4;
    /// @brief Same generous iteration cap as native scalar convergence solves.
    static constexpr size_t MaximumIterations = 20000;
    /// @brief Fixed restart policy; not tuned separately to each geometry.
    static constexpr size_t Restart = 50;
    /// @brief Integrates polynomial patches and rational pyramid basis products.
    static constexpr size_t QuadratureOrder = HelmholtzWorkload<K>::QuadratureOrder;

    /**
     * @brief Measures a fixed assembled operator and records its numerical work.
     * @param state Grid points per axis and benchmark repetitions.
     * @param geometry Unit-box cell geometry.
     * @param shifted Whether the complex-shifted operator is selected.
     */
    static void run(
      benchmark::State& state, Geometry::Polytope::Type geometry, bool shifted)
    {
      using namespace Variational;
      using Convergence::HelmholtzData;
      auto mesh = Convergence::UniformGrid(geometry).makeMesh(state.range(0));
      const size_t dim = mesh.getDimension();
      const HelmholtzData data(
        dim, K == 1 ? HelmholtzData::Field::Affine : HelmholtzData::Field::Quadratic);
      const Complex alpha = HelmholtzWorkload<K>::getMassCoefficient(shifted);
      H1<K, Complex> space(std::integral_constant<size_t, K>{}, mesh);
      TrialFunction u(space);
      TestFunction v(space);
      Problem problem(u, v), reference(u, v);
      HelmholtzWorkload<K>::configure(u, v, problem, data, alpha);
      HelmholtzWorkload<K>::template configure<true>(u, v, reference, data, alpha);
      const Math::Vector<Complex> exact = u.getSolution().getData();
      u.getSolution().getData().setZero();
      reference.assemble();
      problem.assemble();
      auto& system = problem.getLinearSystem();
      const auto& baseline = reference.getLinearSystem();
      const Real rhsScale = std::max(Real(1), baseline.getVector().norm());
      const Real exactScale = std::max(Real(1), exact.norm());
      BilinearForm mass(u, v);
      auto massIntegral = Integral(u, v);
      massIntegral.setOrder(QuadratureOrder);
      mass = massIntegral;
      mass.assemble();
      auto boundary = DirichletBC(u, data.getSolution());
      boundary.assemble();
      const auto& prescribed = std::get<IndexMap<Complex>>(boundary.getDOFs());
      Math::Vector<Complex> freeExact = exact;
      for (const auto& [index, value] : prescribed)
        freeExact(index) = 0;
      Math::Vector<Complex> missingMassAction = alpha * (mass.getOperator() * freeExact);
      for (const auto& [index, value] : prescribed)
        missingMassAction(index) = 0;
      const Real massActionNorm = missingMassAction.norm();
      if (!std::isfinite(massActionNorm) || !(massActionNorm > 0))
      {
        helmholtzSolveFailed = true;
        state.SkipWithError("Manufactured field has no finite nonzero free mass action");
        return;
      }
      const Real omittedMass = massActionNorm / rhsScale;
      const Real omittedMassResidual =
        (system.getOperator() * exact - system.getVector() - missingMassAction).norm() /
        massActionNorm;
      GridFunction probe(space);
      probe = ComplexFunction([dim](const Geometry::Point& p) {
        Real value = 1;
        for (size_t j = 0; j < dim; ++j)
          value *= std::sin(std::acos(Real(-1)) * p(j));
        return Complex(value);
      });
      const Real rayleigh =
        std::real(probe.getData().dot(system.getOperator() * probe.getData())) /
        probe.getData().squaredNorm();
      Solver::GMRES solver(problem);
      solver.setTolerance(SolverTolerance)
        .setMaxIterations(MaximumIterations)
        .setRestart(Restart);
      system.getSolution().setZero();
      solver.solve(system);
      const size_t iterations = solver.getIterationNumber();
      const auto residual = [&]() {
        return (system.getOperator() * system.getSolution() - system.getVector()).norm() /
          rhsScale;
      };
      const auto coefficientError = [&]() {
        return (system.getSolution() - exact).norm() / exactScale;
      };
      const auto valid = [&]() {
        return solver.success() && solver.getIterationNumber() == iterations &&
          std::isfinite(solver.getError()) && std::isfinite(residual()) &&
          residual() < ResidualTolerance && std::isfinite(coefficientError()) &&
          coefficientError() < CoefficientTolerance &&
          (system.getOperator() - baseline.getOperator()).norm() <
          AssemblyTolerance * std::max(Real(1), baseline.getOperator().norm()) &&
          (system.getVector() - baseline.getVector()).norm() <
          AssemblyTolerance * rhsScale;
      };
      const Math::Vector<Complex> wrong = exact.array() + Complex(1, -1);
      const Real sensitivity =
        (system.getOperator() * wrong - system.getVector()).norm() / rhsScale;
      if (!valid() || !(sensitivity > SensitivityFloor) ||
        !(omittedMassResidual > SensitivityFloor) || !std::isfinite(rayleigh) ||
        (shifted && !(rayleigh < 0)))
      {
        helmholtzSolveFailed = true;
        std::cerr << "Helmholtz qualification: residual=" << residual()
                  << " coefficient_error=" << coefficientError()
                  << " wrong_field=" << sensitivity << " omitted_mass=" << omittedMass
                  << " omitted_mass_normalized_residual=" << omittedMassResidual
                  << " rayleigh=" << rayleigh << " solver_success=" << solver.success()
                  << " matrix_difference="
                  << (system.getOperator() - baseline.getOperator()).norm()
                  << " load_difference="
                  << (system.getVector() - baseline.getVector()).norm() << '\n';
        state.SkipWithError(
          "Warmup failed Helmholtz patch, reference or sensitivity guard");
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
        repeatedValid =
          repeatedValid && solver.success() && solver.getIterationNumber() == iterations;
        state.ResumeTiming();
      }
      if (!repeatedValid || !valid())
      {
        helmholtzSolveFailed = true;
        state.SkipWithError(
          "Repeated GMRES changed diagnostics or manufactured solution");
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
      state.counters["cells"] = cells;
      state.counters["rows"] = system.getOperator().rows();
      state.counters["nnz"] = system.getOperator().nonZeros();
      state.counters["local_basis_per_cell"] = space.getDOFs(dim, 0).size();
      state.counters["scalar_bytes"] = sizeof(Complex);
      state.counters["local_entries_per_operator_term"] = localEntries;
      state.counters["quadrature_order"] = QuadratureOrder;
      state.counters["quadrature_points_per_cell"] =
        QF::PolytopeQuadratureFormula::get(QuadratureOrder, geometry).getSize();
      state.counters["solver_reported_iterations"] = iterations;
      state.counters["solver_reported_error"] = solver.getError();
      state.counters["restart"] = Restart;
      state.counters["mass_real"] = alpha.real();
      state.counters["mass_imaginary"] = alpha.imag();
      state.counters["probe_real_rayleigh"] = rayleigh;
      state.counters["omitted_mass_relative_action"] = omittedMass;
      state.counters["omitted_mass_normalized_residual"] = omittedMassResidual;
      state.counters["relative_residual"] = residual();
      state.counters["relative_coefficient_error"] = coefficientError();
    }

    /// @brief Registers both operators on all seven mesh geometries at three sizes.
    static void registerCases()
    {
      HelmholtzWorkload<K>::registerCases("GMRES",
        [](auto& state, auto geometry, auto shifted) { run(state, geometry, shifted); });
    }
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  HelmholtzSolve<1>::registerCases();
  HelmholtzSolve<2>::registerCases();
  HelmholtzSolve<3>::registerCases();
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
    "timed_operation", "GMRES setup and solve; zero initial guess");
  benchmark::AddCustomContext("compiler", __VERSION__);
  const size_t cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || helmholtzSolveFailed ? 1 : 0;
}
