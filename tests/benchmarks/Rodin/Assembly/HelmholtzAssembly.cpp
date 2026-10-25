/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Complete constrained Helmholtz assembly, separately from solving.
 */
#include "Rodin/Assembly.h"
#include "Workloads/HelmholtzWorkload.h"

#ifdef RODIN_USE_OPENMP
#include <omp.h>
#endif

namespace Rodin::Tests::Benchmarks
{
  bool helmholtzAssemblyFailed = false;

  /**
   * @brief Measures repeated full assembly of the shared Helmholtz formulations.
   *
   * The physical forms, exact coefficients, quadrature and grid sizes match
   * @ref HelmholtzWorkload. Setup, interpolation, the original-kernel reference
   * and numerical checks precede timing. Public @c Problem::assemble includes
   * volume kernels, matrix/load insertion and Dirichlet processing; no solve
   * is performed. Repeated calls measure the warmed fixed-layout operation,
   * not process-first allocation or changing coefficients.
   *
   * The manufactured equation @f$Ax_\star=b@f$ and the full matrix/load
   * comparison are checked before and after timing. A wrong coefficient vector
   * must fail the equation; the complex-shifted case must retain a positive
   * imaginary diagonal sum. These are untimed correctness guards, not CTest
   * registrations or performance gates.
   *
   * @tparam K Degree of the complex H1 space.
   */
  template <size_t K>
  class HelmholtzAssembly
  {
    public:
      /// @brief Same original-kernel comparison budget as constrained assembly.
      static constexpr Real AssemblyTolerance = 1e-11;

      /// @brief Same manufactured patch budget as constrained assembly.
      static constexpr Real PatchTolerance = 1e-10;

      /// @brief Same wrong-field rejection floor as constrained assembly.
      static constexpr Real SensitivityFloor = 1e-4;

      /**
       * @brief Measures full problem assembly on a fixed mesh and space.
       * @param state Grid points per axis and benchmark iterations.
       * @param geometry Unit-box reference-cell family.
       * @param shifted Whether the complex-shifted physical operator is selected.
       */
      static void run(
        benchmark::State& state, Geometry::Polytope::Type geometry, bool shifted)
      {
        using namespace Variational;
        using Workload = HelmholtzWorkload<K>;
        using Convergence::HelmholtzData;
        auto mesh = Convergence::UniformGrid(geometry).makeMesh(state.range(0));
        const size_t dim = mesh.getDimension();
        const HelmholtzData data(
          dim, K == 1 ? HelmholtzData::Field::Affine : HelmholtzData::Field::Quadratic);
        const Complex alpha = Workload::getMassCoefficient(shifted);
        H1<K, Complex> space(std::integral_constant<size_t, K>{}, mesh);
        TrialFunction u(space);
        TestFunction v(space);
        Problem problem(u, v), reference(u, v);
        Workload::configure(u, v, problem, data, alpha);
        Workload::template configure<true>(u, v, reference, data, alpha);
        const Math::Vector<Complex> exact = u.getSolution().getData();
        const Math::Vector<Complex> wrong = exact.array() + Complex(1, -1);
        reference.assemble();
        problem.assemble();
        const auto& system = problem.getLinearSystem();
        const auto& baseline = reference.getLinearSystem();
        const Real rhsScale = std::max(Real(1), baseline.getVector().norm());
        const Real matrixScale = std::max(Real(1), baseline.getOperator().norm());
        const auto residual = [&]() {
          return (system.getOperator() * exact - system.getVector()).norm() / rhsScale;
        };
        const auto matrixError = [&]() {
          return (system.getOperator() - baseline.getOperator()).norm() / matrixScale;
        };
        const auto loadError = [&]() {
          return (system.getVector() - baseline.getVector()).norm() / rhsScale;
        };
        const auto imaginaryDiagonal = [&]() {
          Real sum = 0;
          const size_t rows = static_cast<size_t>(system.getOperator().rows());
          for (size_t j = 0; j < rows; ++j)
            sum += system.getOperator().coeff(j, j).imag();
          return sum;
        };
        const auto valid = [&]() {
          const Real sensitivity =
            (system.getOperator() * wrong - system.getVector()).norm() / rhsScale;
          return std::isfinite(residual()) && residual() < PatchTolerance &&
            std::isfinite(matrixError()) && matrixError() < AssemblyTolerance &&
            std::isfinite(loadError()) && loadError() < AssemblyTolerance &&
            sensitivity > SensitivityFloor &&
            (shifted ? imaginaryDiagonal() > 0 : imaginaryDiagonal() == 0);
        };
        if (!valid())
        {
          helmholtzAssemblyFailed = true;
          state.SkipWithError(
            "Initial Helmholtz matrix, load or manufactured patch failed");
          return;
        }
        for (auto _ : state)
        {
          problem.assemble();
          benchmark::DoNotOptimize(system.getOperator().valuePtr());
          benchmark::DoNotOptimize(system.getVector().data());
          benchmark::ClobberMemory();
        }
        if (!valid())
        {
          helmholtzAssemblyFailed = true;
          state.SkipWithError(
            "Repeated Helmholtz assembly changed the qualified operator");
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
        state.counters["quadrature_order"] = Workload::QuadratureOrder;
        state.counters["quadrature_points_per_cell"] =
          QF::QuadratureFormula::get(Workload::QuadratureOrder, geometry).getSize();
        state.counters["mass_real"] = alpha.real();
        state.counters["mass_imaginary"] = alpha.imag();
        state.counters["relative_residual"] = residual();
        state.counters["relative_matrix_error"] = matrixError();
        state.counters["relative_load_error"] = loadError();
        state.counters["imaginary_diagonal_sum"] = imaginaryDiagonal();
      }

      /// @brief Registers the same physical operator matrix as the solver benchmark.
      static void registerCases()
      {
        HelmholtzWorkload<K>::registerCases(
          "Assembly", [](auto& state, auto geometry, auto shifted) {
            run(state, geometry, shifted);
          });
      }
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  HelmholtzAssembly<1>::registerCases();
  HelmholtzAssembly<2>::registerCases();
  HelmholtzAssembly<3>::registerCases();
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
    "timed_operation", "Warmed full Helmholtz problem assembly");
  benchmark::AddCustomContext("compiler", __VERSION__);
  const size_t cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || helmholtzAssemblyFailed ? 1 : 0;
}
