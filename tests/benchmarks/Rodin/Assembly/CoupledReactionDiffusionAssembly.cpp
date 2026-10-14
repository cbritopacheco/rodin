/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Full two-field reaction--diffusion assembly, separately from solving.
 */
#include <cmath>
#include <iostream>
#include "Rodin/Assembly.h"
#include "Workloads/CoupledReactionDiffusionWorkload.h"

#ifdef RODIN_USE_OPENMP
#include <omp.h>
#endif

namespace Rodin::Tests::Benchmarks
{
  bool coupledAssemblyFailed = false;

  /**
   * @brief Measures the public coupled problem assembly and block completion.
   *
   * ## Architecture
   *
   * The shared workload defines two diffusion blocks, two diagonal and two
   * off-diagonal reaction terms, two loads and both Dirichlet traces.
   * Original local operator kernels independently qualify the complete matrix.
   * Loads share their production integration path and are also checked through
   * the analytic manufactured patch, not an independent load algorithm.
   * The manufactured block vector is gathered by the public problem assembly,
   * rather than constructing an assumed component ordering.
   *
   * An omitted-coupling problem retains the original sources and traces. Its
   * residual is normalized by the independently assembled omitted mass action.
   * Only essential output rows are masked: prescribed input coefficients still
   * contribute to the missing physical reaction. Logical boundary DOFs and the
   * problem's test offsets determine these rows. No coordinate tests are used.
   *
   * Repeated full assembly is timed after setup; no solve or numerical guard is
   * timed. The operation includes block offsets, insertion, constraints and
   * initial-guess gathering. It is opt-in and is not a CTest timing gate.
   *
   * @tparam K Degree of both H1 component spaces.
   * @tparam Scalar Real or complex coefficient field.
   * @tparam TestWeighted Whether the second diffusion block uses the equivalent
   * real test-side weighting, exercising generic physical-gradient preparation.
   */
  template <size_t K, class Scalar, bool TestWeighted = false>
  class CoupledReactionDiffusionAssembly
  {
    public:
    /// @brief Existing full-operator and load reference-comparison budget.
    static constexpr Real AssemblyTolerance = 1e-11;
    /// @brief Existing manufactured coefficient-patch residual budget.
    static constexpr Real PatchTolerance = 1e-10;
    /// @brief Wrong-field and omitted-term rejection floor, outside timing.
    static constexpr Real SensitivityFloor = 1e-4;

    /**
     * @brief Measures complete warmed assembly of one coupled discretization.
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
      Workload::template configure<false, TestWeighted>(u, w, v, z, problem, data);
      Workload::template configure<true, TestWeighted>(u, w, v, z, reference, data);
      Workload::template configure<false, TestWeighted>(u, w, v, z, omitted, data, true);
      reference.assemble();
      omitted.assemble();
      problem.assemble();
      const auto& system = problem.getLinearSystem();
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
        coupledAssemblyFailed = true;
        state.SkipWithError(
          "Coupled manufactured field has no finite nonzero reaction action");
        return;
      }
      const auto& wrongSystem = omitted.getLinearSystem();
      const Real omittedResidual =
        (wrongSystem.getOperator() * exact - wrongSystem.getVector()).norm() /
        couplingNorm;
      const Math::Vector<Scalar> wrong = exact.array() + Scalar(1);
      const auto residual = [&]() {
        return (system.getOperator() * exact - system.getVector()).norm() / rhsScale;
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
        return std::isfinite(residual()) && residual() < PatchTolerance &&
          std::isfinite(matrixError()) && matrixError() < AssemblyTolerance &&
          std::isfinite(loadError()) && loadError() < AssemblyTolerance &&
          std::isfinite(omittedResidual) && omittedResidual > SensitivityFloor &&
          sensitivity > SensitivityFloor;
      };
      if (!valid())
      {
        coupledAssemblyFailed = true;
        std::cerr << "Coupled assembly qualification: residual=" << residual()
                  << " matrix_error=" << matrixError() << " load_error=" << loadError()
                  << " omitted_coupling_residual=" << omittedResidual << '\n';
        state.SkipWithError(
          "Initial coupled operator, patch or sensitivity check failed");
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
        coupledAssemblyFailed = true;
        state.SkipWithError("Repeated coupled assembly changed the qualified system");
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
        QF::PolytopeQuadratureFormula::get(Workload::QuadratureOrder, geometry).getSize();
      state.counters["scalar_bytes"] = sizeof(Scalar);
      state.counters["relative_residual"] = residual();
      state.counters["relative_matrix_error"] = matrixError();
      state.counters["relative_load_error"] = loadError();
      state.counters["omitted_coupling_normalized_residual"] = omittedResidual;
    }

    /// @brief Registers the common physical workload for this operation.
    static void registerCases()
    {
      CoupledReactionDiffusionWorkload<K, Scalar>::registerCases(
        TestWeighted ? "GenericGradientAssembly" : "Assembly",
        [](auto& state, auto geometry) { run(state, geometry); });
    }
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  CoupledReactionDiffusionAssembly<1, Rodin::Real>::registerCases();
  CoupledReactionDiffusionAssembly<2, Rodin::Real>::registerCases();
  CoupledReactionDiffusionAssembly<3, Rodin::Real>::registerCases();
  CoupledReactionDiffusionAssembly<1, Rodin::Complex>::registerCases();
  CoupledReactionDiffusionAssembly<2, Rodin::Complex>::registerCases();
  CoupledReactionDiffusionAssembly<3, Rodin::Complex>::registerCases();
  CoupledReactionDiffusionAssembly<1, Rodin::Real, true>::registerCases();
  CoupledReactionDiffusionAssembly<2, Rodin::Real, true>::registerCases();
  CoupledReactionDiffusionAssembly<3, Rodin::Real, true>::registerCases();
  CoupledReactionDiffusionAssembly<1, Rodin::Complex, true>::registerCases();
  CoupledReactionDiffusionAssembly<2, Rodin::Complex, true>::registerCases();
  CoupledReactionDiffusionAssembly<3, Rodin::Complex, true>::registerCases();
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
    "timed_operation", "Warmed full coupled reaction-diffusion assembly");
  benchmark::AddCustomContext("compiler", __VERSION__);
  const size_t cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || coupledAssemblyFailed ? 1 : 0;
}
