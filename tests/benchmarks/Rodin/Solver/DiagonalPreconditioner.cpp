/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Diagonal-preconditioner setup and application on constrained operators.
 */
#include <benchmark/benchmark.h>
#include <Eigen/IterativeLinearSolvers>
#include "Rodin/Assembly.h"
#include "Workloads/ConstrainedWorkload.h"

#ifdef RODIN_USE_OPENMP
#include <omp.h>
#endif

namespace Rodin::Tests::Benchmarks
{
  bool diagonalPreconditionerFailed = false;

  /// @brief Distinct lifetimes of the backend used by Rodin's native CG solver.
  enum class DiagonalOperation
  {
    OwnedCompute,
    RetainedCompute,
    Apply
  };

  /**
   * @brief Separates diagonal setup from iteration-time application.
   *
   * The Eigen backend and lower/upper policy match Solver::CG's current native
   * implementation. This executable uses that backend's preconditioner type;
   * it does not add a source inspection API or a second Rodin solver path.
   * Mesh/space setup, manufactured assembly and reference checks are untimed.
   *
   * ## Lifetimes
   *
   * OwnedCompute includes local construction, compute and destruction.
   * RetainedCompute reuses the already-sized diagonal storage, as in repeated
   * fixed-size CG calls. Apply overwrites a preallocated output and excludes
   * compute. OwnedCompute is not process-first or allocator-cold measurement.
   * Fresh and retained instances are checked before and after the timed loop.
   *
   * ## Numerical guards
   *
   * The independent result is @f$y_j=x_j/A_{jj}@f$. These manufactured SPD
   * systems must have finite, nonzero diagonal entries. Wrong outputs are
   * rejected. Recomputing on @f$1.5A@f$ must produce @f$y/1.5@f$, followed by
   * restoration of the original result. Thus stale retained diagonal data is
   * not accepted. No timing gate or CTest entry is added.
   *
   * @tparam K H1 degree.
   * @tparam Scalar Real or complex coefficient field.
   * @tparam Operation Timed setup/application scope.
   */
  template <size_t K, class Scalar, DiagonalOperation Operation>
  class DiagonalPreconditioner
  {
    public:
      /// @brief Shared algebraic comparison budget used by the sparse profiler.
      static constexpr Real Tolerance = 1e-11;

      /**
       * @brief Measures one preconditioner operation on a fixed FE operator.
       *
       * Diagonal search entries count visited compressed-column entries through
       * the diagonal, not instructions, nonzero arithmetic or memory traffic.
       *
       * @param state Grid points per axis and timed iteration range.
       * @param geometry Reference-cell family.
       * @param physics Manufactured scalar equation selector.
       */
      static void run(
        benchmark::State& state, Geometry::Polytope::Type geometry, unsigned physics)
      {
        using namespace Variational;
        using Workload = ConstrainedWorkload<K, Scalar>;
        using Backend = Eigen::ConjugateGradient<Math::SparseMatrix<Scalar>,
          Eigen::Lower | Eigen::Upper>;
        using Preconditioner = typename Backend::Preconditioner;
        static_assert(
          std::is_same_v<Preconditioner, Eigen::DiagonalPreconditioner<Scalar>>);
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
        problem.assemble();
        reference.assemble();
        const auto& system = problem.getLinearSystem();
        const auto& matrix = system.getOperator();
        const auto& rhs = system.getVector();
        const auto& baseline = reference.getLinearSystem();
        const auto fail = [&](const char* message) {
          diagonalPreconditionerFailed = true;
          state.SkipWithError(message);
        };
        const Real rhsScale = std::max(Real(1), rhs.norm());
        if ((matrix - baseline.getOperator()).norm() >=
            Tolerance * std::max(Real(1), baseline.getOperator().norm()) ||
          (rhs - baseline.getVector()).norm() >= Tolerance * rhsScale ||
          !((matrix * exact - rhs).norm() < Real(1e-10) * rhsScale))
        {
          fail("Manufactured operator/load/patch check failed");
          return;
        }
        Math::Vector<Scalar> x(matrix.rows()), expected(matrix.rows()), y(matrix.rows());
        for (Eigen::Index j = 0; j < x.size(); ++j)
        {
          const Scalar diagonal = matrix.coeff(j, j);
          if (!(std::isfinite(std::abs(diagonal)) && std::abs(diagonal) > Real(0)))
          {
            fail("Manufactured SPD operator has an invalid diagonal");
            return;
          }
          const Real a = Real(1) + Real(j % 17) / 19;
          const Real b = Real(1) / 3 + Real(j % 11) / 13;
          if constexpr (std::is_same_v<Scalar, Complex>)
            x(j) = Scalar(a, b);
          else
            x(j) = a;
          expected(j) = x(j) / diagonal;
        }
        size_t searchEntries = 0;
        for (Eigen::Index outer = 0; outer < matrix.outerSize(); ++outer)
        {
          for (typename Math::SparseMatrix<Scalar>::InnerIterator entry(matrix, outer);
               entry; ++entry)
          {
            ++searchEntries;
            if (entry.index() == outer)
              break;
          }
        }
        const auto valid = [&](const auto& result, const auto& target) {
          return (result - target).norm() < Tolerance * std::max(Real(1), target.norm());
        };
        Preconditioner retained;
        retained.compute(matrix);
        y = retained.solve(x);
        const Math::Vector<Scalar> wrong = expected.array() + Scalar(1);
        const auto freshValid = [&]() {
          Preconditioner fresh;
          fresh.compute(matrix);
          const Math::Vector<Scalar> result = fresh.solve(x);
          return fresh.info() == Eigen::Success && valid(result, expected);
        };
        if (retained.info() != Eigen::Success || !valid(y, expected) ||
          valid(wrong, expected) || !freshValid())
        {
          fail("Independent diagonal or wrong-output guard failed");
          return;
        }
        for (auto _ : state)
        {
          if constexpr (Operation == DiagonalOperation::OwnedCompute)
          {
            Preconditioner owned;
            owned.compute(matrix);
            benchmark::DoNotOptimize(&owned);
            benchmark::ClobberMemory();
          }
          else if constexpr (Operation == DiagonalOperation::RetainedCompute)
          {
            retained.compute(matrix);
            benchmark::DoNotOptimize(&retained);
            benchmark::ClobberMemory();
          }
          else
          {
            benchmark::DoNotOptimize(x.data());
            y.noalias() = retained.solve(x);
            benchmark::DoNotOptimize(y.data());
            benchmark::ClobberMemory();
          }
        }
        if constexpr (Operation != DiagonalOperation::Apply)
          y = retained.solve(x);
        bool ok = retained.info() == Eigen::Success && valid(y, expected) && freshValid();
        const Math::SparseMatrix<Scalar> changed = Real(1.5) * matrix;
        retained.compute(changed);
        y = retained.solve(x);
        const Math::Vector<Scalar> changedExpected = expected / Real(1.5);
        ok = ok && valid(y, changedExpected);
        retained.compute(matrix);
        y = retained.solve(x);
        if (!ok || !valid(y, expected))
          fail("Repeated compute/apply or changed-diagonal check failed");
        state.counters["dimension"] = dim;
        state.counters["degree"] = K;
        state.counters["cells"] = mesh.getPolytopeCount(dim);
        state.counters["rows"] = matrix.rows();
        state.counters["nnz"] = matrix.nonZeros();
        state.counters["diagonal_search_entries"] = searchEntries;
        state.counters["inverse_diagonal_bytes"] = matrix.rows() * sizeof(Scalar);
        state.counters["relative_error"] =
          (y - expected).norm() / std::max(Real(1), expected.norm());
      }
  };

  /// @brief Registers a setup/application scope over the shared FE workload axes.
  template <size_t K, class Scalar, DiagonalOperation Operation>
  void registerDiagonalPreconditioner(const char* operation)
  {
    using G = Geometry::Polytope::Type;
    for (auto geometry : {G::Segment, G::Triangle, G::Quadrilateral, G::Tetrahedron,
           G::Pyramid, G::Hexahedron, G::Wedge})
    {
      for (unsigned physics : {0u, 1u, 2u})
      {
        const std::string name = std::string("DiagonalPreconditioner/") + operation +
          "/" +
          (physics == 0      ? "Poisson"
              : physics == 1 ? "Conductivity"
                             : "ReactionDiffusion") +
          "/H1/P" + std::to_string(K) + "/" +
          (std::is_same_v<Scalar, Complex> ? "Complex/" : "Real/") +
          std::string(Convergence::UniformGrid::getGeometryName(geometry));
        benchmark::RegisterBenchmark(name.c_str(),
          [geometry, physics](auto& state) {
            DiagonalPreconditioner<K, Scalar, Operation>::run(state, geometry, physics);
          })
          ->Arg(3)
          ->Arg(5)
          ->Arg(9)
          ->UseRealTime();
      }
    }
  }

  /// @brief Registers three preconditioner lifetimes for one field and degree.
  template <size_t K, class Scalar>
  void registerDiagonalPreconditioner()
  {
    registerDiagonalPreconditioner<K, Scalar, DiagonalOperation::OwnedCompute>(
      "OwnedCompute");
    registerDiagonalPreconditioner<K, Scalar, DiagonalOperation::RetainedCompute>(
      "RetainedCompute");
    registerDiagonalPreconditioner<K, Scalar, DiagonalOperation::Apply>("Apply");
  }
}

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  registerDiagonalPreconditioner<1, Rodin::Real>();
  registerDiagonalPreconditioner<2, Rodin::Real>();
  registerDiagonalPreconditioner<3, Rodin::Real>();
  registerDiagonalPreconditioner<1, Rodin::Complex>();
  registerDiagonalPreconditioner<2, Rodin::Complex>();
  registerDiagonalPreconditioner<3, Rodin::Complex>();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext(
    "timed_operation", "Eigen CG diagonal setup/apply; assembly excluded");
  benchmark::AddCustomContext("compiler", __VERSION__);
#ifdef RODIN_USE_OPENMP
  benchmark::AddCustomContext("assembly_backend", "Eigen/OpenMP");
  benchmark::AddCustomContext(
    "openmp_max_threads", std::to_string(omp_get_max_threads()));
#else
  benchmark::AddCustomContext("assembly_backend", "Eigen/Sequential");
#endif
  const size_t cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || diagonalPreconditionerFailed ? 1 : 0;
}
