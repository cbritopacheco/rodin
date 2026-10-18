/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Sparse products and reductions on manufactured finite-element systems.
 */
#include <benchmark/benchmark.h>
#include <complex>
#include <limits>
#include "Rodin/Assembly.h"
#include "Workloads/ConstrainedWorkload.h"

#ifdef RODIN_USE_OPENMP
#include <omp.h>
#endif

namespace Rodin::Tests::Benchmarks
{
  bool sparseAlgebraFailed = false;

  /// @brief Distinct source operations, fixed at compile time in timed loops.
  enum class SparseOperation
  {
    Product,
    TransposeProduct,
    Dot,
    SquaredNorm
  };

  /**
   * @brief Isolates the Math operations used inside constrained CG iterations.
   *
   * Matrices are the actual constrained Poisson, conductivity and scalar
   * reaction--diffusion operators, not synthetic dimension-dependent stencils.
   * Mesh/space construction, original-kernel comparison, matrix assembly,
   * vector allocation and numerical checks precede timing. Products overwrite
   * a preallocated output; reductions return a scalar. These scopes exclude
   * preconditioner setup, solver control and assembly.
   *
   * ## Numerical guards
   *
   * An independent serial traversal accumulates products and Hermitian dot
   * products using long double (whose precision is reported, not assumed).
   * Non-collinear complex vector components
   * exercise conjugation. Manufactured coefficients separately check the
   * assembled equation. Wrong outputs must be rejected. No timing assertion,
   * CTest entry, storage-layout change or relaxed arithmetic is introduced.
   *
   * @tparam K Degree of the H1 field.
   * @tparam Scalar Real or complex coefficient field.
   * @tparam Operation Source operation to measure.
   */
  template <size_t K, class Scalar, SparseOperation Operation>
  class SparseAlgebra
  {
    public:
      static constexpr Real AlgebraTolerance = 1e-11;
      static constexpr Real PatchTolerance = 1e-10;
      using Wide = std::conditional_t<std::is_same_v<Scalar, Complex>,
        std::complex<long double>, long double>;

      /**
       * @brief Measures one operation on a retained assembled matrix and vectors.
       *
       * Compiler barriers retain input reads. Reported matrix storage counts
       * describe compressed values/indices, not resident memory or traffic.
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
          sparseAlgebraFailed = true;
          state.SkipWithError(message);
        };
        const Real rhsScale = std::max(Real(1), rhs.norm());
        if (!matrix.isCompressed() ||
          (matrix - baseline.getOperator()).norm() >=
            AlgebraTolerance * std::max(Real(1), baseline.getOperator().norm()) ||
          (rhs - baseline.getVector()).norm() >= AlgebraTolerance * rhsScale ||
          !((matrix * exact - rhs).norm() < PatchTolerance * rhsScale))
        {
          fail("Assembled operator/load/manufactured patch check failed");
          return;
        }

        Math::Vector<Scalar> x(matrix.cols()), z(matrix.cols()), y(matrix.rows());
        for (Eigen::Index j = 0; j < x.size(); ++j)
        {
          const Real a = Real(1) + Real(j % 17) / 19;
          const Real b = Real(1) / 3 + Real(j % 11) / 13;
          if constexpr (std::is_same_v<Scalar, Complex>)
          {
            x(j) = Scalar(a, b);
            z(j) = Scalar(b, -a / 2);
          }
          else
          {
            x(j) = a;
            z(j) = b;
          }
        }
        Math::Vector<Wide> productReference(matrix.rows());
        productReference.setZero();
        for (Eigen::Index outer = 0; outer < matrix.outerSize(); ++outer)
        {
          for (typename Math::SparseMatrix<Scalar>::InnerIterator entry(matrix, outer);
               entry; ++entry)
          {
            if constexpr (Operation == SparseOperation::TransposeProduct)
              productReference(entry.col()) += Wide(entry.value()) * Wide(x(entry.row()));
            else
              productReference(entry.row()) += Wide(entry.value()) * Wide(x(entry.col()));
          }
        }
        Wide dotReference = 0;
        long double normReference = 0;
        for (Eigen::Index j = 0; j < x.size(); ++j)
        {
          if constexpr (std::is_same_v<Scalar, Complex>)
          {
            dotReference += std::conj(Wide(x(j))) * Wide(z(j));
            normReference += std::norm(Wide(x(j)));
          }
          else
          {
            dotReference += Wide(x(j)) * Wide(z(j));
            normReference += Wide(x(j)) * Wide(x(j));
          }
        }
        const Math::Vector<Scalar> expected = productReference.template cast<Scalar>();
        const Scalar expectedDot = Scalar(dotReference);
        const Real expectedNorm = Real(normReference);
        const auto productValid = [&](const auto& value) {
          return (value - expected).norm() <
            AlgebraTolerance * std::max(Real(1), expected.norm());
        };
        const auto dotValid = [&](Scalar value) {
          return std::abs(value - expectedDot) <
            AlgebraTolerance * std::max(Real(1), std::abs(expectedDot));
        };
        const auto normValid = [&](Real value) {
          return std::abs(value - expectedNorm) <
            AlgebraTolerance * std::max(Real(1), expectedNorm);
        };
        const auto product = [&]() {
          if constexpr (Operation == SparseOperation::TransposeProduct)
            y.noalias() = matrix.transpose() * x;
          else
            y.noalias() = matrix * x;
        };
        product();
        Scalar dot = x.dot(z);
        Real norm = x.squaredNorm();
        const Math::Vector<Scalar> wrong = expected.array() + Scalar(1);
        if (!productValid(y) || !dotValid(dot) || !normValid(norm) ||
          productValid(wrong) || dotValid(Scalar(0)) || normValid(Real(0)))
        {
          fail("Independent product/reduction or wrong-output guard failed");
          return;
        }
        for (auto _ : state)
        {
          benchmark::DoNotOptimize(x.data());
          if constexpr (Operation == SparseOperation::Product ||
            Operation == SparseOperation::TransposeProduct)
          {
            product();
            benchmark::DoNotOptimize(y.data());
          }
          else if constexpr (Operation == SparseOperation::Dot)
          {
            benchmark::DoNotOptimize(z.data());
            dot = x.dot(z);
            benchmark::DoNotOptimize(dot);
          }
          else
          {
            norm = x.squaredNorm();
            benchmark::DoNotOptimize(norm);
          }
          benchmark::ClobberMemory();
        }
        if (!productValid(y) || !dotValid(dot) || !normValid(norm))
          fail("Repeated operation failed the independent scalar reference");
        state.counters["dimension"] = dim;
        state.counters["degree"] = K;
        state.counters["cells"] = mesh.getPolytopeCount(dim);
        state.counters["rows"] = matrix.rows();
        state.counters["nnz"] = matrix.nonZeros();
        state.counters["scalar_bytes"] = sizeof(Scalar);
        state.counters["reference_mantissa_bits"] =
          std::numeric_limits<long double>::digits;
        state.counters["compressed_storage_bytes"] =
          matrix.nonZeros() * (sizeof(Scalar) + sizeof(Math::SparseIndex)) +
          (matrix.outerSize() + 1) * sizeof(Math::SparseIndex);
        state.counters["vector_entries"] = x.size();
        state.counters["product_relative_error"] =
          (y - expected).norm() / std::max(Real(1), expected.norm());
        state.counters["dot_relative_error"] =
          std::abs(dot - expectedDot) / std::max(Real(1), std::abs(expectedDot));
        state.counters["norm_relative_error"] =
          std::abs(norm - expectedNorm) / std::max(Real(1), expectedNorm);
      }
  };

  /// @brief Registers the assembly workload axes under their Math operation.
  template <size_t K, class Scalar, SparseOperation Operation>
  void registerSparseAlgebra(const char* operation)
  {
    using G = Geometry::Polytope::Type;
    for (auto geometry : {G::Segment, G::Triangle, G::Quadrilateral, G::Tetrahedron,
           G::Pyramid, G::Hexahedron, G::Wedge})
    {
      for (unsigned physics : {0u, 1u, 2u})
      {
        const std::string name = std::string("SparseAlgebra/") + operation + "/" +
          (physics == 0      ? "Poisson"
              : physics == 1 ? "Conductivity"
                             : "ReactionDiffusion") +
          "/H1/P" + std::to_string(K) + "/" +
          (std::is_same_v<Scalar, Complex> ? "Complex/" : "Real/") +
          std::string(Convergence::UniformGrid::getGeometryName(geometry));
        benchmark::RegisterBenchmark(name.c_str(),
          [geometry, physics](auto& state) {
            SparseAlgebra<K, Scalar, Operation>::run(state, geometry, physics);
          })
          ->Arg(3)
          ->Arg(5)
          ->Arg(9)
          ->UseRealTime();
      }
    }
  }

  /// @brief Registers four source operations for one scalar field and degree.
  template <size_t K, class Scalar>
  void registerSparseAlgebra()
  {
    registerSparseAlgebra<K, Scalar, SparseOperation::Product>("Product");
    registerSparseAlgebra<K, Scalar, SparseOperation::TransposeProduct>(
      "TransposeProduct");
    registerSparseAlgebra<K, Scalar, SparseOperation::Dot>("Dot");
    registerSparseAlgebra<K, Scalar, SparseOperation::SquaredNorm>("SquaredNorm");
  }
}

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
  registerSparseAlgebra<1, Rodin::Real>();
  registerSparseAlgebra<2, Rodin::Real>();
  registerSparseAlgebra<3, Rodin::Real>();
  registerSparseAlgebra<1, Rodin::Complex>();
  registerSparseAlgebra<2, Rodin::Complex>();
  registerSparseAlgebra<3, Rodin::Complex>();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext(
    "timed_operation", "Math product/reduction; assembly excluded");
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
  return cases == 0 || sparseAlgebraFailed ? 1 : 0;
}
