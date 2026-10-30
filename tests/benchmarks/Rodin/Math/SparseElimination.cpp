/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Serial sparse elimination by system size and constraint distribution.
 */
#include <benchmark/benchmark.h>
#include <algorithm>
#include <array>
#include <ctime>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include "Rodin/Geometry.h"
#include "Rodin/Variational.h"
#include "Rodin/Variational/H1.h"

#ifdef RODIN_USE_OPENMP
#include <omp.h>
#endif

namespace Rodin::Tests::Benchmarks
{
  bool sparseEliminationFailed = false;

  /**
   * @brief Measures the public elimination operation on assembled systems.
   *
   * ## Architecture
   *
   * Scalar stiffness-plus-mass matrices are assembled on affine unit-box meshes.
   * A deterministic coefficient vector defines an exact algebraic load and
   * nonzero prescribed values. Random and contiguous subsets distinguish the
   * effects of count and placement. Every measured call receives restored
   * original coefficients and load in the same compressed storage allocation.
   * Assembly, constraint construction, restoration, and numerical guards are
   * excluded. Algorithm-owned scratch remains part of the public operation.
   *
   * An independent diagonal-projector construction checks the resulting matrix
   * and load, and the prescribed coefficient vector checks the residual. These
   * FEM fixtures retain their assembled sparsity. Separate direct-system
   * benchmarks exercise asymmetric storage and missing diagonals. No solver,
   * prepared transpose index, or parallel elimination is introduced. Timings use the calling thread's CPU clock, including two clock
   * reads per sample, rather than wall time under unrelated system workloads.
   *
   * @tparam Degree Polynomial degree, with the dedicated P1 space at degree one.
   */
  template <size_t Degree>
  class SparseElimination
  {
    public:
      /// @brief Compressed scalar operator representation.
      using Matrix = Math::SparseMatrix<Real>;

      /// @brief Scalar coefficient and load representation.
      using Vector = Math::Vector<Real>;

      /// @brief Native sparse system whose public operation is measured.
      using System = Math::LinearSystem<Matrix, Vector>;

      /// @brief Relative algebraic comparison tolerance.
      static constexpr Real AlgebraTolerance = 1e-11;

      /// @brief Relative exact-vector residual tolerance.
      static constexpr Real ResidualTolerance = 1e-10;

      /// @brief Reproducible ordering of randomly selected constraints.
      static constexpr unsigned SubsetSeed = 1942;

      /**
       * @brief Reads CPU seconds consumed by the calling thread.
       * @return Monotonic thread CPU time in seconds.
       */
      static double getCPUTime()
      {
        timespec value{};
        if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &value) != 0)
        {
          Alert::Exception() << "Unable to read the thread CPU clock." << Alert::Raise;
        }
        return static_cast<double>(value.tv_sec) +
          static_cast<double>(value.tv_nsec) / 1e9;
      }

      /**
       * @brief Measures one system size, percentage, and placement.
       * @param state Grid points per axis, eliminated percentage, and placement.
       * @param geometry Reference-cell family.
       */
      static void run(benchmark::State& state, Geometry::Polytope::Type geometry)
      {
        using namespace Geometry;
        using namespace Variational;
        const size_t dimension = Polytope::Traits(geometry).getDimension();
        const size_t points = state.range(0);
        Array<size_t> shape(dimension);
        std::fill(shape.begin(), shape.end(), points);
        auto mesh = LocalMesh::UniformGrid(geometry, shape);
        mesh.scale(Real(1) / Real(points - 1));
        for (size_t high = 1; high <= dimension; ++high)
        {
          for (size_t low = 0; low < high; ++low)
            mesh.getConnectivity().compute(high, low);
        }
        auto space = [&]() {
          if constexpr (Degree == 1)
            return P1<Real>(mesh);
          else
            return H1<Degree, Real>(std::integral_constant<size_t, Degree>{}, mesh);
        }();
        TrialFunction u(space);
        TestFunction v(space);
        BilinearForm form(u, v);
        form = Integral(Grad(u), Grad(v)) + Integral(u, v);
        form.assemble();
        Matrix original = form.getOperator();
        original.makeCompressed();
        const size_t count = original.cols();
        const size_t eliminated = std::max(size_t(1), count * state.range(1) / 100);
        Vector exact(count);
        for (size_t j = 0; j < count; ++j)
          exact[j] = Real(1) + Real(j % 19) / 23;
        const Vector load = original * exact;
        IndexVector permutation(count);
        std::iota(permutation.begin(), permutation.end(), Index(0));
        std::mt19937 random(SubsetSeed);
        std::shuffle(permutation.begin(), permutation.end(), random);
        Array<unsigned char> selected(count);
        std::fill(selected.begin(), selected.end(), 0);
        if (state.range(2) == 2)
        {
          mesh.getConnectivity().compute(dimension - 1, dimension);
          auto boundary = DirichletBC(u, RealFunction(Real(1)));
          boundary.assemble();
          for (const auto& [index, value] :
            std::get<typename decltype(boundary)::ValueDOFs>(boundary.getDOFs()))
            selected[index] = 1;
        }
        else
        {
          for (size_t j = 0; j < eliminated; ++j)
          {
            const Index index =
              state.range(2) == 0 ? permutation[j] : (count - eliminated) / 2 + j;
            selected[index] = 1;
          }
        }
        IndexMap<Real> fixed;
        Vector prescribed = Vector::Zero(count);
        Eigen::DiagonalMatrix<Real, Eigen::Dynamic> free(count);
        free.diagonal().setOnes();
        std::vector<Eigen::Triplet<Real>> diagonal;
        diagonal.reserve(eliminated);
        for (size_t j = 0; j < count; ++j)
        {
          if (selected[j])
          {
            fixed.emplace(j, exact[j]);
            prescribed[j] = exact[j];
            free.diagonal()[j] = 0;
            diagonal.emplace_back(j, j, 1);
          }
        }
        Matrix constrained(count, count);
        constrained.setFromTriplets(diagonal.begin(), diagonal.end());
        const Matrix expected = free * original * free + constrained;
        Vector expectedLoad = load - original * prescribed;
        for (const auto& [j, value] : fixed)
          expectedLoad[j] = value;
        System working;
        working.getOperator() = original;
        working.getVector() = load;
        const auto restore = [&]() {
          std::copy_n(
            original.valuePtr(), original.nonZeros(), working.getOperator().valuePtr());
          std::copy_n(load.data(), load.size(), working.getVector().data());
        };
        const auto valid = [&]() {
          return working.getOperator().isCompressed() &&
            working.getOperator().nonZeros() == original.nonZeros() &&
            (working.getOperator() - expected).norm() <
            AlgebraTolerance * std::max(Real(1), expected.norm()) &&
            (working.getVector() - expectedLoad).norm() <
            AlgebraTolerance * std::max(Real(1), expectedLoad.norm()) &&
            (working.getOperator() * exact - working.getVector()).norm() <
            ResidualTolerance * std::max(Real(1), load.norm());
        };
        working.eliminate(fixed);
        if (!valid())
        {
          sparseEliminationFailed = true;
          state.SkipWithError("Sparse elimination failed matrix/load/residual guards");
          return;
        }
        for (auto _ : state)
        {
          state.PauseTiming();
          restore();
          benchmark::ClobberMemory();
          state.ResumeTiming();
          const double start = getCPUTime();
          working.eliminate(fixed);
          const double elapsed = getCPUTime() - start;
          benchmark::DoNotOptimize(working.getOperator().valuePtr());
          benchmark::ClobberMemory();
          state.SetIterationTime(elapsed);
        }
        if (!valid())
        {
          sparseEliminationFailed = true;
          state.SkipWithError("Repeated sparse elimination failed numerical guards");
        }
        state.counters["degree"] = Degree;
        state.counters["dimension"] = dimension;
        state.counters["cells"] = mesh.getCellCount();
        state.counters["local_basis"] =
          space.getFiniteElement(dimension, mesh.getCell()->getIndex()).getCount();
        state.counters["dofs"] = count;
        state.counters["nnz"] = original.nonZeros();
        state.counters["eliminated"] = fixed.size();
        state.counters["eliminated_fraction"] = Real(fixed.size()) / Real(count);
      }
  };

  /**
   * @brief Registers size and constraint axes for one field degree.
   * @tparam Degree P1 or degree-three H1 field.
   */
  template <size_t Degree>
  void registerSparseElimination()
  {
    using Geometry = Rodin::Geometry::Polytope::Type;
    for (auto geometry :
      {Geometry::Segment, Geometry::Triangle, Geometry::Quadrilateral,
        Geometry::Tetrahedron, Geometry::Hexahedron, Geometry::Pyramid, Geometry::Wedge})
    {
      const std::string geometryName = [&]() -> std::string {
        switch (geometry)
        {
          case Geometry::Segment:
            return "Segment";
          case Geometry::Triangle:
            return "Triangle";
          case Geometry::Quadrilateral:
            return "Quadrilateral";
          case Geometry::Tetrahedron:
            return "Tetrahedron";
          case Geometry::Hexahedron:
            return "Hexahedron";
          case Geometry::Pyramid:
            return "Pyramid";
          case Geometry::Wedge:
            return "Wedge";
          default:
            return "Unknown";
        }
      }();
      const std::string name =
        std::string("SparseElimination/P") + std::to_string(Degree) + "/" + geometryName;
      auto* operation = benchmark::RegisterBenchmark(name.c_str(),
        [geometry](auto& state) { SparseElimination<Degree>::run(state, geometry); });
      const std::array<int64_t, 4> points = [&]() -> std::array<int64_t, 4> {
        if (geometry == Geometry::Segment)
        {
          if constexpr (Degree == 1)
            return {33, 1001, 10001, 100001};
          else
            return {11, 335, 3335, 33335};
        }
        if constexpr (Degree == 1)
        {
          return Rodin::Geometry::Polytope::Traits(geometry).getDimension() == 2
            ? std::array<int64_t, 4>{33, 101, 225, 317}
            : std::array<int64_t, 4>{11, 22, 37, 47};
        }
        else
        {
          return Rodin::Geometry::Polytope::Traits(geometry).getDimension() == 2
            ? std::array<int64_t, 4>{11, 34, 75, 106}
            : std::array<int64_t, 4>{4, 8, 13, 16};
        }
      }();
      for (auto pointCount : points)
      {
        for (int64_t percent : {1, 2, 5, 10, 40, 50, 80, 90})
        {
          for (int64_t contiguous : {0, 1})
            operation->Args({pointCount, percent, contiguous});
        }
      }
      for (auto pointCount : points)
        operation->Args({pointCount, 0, 2});
      operation->ArgNames({"points", "percent", "placement"})->UseManualTime();
    }
  }

  /**
   * @brief Measures arbitrary directed CSC storage, including missing diagonals.
   *
   * # Architecture
   *
   * Upper-band columns have no guaranteed stored transpose. Every seventh
   * diagonal and every seventeenth column are omitted in the second pattern.
   * Prescribed coefficients define an algebraically exact solution. Independent
   * projectors check the entire result before and after measurement.
   *
   * # Implementation
   *
   * Each invocation receives an exact-capacity copy of the original CSC matrix
   * outside timing. Scratch and missing-diagonal growth belong to the measured
   * call; assembly, restoration and verification do not. The calling thread's
   * CPU time is recorded, with no worker pool.
   *
   * @param state Size, band width, eliminated percentage and missing-entry flag.
   */
  void directSparseElimination(benchmark::State& state)
  {
    using Matrix = Math::SparseMatrix<Real>;
    using Vector = Math::Vector<Real>;
    const size_t count = state.range(0);
    const size_t width = state.range(1);
    const bool missing = state.range(3);
    std::vector<Math::SparseTriplet<Real>> entries;
    for (size_t column = 0; column < count; ++column)
    {
      if (missing && column % 17 == 0)
        continue;
      for (size_t row = column > width ? column - width : 0; row <= column; ++row)
      {
        if (missing && row == column && column % 7 == 0)
          continue;
        entries.emplace_back(row, column, row == column ? Real(4) : Real(0.25));
      }
    }
    Matrix original(count, count);
    original.setFromTriplets(entries.begin(), entries.end());
    Vector exact(count);
    for (size_t index = 0; index < count; ++index)
      exact[index] = Real(1) + Real(index % 19) / 23;
    const Vector load = original * exact;
    IndexVector indices(count);
    std::iota(indices.begin(), indices.end(), Index(0));
    std::mt19937 random(SparseElimination<1>::SubsetSeed);
    std::shuffle(indices.begin(), indices.end(), random);
    const size_t eliminated = std::max(size_t(1), count * state.range(2) / 100);
    IndexMap<Real> fixed;
    Vector prescribed = Vector::Zero(count);
    Eigen::DiagonalMatrix<Real, Eigen::Dynamic> free(count);
    free.diagonal().setOnes();
    std::vector<Math::SparseTriplet<Real>> diagonal;
    for (size_t position = 0; position < eliminated; ++position)
    {
      const Index index = indices[position];
      fixed.emplace(index, exact[index]);
      prescribed[index] = exact[index];
      free.diagonal()[index] = 0;
      diagonal.emplace_back(index, index, 1);
    }
    Matrix identity(count, count);
    identity.setFromTriplets(diagonal.begin(), diagonal.end());
    const Matrix expected = free * original * free + identity;
    const Vector expectedLoad = free * (load - original * prescribed) + prescribed;
    Math::LinearSystem<Matrix, Vector> system;
    const auto restore = [&]() {
      Matrix fresh = original;
      system.getOperator().swap(fresh);
      system.getVector() = load;
    };
    const auto valid = [&]() {
      return system.getOperator().isCompressed() &&
        (system.getOperator() - expected).norm() <
        SparseElimination<1>::AlgebraTolerance * std::max(Real(1), expected.norm()) &&
        (system.getVector() - expectedLoad).norm() <
        SparseElimination<1>::AlgebraTolerance * std::max(Real(1), expectedLoad.norm()) &&
        (system.getOperator() * exact - system.getVector()).norm() <
        SparseElimination<1>::ResidualTolerance * std::max(Real(1), load.norm());
    };
    restore();
    system.eliminate(fixed);
    if (!valid())
    {
      sparseEliminationFailed = true;
      state.SkipWithError("Directed CSC elimination failed numerical guards");
      return;
    }
    for (auto _ : state)
    {
      state.PauseTiming();
      restore();
      benchmark::ClobberMemory();
      state.ResumeTiming();
      const double start = SparseElimination<1>::getCPUTime();
      system.eliminate(fixed);
      const double elapsed = SparseElimination<1>::getCPUTime() - start;
      benchmark::DoNotOptimize(system.getOperator().valuePtr());
      benchmark::ClobberMemory();
      state.SetIterationTime(elapsed);
    }
    if (!valid())
    {
      sparseEliminationFailed = true;
      state.SkipWithError("Repeated directed CSC elimination failed numerical guards");
    }
    state.counters["dofs"] = count;
    state.counters["nnz"] = original.nonZeros();
    state.counters["eliminated"] = fixed.size();
    state.counters["inserted_diagonals"] =
      system.getOperator().nonZeros() - original.nonZeros();
  }

  /// @brief Registers direct arbitrary CSC systems up to one hundred thousand rows.
  void registerDirectSparseElimination()
  {
    auto* operation =
      benchmark::RegisterBenchmark("SparseElimination/Direct", directSparseElimination);
    for (int64_t count : {1000, 10000, 100000})
    {
      for (int64_t width : {8, 64})
      {
        for (int64_t percent : {1, 2, 5, 10})
        {
          for (int64_t missing : {0, 1})
            operation->Args({count, width, percent, missing});
        }
      }
    }
    operation->ArgNames({"dofs", "width", "percent", "missing"})->UseManualTime();
  }

}

int main(int argc, char** argv)
{
  using namespace Rodin::Tests::Benchmarks;
#ifdef RODIN_USE_OPENMP
  omp_set_dynamic(0);
  omp_set_num_threads(1);
#endif
  Eigen::setNbThreads(1);
  registerSparseElimination<1>();
  registerSparseElimination<3>();
  registerDirectSparseElimination();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext(
    "timed_operation", "public sparse eliminate; thread CPU seconds");
  benchmark::AddCustomContext("restoration", "untimed; resident compressed allocation");
  benchmark::AddCustomContext("threads", "1");
  benchmark::AddCustomContext("compiler", __VERSION__);
  const size_t count = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return count == 0 || sparseEliminationFailed ? 1 : 0;
}
