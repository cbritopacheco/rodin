/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Construction and ownership transfer of a coupled variational problem.
 */
#include <cmath>
#include "Rodin/Assembly.h"
#include "Workloads/CoupledReactionDiffusionWorkload.h"

namespace Rodin::Tests::Benchmarks
{
  bool problemLifecycleFailed = false;

  /// @brief Separate construction and ownership-transfer lifetimes.
  enum class ProblemLifetime
  {
    Empty,
    Defined,
    Copy,
    Move,
    CopyAssignment,
    MoveAssignment
  };

  /**
   * @brief Measures public problem lifetimes, excluding mesh/space construction.
   *
   * ## Architecture
   *
   * Degree-two real/complex coupled reaction--diffusion fixtures supply the
   * same formulation as the assembly/solve profilers. Sparse and dense native
   * systems are distinguished. Empty and Defined include construction and
   * destruction, with Defined also cloning the body. Copy includes copying
   * the assembled system and owned body. Move starts from a fresh untimed
   * copy; moved-result destruction is timed, moved-from cleanup is not.
   * These are warmed allocator lifetimes, not process-first allocation.
   * Assignment retains an existing destination: copy assignment may reuse
   * fixed-size system buffers; move assignment replaces the previous value.
   * Fresh staging copies and their moved-from cleanup remain untimed.
   *
   * Numerical qualification assembles a transferred body before and after
   * measurement and compares the operator/load exactly with the fixture.
   * Empty construction is qualified by subsequently assigning that body.
   * A wrong empty body must fail the logical integrator-count guard. No PDE
   * assembly, interpolation or numerical guard is included in lifecycle timing.
   *
   * @tparam Scalar Real or complex coefficient field.
   * @tparam Matrix Sparse or dense operator representation.
   * @tparam Lifetime Public operation measured.
   */
  template <class Scalar, class Matrix, ProblemLifetime Lifetime>
  class ProblemLifecycle
  {
    public:
    /// @brief Existing manufactured coefficient-patch residual budget.
    static constexpr Real PatchTolerance = 1e-10;

    /**
     * @brief Measures one lifetime on fixed spaces and a qualified formulation.
     * @param state Grid points per axis and repetition protocol.
     * @param geometry Reference-cell family of the unit-box mesh.
     */
    static void run(benchmark::State& state, Geometry::Polytope::Type geometry)
    {
      using namespace Variational;
      using Workload = CoupledReactionDiffusionWorkload<2, Scalar>;
      using System = Math::LinearSystem<Matrix, Math::Vector<Scalar>>;
      auto mesh = Convergence::UniformGrid(geometry).makeMesh(state.range(0));
      H1<2, Scalar> space(std::integral_constant<size_t, 2>{}, mesh);
      TrialFunction u(space), w(space);
      TestFunction v(space), z(space);
      using ProblemType =
        Problem<System, decltype(u), decltype(w), decltype(v), decltype(z)>;
      ProblemType source(u, w, v, z);
      const Convergence::ReactionDiffusionData data(
        mesh.getDimension(), Convergence::ReactionDiffusionData::Field::Quadratic);
      Workload::configure(u, w, v, z, source, data);
      source.assemble();
      const auto& system = source.getLinearSystem();
      const Math::Vector<Scalar> exact = system.getSolution();
      const Real residual = (system.getOperator() * exact - system.getVector()).norm() /
        std::max(Real(1), system.getVector().norm());
      const auto bodyMatches = [&](const auto& candidate) {
        return candidate.getBody().getLocalBFIs().size() ==
          source.getBody().getLocalBFIs().size() &&
          candidate.getBody().getLFIs().size() == source.getBody().getLFIs().size();
      };
      const auto makeCandidate = [&]() -> ProblemType {
        if constexpr (Lifetime == ProblemLifetime::Empty ||
          Lifetime == ProblemLifetime::Defined)
        {
          ProblemType candidate(u, w, v, z);
          candidate = source.getBody();
          candidate.assemble();
          return candidate;
        }
        else if constexpr (Lifetime == ProblemLifetime::Copy)
          return ProblemType(source);
        else if constexpr (Lifetime == ProblemLifetime::Move)
        {
          auto staging(source);
          return ProblemType(std::move(staging));
        }
        else
        {
          ProblemType candidate(u, w, v, z);
          if constexpr (Lifetime == ProblemLifetime::CopyAssignment)
            candidate = source;
          else
          {
            auto staging(source);
            candidate = std::move(staging);
          }
          return candidate;
        }
      };
      const auto qualify = [&]() {
        auto candidate = makeCandidate();
        if (!bodyMatches(candidate))
          return false;
        candidate.assemble();
        const auto& result = candidate.getLinearSystem();
        return (result.getOperator() - system.getOperator()).norm() == 0 &&
          (result.getVector() - system.getVector()).norm() == 0 &&
          candidate.getTrialOffsets() == source.getTrialOffsets() &&
          candidate.getTestOffsets() == source.getTestOffsets();
      };
      ProblemType wrong(u, w, v, z);
      if (!std::isfinite(residual) || residual >= PatchTolerance || bodyMatches(wrong) ||
        !qualify())
      {
        problemLifecycleFailed = true;
        state.SkipWithError("Initial lifecycle or manufactured fixture guard failed");
        return;
      }
      Optional<ProblemType> staging;
      Optional<ProblemType> destination;
      if constexpr (Lifetime == ProblemLifetime::CopyAssignment ||
        Lifetime == ProblemLifetime::MoveAssignment)
        destination.emplace(u, w, v, z);
      for (auto _ : state)
      {
        if constexpr (Lifetime == ProblemLifetime::Empty ||
          Lifetime == ProblemLifetime::Defined)
        {
          ProblemType candidate(u, w, v, z);
          if constexpr (Lifetime == ProblemLifetime::Defined)
            candidate = source.getBody();
          benchmark::DoNotOptimize(&candidate);
          benchmark::ClobberMemory();
        }
        else if constexpr (Lifetime == ProblemLifetime::Copy)
        {
          auto candidate(source);
          benchmark::DoNotOptimize(&candidate);
          benchmark::ClobberMemory();
        }
        else if constexpr (Lifetime == ProblemLifetime::CopyAssignment)
        {
          *destination = source;
          benchmark::DoNotOptimize(&*destination);
          benchmark::ClobberMemory();
        }
        else
        {
          state.PauseTiming();
          staging.emplace(source);
          state.ResumeTiming();
          {
            if constexpr (Lifetime == ProblemLifetime::Move)
            {
              auto candidate(std::move(*staging));
              benchmark::DoNotOptimize(&candidate);
              benchmark::ClobberMemory();
            }
            else
            {
              *destination = std::move(*staging);
              benchmark::DoNotOptimize(&*destination);
              benchmark::ClobberMemory();
            }
          }
          state.PauseTiming();
          staging.reset();
          state.ResumeTiming();
        }
      }
      bool destinationMatches = true;
      if constexpr (Lifetime == ProblemLifetime::CopyAssignment ||
        Lifetime == ProblemLifetime::MoveAssignment)
      {
        destinationMatches = bodyMatches(*destination) &&
          (destination->getLinearSystem().getOperator() - system.getOperator()).norm() ==
            0 &&
          (destination->getLinearSystem().getVector() - system.getVector()).norm() == 0;
        if (destinationMatches)
        {
          destination->assemble();
          destinationMatches =
            (destination->getLinearSystem().getOperator() - system.getOperator())
                .norm() == 0 &&
            (destination->getLinearSystem().getVector() - system.getVector()).norm() == 0;
        }
      }
      if (!destinationMatches || !qualify())
      {
        problemLifecycleFailed = true;
        state.SkipWithError("Final lifecycle guard changed the operator or load");
      }
      state.counters["dimension"] = mesh.getDimension();
      state.counters["degree"] = 2;
      state.counters["cells"] = mesh.getPolytopeCount(mesh.getDimension());
      state.counters["rows"] = system.getOperator().rows();
      if constexpr (std::is_base_of_v<Eigen::SparseMatrixBase<Matrix>, Matrix>)
        state.counters["stored_entries"] = system.getOperator().nonZeros();
      else
        state.counters["stored_entries"] = system.getOperator().size();
      state.counters["problem_object_bytes"] = sizeof(ProblemType);
      state.counters["scalar_bytes"] = sizeof(Scalar);
      state.counters["relative_patch_residual"] = residual;
    }

    /// @brief Registers seven mesh families and three independent size levels.
    static void registerCases()
    {
      const char* operation = [] {
        if constexpr (Lifetime == ProblemLifetime::Empty)
          return "Empty";
        else if constexpr (Lifetime == ProblemLifetime::Defined)
          return "Defined";
        else if constexpr (Lifetime == ProblemLifetime::Copy)
          return "Copy";
        else if constexpr (Lifetime == ProblemLifetime::Move)
          return "Move";
        else if constexpr (Lifetime == ProblemLifetime::CopyAssignment)
          return "CopyAssignment";
        else
          return "MoveAssignment";
      }();
      using Geometry::Polytope;
      for (const auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
             Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
             Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
      {
        const std::string name = std::string("ProblemLifecycle/") + operation +
          "/H1/P2/" + (std::is_same_v<Scalar, Complex> ? "Complex/" : "Real/") +
          (std::is_base_of_v<Eigen::SparseMatrixBase<Matrix>, Matrix> ? "Sparse/"
                                                                      : "Dense/") +
          std::string(Convergence::UniformGrid::getGeometryName(geometry));
        benchmark::RegisterBenchmark(
          name.c_str(), [geometry](auto& state) { run(state, geometry); })
          ->Arg(2)
          ->Arg(3)
          ->Arg(4)
          ->UseRealTime();
      }
    }
  };

  /**
   * @brief Registers the separate public lifetimes for one representation.
   * @tparam Scalar Real or complex coefficient field.
   * @tparam Matrix Sparse or dense native operator representation.
   */
  template <class Scalar, class Matrix>
  void registerProblemLifetimes()
  {
    ProblemLifecycle<Scalar, Matrix, ProblemLifetime::Empty>::registerCases();
    ProblemLifecycle<Scalar, Matrix, ProblemLifetime::Defined>::registerCases();
    ProblemLifecycle<Scalar, Matrix, ProblemLifetime::Copy>::registerCases();
    ProblemLifecycle<Scalar, Matrix, ProblemLifetime::Move>::registerCases();
    ProblemLifecycle<Scalar, Matrix, ProblemLifetime::CopyAssignment>::registerCases();
    ProblemLifecycle<Scalar, Matrix, ProblemLifetime::MoveAssignment>::registerCases();
  }
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Tests::Benchmarks;
  registerProblemLifetimes<Real, Math::SparseMatrix<Real>>();
  registerProblemLifetimes<Complex, Math::SparseMatrix<Complex>>();
  registerProblemLifetimes<Real, Math::Matrix<Real>>();
  registerProblemLifetimes<Complex, Math::Matrix<Complex>>();
  benchmark::Initialize(&argc, argv);
  if (benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  benchmark::AddCustomContext("compiler", __VERSION__);
  benchmark::AddCustomContext("timed_operation", "Problem lifetime; warmed allocator");
  const auto cases = benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return cases == 0 || problemLifecycleFailed ? 1 : 0;
}
