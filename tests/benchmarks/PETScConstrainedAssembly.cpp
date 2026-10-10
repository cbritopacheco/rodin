/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief Full constrained PETSc local/MPI assembly, with aligned timing. */
#include <chrono>
#include <iostream>
#include <boost/mpi.hpp>
#include "Rodin/PETSc.h"
#include "Rodin/MPI/Variational/H1/H1.h"
#include "Rodin/MPI/Geometry/Sharder.h"
#include "Rodin/Geometry/BalancedCompactPartitioner.h"
#include "ConstrainedWorkload.h"
#ifdef RODIN_USE_OPENMP
#include <omp.h>
#endif

namespace Rodin::Tests::Benchmarks
{
  bool constrainedCheckFailed = false;

  /**
   * @brief Measures complete PETSc problem assembly with independent patch checks.
   * Fixed iterations align collective calls. The maximum rank duration includes
   * constraint assembly and matrix/vector completion but excludes the preceding
   * barrier, timing reduction, interpolation and numerical checks.
   */
  template <size_t K, bool Distributed>
  class PETScConstrainedAssembly
  {
    public:
      static void run(benchmark::State& state, Geometry::Polytope::Type geometry,
        unsigned physics, const Context::MPI& context)
      {
        using namespace Variational;
        using Workload = ConstrainedWorkload<K, PetscScalar>;
        const auto& comm = context.getCommunicator();
        auto mesh = [&]() {
          if constexpr (Distributed)
          {
            Geometry::Sharder<Context::MPI> sharder(context);
            if (comm.rank() == 0)
            {
              auto root = Convergence::UniformGrid(geometry).makeMesh(state.range(0));
              const size_t dimension = root.getDimension();
              for (size_t d = 0; d <= dimension; ++d)
                for (size_t dp = 0; dp <= dimension; ++dp)
                  root.getConnectivity().compute(d, dp);
              Geometry::BalancedCompactPartitioner partitioner(root);
              partitioner.partition(static_cast<size_t>(comm.size()));
              sharder.shard(partitioner);
            }
            sharder.scatter(0);
            return sharder.gather(0);
          }
          else
            return Convergence::UniformGrid(geometry).makeMesh(state.range(0));
        }();
        const size_t dim = Geometry::Polytope::Traits(geometry).getDimension();
        H1<K, PetscScalar, decltype(mesh)> space(
          std::integral_constant<size_t, K>{}, mesh);
        PETSc::Variational::TrialFunction u(space);
        PETSc::Variational::TestFunction v(space);
        Problem problem(u, v), baseline(u, v);
        Real operatorScale = 1, boundaryScale = 1;
        Workload::configure(u, v, problem, dim, physics, operatorScale, boundaryScale);
        Workload::template configure<true>(
          u, v, baseline, dim, physics, operatorScale, boundaryScale);
        baseline.assemble();
        problem.assemble();
        const auto valid = [&]() {
          auto& system = problem.getLinearSystem();
          auto& reference = baseline.getLinearSystem();
          bool callsSucceeded = true;
          const auto check = [&](PetscErrorCode ierr) {
            assert(ierr == PETSC_SUCCESS);
            callsSucceeded = callsSucceeded && ierr == PETSC_SUCCESS;
          };
          PetscReal rhsNorm = 0, baselineNorm = 0;
          check(VecNorm(system.getVector(), NORM_2, &rhsNorm));
          check(MatNorm(reference.getOperator(), NORM_FROBENIUS, &baselineNorm));
          Vec residual = nullptr, wrong = nullptr;
          check(VecDuplicate(system.getVector(), &residual));
          check(VecDuplicate(u.getSolution().getData(), &wrong));
          check(MatMult(system.getOperator(), u.getSolution().getData(), residual));
          check(VecAXPY(residual, -1, system.getVector()));
          PetscReal patchNorm = 0, sensitivity = 0, loadNorm = 0, matrixNorm = 0;
          check(VecNorm(residual, NORM_2, &patchNorm));
          check(VecCopy(u.getSolution().getData(), wrong));
          check(VecShift(wrong, 1));
          check(MatMult(system.getOperator(), wrong, residual));
          check(VecAXPY(residual, -1, system.getVector()));
          check(VecNorm(residual, NORM_2, &sensitivity));
          check(VecCopy(system.getVector(), residual));
          check(VecAXPY(residual, -1, reference.getVector()));
          check(VecNorm(residual, NORM_2, &loadNorm));
          Mat difference = nullptr;
          check(MatDuplicate(system.getOperator(), MAT_COPY_VALUES, &difference));
          check(
            MatAXPY(difference, -1, reference.getOperator(), DIFFERENT_NONZERO_PATTERN));
          check(MatNorm(difference, NORM_FROBENIUS, &matrixNorm));
          check(MatDestroy(&difference));
          check(VecDestroy(&wrong));
          check(VecDestroy(&residual));
          const Real denominator = std::max(Real(1), Real(rhsNorm));
          bool ok = callsSucceeded && std::isfinite(patchNorm) &&
            patchNorm / denominator < 1e-10 && sensitivity / denominator > 1e-4 &&
            loadNorm / denominator < 1e-11 &&
            matrixNorm < 1e-11 * std::max(Real(1), Real(baselineNorm));
          if (!ok)
            std::cerr << "rank=" << comm.rank() << " patch=" << patchNorm
                      << " sensitivity=" << sensitivity << " load=" << loadNorm
                      << " matrix=" << matrixNorm << " rhs=" << rhsNorm << '\n';
          if constexpr (Distributed)
            ok = boost::mpi::all_reduce(comm, ok, std::logical_and<bool>());
          return ok;
        };
        problem.assemble();
        if (!valid())
        {
          constrainedCheckFailed = true;
          state.SkipWithError("PETSc constrained matrix/load/patch check failed");
          return;
        }
        for (auto _ : state)
        {
          if constexpr (Distributed)
            comm.barrier();
          const auto start = std::chrono::steady_clock::now();
          problem.assemble();
          const auto stop = std::chrono::steady_clock::now();
          double elapsed = std::chrono::duration<double>(stop - start).count();
          if constexpr (Distributed)
            elapsed =
              boost::mpi::all_reduce(comm, elapsed, boost::mpi::maximum<double>());
          state.SetIterationTime(elapsed);
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
          state.SkipWithError(
            "Timed PETSc assembly or coefficient/boundary update failed");
        }
        double owned = 0;
        if constexpr (Distributed)
        {
          for (Index i = 0; i < mesh.getShard().getCellCount(); ++i)
            if (mesh.getShard().isOwned(dim, i))
              ++owned;
        }
        else
          owned = mesh.getPolytopeCount(dim);
        double cells = owned, minimum = owned, maximum = owned;
        if constexpr (Distributed)
        {
          cells = boost::mpi::all_reduce(comm, owned, std::plus<double>());
          minimum = boost::mpi::all_reduce(comm, owned, boost::mpi::minimum<double>());
          maximum = boost::mpi::all_reduce(comm, owned, boost::mpi::maximum<double>());
        }
        MatInfo info;
        const auto ierr =
          MatGetInfo(problem.getLinearSystem().getOperator(), MAT_GLOBAL_SUM, &info);
        assert(ierr == PETSC_SUCCESS);
        if (ierr != PETSC_SUCCESS)
        {
          constrainedCheckFailed = true;
          state.SkipWithError("PETSc matrix metadata query failed");
          return;
        }
        const double ranks = Distributed ? comm.size() : 1;
        state.counters["cells"] = cells;
        state.counters["owned_cells_min"] = minimum;
        state.counters["owned_cells_max"] = maximum;
        state.counters["cell_imbalance"] = ranks * maximum / cells;
        state.counters["ranks"] = ranks;
        state.counters["dofs"] = space.getSize();
        state.counters["nnz"] = info.nz_used;
        state.counters["degree"] = K;
        state.counters["quadrature_order"] = Workload::quadratureOrder;
        state.counters["cells_per_second"] =
          benchmark::Counter(cells, benchmark::Counter::kIsIterationInvariantRate);
      }
  };

  class SilentReporter final : public benchmark::BenchmarkReporter
  {
    public:
      bool ReportContext(const Context&) override
      {
        return true;
      }
      void ReportRuns(const std::vector<Run>&) override {}
  };

  template <size_t K, bool Distributed>
  void registerConstrained(const Context::MPI& context)
  {
    ConstrainedWorkload<K, PetscScalar>::registerCases(
      [&context](auto& state, auto geometry, auto physics) {
        PETScConstrainedAssembly<K, Distributed>::run(state, geometry, physics, context);
      },
      true);
  }
}

int main(int argc, char** argv)
{
  if (PetscInitialize(nullptr, nullptr, nullptr, nullptr) != PETSC_SUCCESS)
    return 1;
  int result = 0;
  {
    boost::mpi::environment environment(argc, argv);
    boost::mpi::communicator world;
    Rodin::Context::MPI context(environment, world);
    bool distributed = false;
    int retained = 1;
    for (int i = 1; i < argc; ++i)
    {
      const std::string_view argument(argv[i]);
      if (argument == "--assembly_mpi")
        distributed = true;
      else if (argument.starts_with("--benchmark_enable_random_interleaving") ||
        argument.starts_with("--benchmark_min_warmup_time"))
        result = 1;
      else if (world.rank() == 0 || !argument.starts_with("--benchmark_out="))
        argv[retained++] = argv[i];
    }
    argc = retained;
    argv[argc] = nullptr;
    if (!distributed && world.size() != 1)
      result = 1;
    if (result == 0)
    {
      using namespace Rodin::Tests::Benchmarks;
      if (distributed)
      {
        registerConstrained<1, true>(context);
        registerConstrained<2, true>(context);
        registerConstrained<3, true>(context);
      }
      else
      {
        registerConstrained<1, false>(context);
        registerConstrained<2, false>(context);
        registerConstrained<3, false>(context);
      }
      benchmark::Initialize(&argc, argv);
      if (benchmark::ReportUnrecognizedArguments(argc, argv))
        result = 1;
      else
      {
        benchmark::AddCustomContext("assembly_backend",
          distributed ? "PETSc/MPI" :
#ifdef RODIN_USE_OPENMP
                      "PETSc/OpenMP"
#else
                      "PETSc/Sequential"
#endif
        );
        benchmark::AddCustomContext("compiler", __VERSION__);
#ifdef RODIN_USE_OPENMP
        benchmark::AddCustomContext(
          "openmp_max_threads", std::to_string(omp_get_max_threads()));
#endif
        SilentReporter silent;
        const size_t cases = world.rank() == 0
          ? benchmark::RunSpecifiedBenchmarks()
          : benchmark::RunSpecifiedBenchmarks(&silent);
        result = cases == 0 || constrainedCheckFailed ? 1 : 0;
      }
      benchmark::Shutdown();
    }
    result = boost::mpi::all_reduce(world, result, boost::mpi::maximum<int>());
  }
  if (PetscFinalize() != PETSC_SUCCESS)
    return 1;
  return result;
}
