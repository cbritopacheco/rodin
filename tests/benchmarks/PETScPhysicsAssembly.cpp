/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/** @file @brief PETSc local/MPI warmed physics assembly with aligned iterations. */
#include <chrono>
#include <iostream>
#include <benchmark/benchmark.h>
#include <boost/mpi.hpp>
#include "Rodin/PETSc.h"
#include "Rodin/MPI/Variational/H1/H1.h"
#include "Rodin/MPI/Geometry/Sharder.h"
#include "Rodin/Geometry/BalancedCompactPartitioner.h"
#include "../convergence/Convergence.h"
#include "PhysicsForm.h"

#ifdef RODIN_USE_OPENMP
#include <omp.h>
#endif

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Benchmarks
{
  bool petscCheckFailed = false;

  /**
   * @brief Complete warmed assembly measured independently of setup and oracles.
   * Architecture: identical forms are configured through PhysicsForm. MPI meshes
   * use explicit incidence and root partitioning; only owned cells contribute
   * to workload counts. Fixed iteration counts align collectives. Assembly is
   * timed between steady-clock readings; MPI reports maximum rank duration.
   * Synchronization, timing reduction, and numerical checks are outside timing.
   */
  template <size_t Degree, bool Elasticity, bool Distributed>
  class PETScPhysicsAssembly
  {
    public:
      static void run(benchmark::State& state, Polytope::Type geometry, bool conductivity,
        const Context::MPI& context)
      {
        const auto& comm = context.getCommunicator();
        auto mesh = [&]() {
          if constexpr (!Distributed)
            return Convergence::UniformGrid(geometry).makeMesh(state.range(0));
          else
          {
            Sharder<Context::MPI> sharder(context);
            if (comm.rank() == 0)
            {
              auto root = Convergence::UniformGrid(geometry).makeMesh(state.range(0));
              const size_t dim = root.getDimension();
              root.getConnectivity().compute(dim, dim);
              root.getConnectivity().compute(dim, 0);
              root.getConnectivity().compute(dim, dim - 1);
              root.getConnectivity().compute(dim - 1, dim);
              root.getConnectivity().compute(dim - 1, 0);
              BalancedCompactPartitioner partitioner(root);
              partitioner.partition(static_cast<size_t>(comm.size()));
              sharder.shard(partitioner);
              sharder.scatter(0);
            }
            return sharder.gather(0);
          }
        }();
        // The declared domain dimension is rank-independent; an empty shard
        // need not carry a positive-dimensional local mesh.
        const size_t dim = Polytope::Traits(geometry).getDimension();
        using Range = std::conditional_t<Elasticity, Math::SpatialVector<Real>, Real>;
        auto space = [&]() {
          if constexpr (Elasticity)
            return H1<Degree, Range, decltype(mesh)>(
              std::integral_constant<size_t, Degree>{}, mesh, dim);
          else
            return H1<Degree, Range, decltype(mesh)>(
              std::integral_constant<size_t, Degree>{}, mesh);
        }();
        PETSc::Variational::TrialFunction u(space);
        PETSc::Variational::TestFunction v(space);
        BilinearForm form(u, v);
        const Real expected =
          PhysicsForm<Elasticity>::configure(u, v, form, dim, conductivity);
        form.assemble();
        const auto baseline = form;
        PetscReal baselineNorm = 0;
        auto ierr = MatNorm(baseline.getOperator(), NORM_FROBENIUS, &baselineNorm);
        assert(ierr == PETSC_SUCCESS);
        bool callsSucceeded = ierr == PETSC_SUCCESS;
        const auto valid = [&]() {
          const auto energy = form(u.getSolution(), u.getSolution());
          auto difference = form;
          auto error = MatAXPY(
            difference.getOperator(), -1, baseline.getOperator(), SAME_NONZERO_PATTERN);
          assert(error == PETSC_SUCCESS);
          bool ok = error == PETSC_SUCCESS;
          PetscReal norm = 0;
          error = MatNorm(difference.getOperator(), NORM_FROBENIUS, &norm);
          assert(error == PETSC_SUCCESS);
          ok = ok && error == PETSC_SUCCESS && callsSucceeded &&
            std::isfinite(std::abs(energy)) && std::abs(energy / expected - 1) < 1e-9 &&
            norm < 1e-12 * std::max(Real(1), baselineNorm);
          if (!ok)
            std::cerr << "rank=" << comm.rank() << " energy=" << energy
                      << " expected=" << expected << " difference=" << norm
                      << " baseline_norm=" << baselineNorm << '\n';
          if constexpr (Distributed)
            ok = boost::mpi::all_reduce(comm, ok, std::logical_and<bool>());
          return ok;
        };
        form.assemble();
        if (!valid())
        {
          petscCheckFailed = true;
          state.SkipWithError("PETSc energy or replacement check failed");
          return;
        }
        for (auto _ : state)
        {
          if constexpr (Distributed)
            comm.barrier();
          const auto start = std::chrono::steady_clock::now();
          form.assemble();
          const auto stop = std::chrono::steady_clock::now();
          double elapsed = std::chrono::duration<double>(stop - start).count();
          if constexpr (Distributed)
            elapsed =
              boost::mpi::all_reduce(comm, elapsed, boost::mpi::maximum<double>());
          state.SetIterationTime(elapsed);
        }
        if (!valid())
        {
          petscCheckFailed = true;
          state.SkipWithError("Timed PETSc assembly changed operator or energy");
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
        double cells = owned;
        double maximum = owned;
        double minimum = owned;
        if constexpr (Distributed)
        {
          cells = boost::mpi::all_reduce(comm, owned, std::plus<double>());
          maximum = boost::mpi::all_reduce(comm, owned, boost::mpi::maximum<double>());
          minimum = boost::mpi::all_reduce(comm, owned, boost::mpi::minimum<double>());
        }
        MatInfo info;
        ierr = MatGetInfo(form.getOperator(), MAT_GLOBAL_SUM, &info);
        assert(ierr == PETSC_SUCCESS);
        if (ierr != PETSC_SUCCESS)
        {
          petscCheckFailed = true;
          state.SkipWithError("PETSc matrix metadata query failed");
          return;
        }
        const double ranks = Distributed ? comm.size() : 1;
        state.counters["cells"] = cells;
        state.counters["owned_cells_min"] = minimum;
        state.counters["owned_cells_max"] = maximum;
        state.counters["cell_imbalance"] = maximum * ranks / cells;
        state.counters["dofs"] = space.getSize();
        state.counters["nnz"] = info.nz_used;
        state.counters["ranks"] = ranks;
        state.counters["degree"] = Degree;
        state.counters["quadrature_order"] = 6;
        state.counters["components"] = Elasticity ? dim : 1;
        state.counters["cells_per_second"] =
          benchmark::Counter(cells, benchmark::Counter::kIsIterationInvariantRate);
      }
  };

  /// @brief Suppresses non-root reporting without suppressing benchmark execution.
  class SilentReporter final : public benchmark::BenchmarkReporter
  {
    public:
      bool ReportContext(const Context&) override
      {
        return true;
      }
      void ReportRuns(const std::vector<Run>&) override {}
  };

  template <size_t Degree, bool Distributed>
  void registerPETSc(const Context::MPI& context)
  {
    for (auto geometry : {Polytope::Type::Segment, Polytope::Type::Triangle,
           Polytope::Type::Quadrilateral, Polytope::Type::Tetrahedron,
           Polytope::Type::Pyramid, Polytope::Type::Hexahedron, Polytope::Type::Wedge})
    {
      const std::string suffix = "/P" + std::to_string(Degree) + "/" +
        std::string(Convergence::UniformGrid::getGeometryName(geometry));
      for (bool conductivity : {false, true})
      {
        const std::string name =
          std::string(conductivity ? "Conductivity" : "Poisson") + suffix;
        benchmark::RegisterBenchmark(name.c_str(),
          [&, geometry, conductivity](auto& state) {
            PETScPhysicsAssembly<Degree, false, Distributed>::run(
              state, geometry, conductivity, context);
          })
          ->Arg(3)
          ->Arg(5)
          ->Arg(9)
          ->Iterations(3)
          ->UseManualTime();
      }
      const std::string name = "LinearElasticity" + suffix;
      benchmark::RegisterBenchmark(name.c_str(),
        [&, geometry](auto& state) {
          PETScPhysicsAssembly<Degree, true, Distributed>::run(
            state, geometry, false, context);
        })
        ->Arg(3)
        ->Arg(5)
        ->Arg(9)
        ->Iterations(3)
        ->UseManualTime();
    }
  }
}

int main(int argc, char** argv)
{
  // Benchmark arguments belong to Google Benchmark, not PETSc's options database.
  if (PetscInitialize(nullptr, nullptr, nullptr, nullptr) != PETSC_SUCCESS)
    return 1;
  int result = 0;
  {
    boost::mpi::environment environment(argc, argv);
    boost::mpi::communicator world;
    Context::MPI context(environment, world);
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
    if (result != 0 && world.rank() == 0)
      std::cerr << "Local mode requires one rank; random interleaving and adaptive "
                   "warmup are unsupported.\n";
    if (result == 0)
    {
      using namespace Rodin::Tests::Benchmarks;
      if (distributed)
      {
        registerPETSc<1, true>(context);
        registerPETSc<2, true>(context);
      }
      else
      {
        registerPETSc<1, false>(context);
        registerPETSc<2, false>(context);
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
        benchmark::AddCustomContext("petsc_version",
          std::to_string(PETSC_VERSION_MAJOR) + "." +
            std::to_string(PETSC_VERSION_MINOR) + "." +
            std::to_string(PETSC_VERSION_SUBMINOR));
#ifdef RODIN_USE_OPENMP
        benchmark::AddCustomContext(
          "openmp_max_threads", std::to_string(omp_get_max_threads()));
#endif
        SilentReporter silent;
        const size_t cases = world.rank() == 0
          ? benchmark::RunSpecifiedBenchmarks()
          : benchmark::RunSpecifiedBenchmarks(&silent);
        result = cases == 0 || petscCheckFailed ? 1 : 0;
      }
      benchmark::Shutdown();
    }
    result = boost::mpi::all_reduce(world, result, boost::mpi::maximum<int>());
  }
  if (PetscFinalize() != PETSC_SUCCESS)
    return 1;
  return result;
}
