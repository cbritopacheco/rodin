/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/** @file @brief PETSc local/MPI mass and reaction operator/load assembly. */
#include <chrono>
#include <iostream>
#include <boost/mpi.hpp>
#include "Rodin/PETSc.h"
#include "Rodin/MPI/Variational.h"
#include "Rodin/MPI/Geometry/Sharder.h"
#include "Rodin/Geometry/BalancedCompactPartitioner.h"
#include "MassWorkload.h"
#ifdef RODIN_USE_OPENMP
#include <omp.h>
#endif

namespace Rodin::Tests::Benchmarks
{
  bool petscMassCheckFailed = false;

  /**
   * @brief Measures full operator/load assembly with aligned MPI collectives.
   * Architecture: explicit root partitioning precedes construction of fixed
   * spaces. Untimed matrix, load, and constant-field identities certify each
   * output. Three iterations align collectives; maximum-rank elapsed assembly
   * time excludes barriers, reductions and checks. Reaction coefficients are
   * doubled after timing to verify replacement and cache invalidation.
   */
  template <size_t K, bool Vector, MassFamily Family, bool Distributed>
  class PETScMassAssembly
  {
    public:
      static void run(benchmark::State& state, Geometry::Polytope::Type geometry,
        bool reaction, bool rhs, const Context::MPI& context)
      {
        using namespace Geometry;
        using namespace Variational;
        using Scalar = PetscScalar;
        using Workload = MassWorkload<K, Scalar, Vector, Family>;
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
              const size_t d = root.getDimension();
              root.getConnectivity().compute(d, d);
              root.getConnectivity().compute(d, 0);
              root.getConnectivity().compute(d, d - 1);
              root.getConnectivity().compute(d - 1, d);
              root.getConnectivity().compute(d - 1, 0);
              BalancedCompactPartitioner partitioner(root);
              partitioner.partition(static_cast<size_t>(comm.size()));
              sharder.shard(partitioner);
              sharder.scatter(0);
            }
            return sharder.gather(0);
          }
        }();
        const size_t dim = Polytope::Traits(geometry).getDimension();
        auto space = Workload::makeSpace(mesh, dim);
        PETSc::Variational::TrialFunction u(space);
        PETSc::Variational::TestFunction v(space);
        BilinearForm matrix(u, v);
        LinearForm load(v);
        Real scale = 1;
        const Scalar expected =
          Workload::configure(u, v, matrix, load, dim, reaction, scale);
        matrix.assemble();
        load.assemble();
        BilinearForm reference(u, v);
        LinearForm referenceLoad(v);
        Workload::template configure<true>(
          u, v, reference, referenceLoad, dim, reaction, scale);
        reference.assemble();
        const auto baselineLoad = load;
        PetscReal baselineNorm = 0;
        auto ierr = MatNorm(reference.getOperator(), NORM_FROBENIUS, &baselineNorm);
        assert(ierr == PETSC_SUCCESS);
        bool callsSucceeded = ierr == PETSC_SUCCESS;
        const auto valid = [&]() {
          const Scalar energy = matrix(u.getSolution(), u.getSolution());
          auto difference = matrix;
          auto error = MatAXPY(difference.getOperator(), -scale, reference.getOperator(),
            SAME_NONZERO_PATTERN);
          assert(error == PETSC_SUCCESS);
          bool ok = error == PETSC_SUCCESS && callsSucceeded;
          PetscReal matrixError = 0;
          error = MatNorm(difference.getOperator(), NORM_FROBENIUS, &matrixError);
          assert(error == PETSC_SUCCESS);
          ok = ok && error == PETSC_SUCCESS;
          auto residual = load;
          error = MatMult(
            matrix.getOperator(), u.getSolution().getData(), residual.getVector());
          assert(error == PETSC_SUCCESS);
          ok = ok && error == PETSC_SUCCESS;
          error = VecAXPY(residual.getVector(), -1, load.getVector());
          assert(error == PETSC_SUCCESS);
          ok = ok && error == PETSC_SUCCESS;
          PetscReal actionError = 0, loadNorm = 0, loadError = 0;
          error = VecNorm(residual.getVector(), NORM_2, &actionError);
          assert(error == PETSC_SUCCESS);
          ok = ok && error == PETSC_SUCCESS;
          error = VecNorm(load.getVector(), NORM_2, &loadNorm);
          assert(error == PETSC_SUCCESS);
          ok = ok && error == PETSC_SUCCESS;
          auto loadDifference = load;
          error = VecAXPY(loadDifference.getVector(), -scale, baselineLoad.getVector());
          assert(error == PETSC_SUCCESS);
          ok = ok && error == PETSC_SUCCESS;
          error = VecNorm(loadDifference.getVector(), NORM_2, &loadError);
          assert(error == PETSC_SUCCESS);
          ok = ok && error == PETSC_SUCCESS && std::isfinite(std::abs(energy)) &&
            std::abs(energy / (scale * expected) - Scalar(1)) < 1e-10 &&
            matrixError < 1e-11 * std::max(Real(1), scale * baselineNorm) &&
            actionError < 1e-11 * std::max(Real(1), loadNorm) &&
            loadError < 1e-11 * std::max(Real(1), loadNorm);
          if (!ok)
            std::cerr << "rank=" << comm.rank() << " energy=" << energy
                      << " expected=" << scale * expected
                      << " matrix_error=" << matrixError
                      << " action_error=" << actionError << " load_error=" << loadError
                      << '\n';
          if constexpr (Distributed)
            ok = boost::mpi::all_reduce(comm, ok, std::logical_and<bool>());
          return ok;
        };
        matrix.assemble();
        load.assemble();
        if (!valid())
        {
          petscMassCheckFailed = true;
          state.SkipWithError(
            "PETSc mass/reaction operator or projection-load oracle failed");
          return;
        }
        for (auto _ : state)
        {
          if constexpr (Distributed)
            comm.barrier();
          const auto start = std::chrono::steady_clock::now();
          if (rhs)
            load.assemble();
          else
            matrix.assemble();
          const auto stop = std::chrono::steady_clock::now();
          double elapsed = std::chrono::duration<double>(stop - start).count();
          if constexpr (Distributed)
            elapsed =
              boost::mpi::all_reduce(comm, elapsed, boost::mpi::maximum<double>());
          state.SetIterationTime(elapsed);
        }
        if (!valid())
        {
          petscMassCheckFailed = true;
          state.SkipWithError("Timed PETSc assembly changed operator or load");
        }
        if (reaction)
        {
          scale = 2;
          matrix.assemble();
          load.assemble();
          if (!valid())
          {
            petscMassCheckFailed = true;
            state.SkipWithError("PETSc coefficient update did not replace operator/load");
          }
        }
        double owned = 0;
        if constexpr (Distributed)
        {
          for (Index i = 0; i < mesh.getShard().getCellCount(); ++i)
          {
            if (mesh.getShard().isOwned(dim, i))
              ++owned;
          }
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
        ierr = MatGetInfo(matrix.getOperator(), MAT_GLOBAL_SUM, &info);
        assert(ierr == PETSC_SUCCESS);
        if (ierr != PETSC_SUCCESS)
        {
          petscMassCheckFailed = true;
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
        state.counters["degree"] = K;
        state.counters["components"] = Vector ? dim : 1;
        state.counters["quadrature_order"] = 8;
        state.counters["quadrature_points_per_cell"] =
          QF::PolytopeQuadratureFormula::get(8, geometry).getSize();
      }

      static void registerCases(const Context::MPI& context)
      {
        MassWorkload<K, PetscScalar, Vector, Family>::registerCases(
          [&context](auto& state, auto geometry, bool reaction, bool rhs) {
            run(state, geometry, reaction, rhs, context);
          },
          true);
      }
  };

  template <bool Vector, bool Distributed>
  void registerPETScMass(const Context::MPI& context)
  {
    PETScMassAssembly<1, Vector, MassFamily::H1, Distributed>::registerCases(context);
    PETScMassAssembly<2, Vector, MassFamily::H1, Distributed>::registerCases(context);
    PETScMassAssembly<3, Vector, MassFamily::H1, Distributed>::registerCases(context);
    PETScMassAssembly<1, Vector, MassFamily::P1, Distributed>::registerCases(context);
    PETScMassAssembly<0, Vector, MassFamily::P0, Distributed>::registerCases(context);
    PETScMassAssembly<0, Vector, MassFamily::P0g, Distributed>::registerCases(context);
  }

  class MassSilentReporter final : public benchmark::BenchmarkReporter
  {
    public:
      bool ReportContext(const Context&) override
      {
        return true;
      }
      void ReportRuns(const std::vector<Run>&) override {}
  };
}

int main(int argc, char** argv)
{
  using namespace Rodin;
  using namespace Rodin::Tests::Benchmarks;
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
      if (distributed)
      {
        registerPETScMass<false, true>(context);
        registerPETScMass<true, true>(context);
      }
      else
      {
        registerPETScMass<false, false>(context);
        registerPETScMass<true, false>(context);
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
#ifdef PETSC_USE_COMPLEX
        benchmark::AddCustomContext("petsc_scalar", "Complex");
#else
        benchmark::AddCustomContext("petsc_scalar", "Real");
#endif
#ifdef RODIN_USE_OPENMP
        benchmark::AddCustomContext(
          "openmp_max_threads", std::to_string(omp_get_max_threads()));
#endif
        benchmark::AddCustomContext("compiler", __VERSION__);
        benchmark::AddCustomContext("petsc_version",
          std::to_string(PETSC_VERSION_MAJOR) + "." +
            std::to_string(PETSC_VERSION_MINOR) + "." +
            std::to_string(PETSC_VERSION_SUBMINOR));
        MassSilentReporter silent;
        const auto cases = world.rank() == 0 ? benchmark::RunSpecifiedBenchmarks()
                                             : benchmark::RunSpecifiedBenchmarks(&silent);
        result = !cases || petscMassCheckFailed ? 1 : 0;
      }
      benchmark::Shutdown();
    }
    result = boost::mpi::all_reduce(world, result, boost::mpi::maximum<int>());
  }
  if (PetscFinalize() != PETSC_SUCCESS)
    return 1;
  return result;
}
