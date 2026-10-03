/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Distributed finite-element-space construction benchmarks.
 *
 * Mesh creation is excluded from the timed region. Each iteration times the
 * slowest rank's P0 or vector P1 construction, including the ghost-DOF
 * exchange and reverse-index construction. Fixed iteration counts keep MPI
 * collectives aligned across ranks. Run with two or more MPI ranks.
 */

#include <chrono>

#include <benchmark/benchmark.h>
#include <boost/mpi/collectives.hpp>

#include <Rodin/MPI/Geometry/Mesh.h>
#include <Rodin/MPI/Variational/P0.h>
#include <Rodin/MPI/Variational/P1.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  template <class Range, template <class, class> class Space>
  void measureConstruction(
    benchmark::State& state, const Context::MPI& context, size_t vectorDimension)
  {
    const size_t n = static_cast<size_t>(state.range(0));
    const auto mesh =
      Mesh<Context::MPI>::UniformGrid(context, Polytope::Type::Triangle, {n, n});
    const auto& comm = context.getCommunicator();
    for (auto _ : state)
    {
      comm.barrier();
      const auto start = std::chrono::steady_clock::now();
      Space<Range, Mesh<Context::MPI>> fes(mesh, vectorDimension);
      const auto stop = std::chrono::steady_clock::now();
      benchmark::DoNotOptimize(fes.getSize());
      const double elapsed = std::chrono::duration<double>(stop - start).count();
      const double slowest =
        boost::mpi::all_reduce(comm, elapsed, boost::mpi::maximum<double>());
      state.SetIterationTime(slowest);
    }
    state.counters["local_cells"] = static_cast<double>(mesh.getShard().getCellCount());
    state.counters["ranks"] = static_cast<double>(comm.size());
  }
}

int main(int argc, char** argv)
{
  boost::mpi::environment environment(argc, argv);
  boost::mpi::communicator world;
  Context::MPI context(environment, world);

  auto add = [&](const char* name, auto function) {
    benchmark::RegisterBenchmark(name, function)
      ->Arg(32)
      ->Arg(64)
      ->Arg(128)
      ->Iterations(20)
      ->UseManualTime();
  };
  add("MPIConstruction/P0Scalar",
    [&](benchmark::State& state) { measureConstruction<Real, P0>(state, context, 1); });
  add("MPIConstruction/P0Vector2", [&](benchmark::State& state) {
    measureConstruction<Math::SpatialVector<Real>, P0>(state, context, 2);
  });
  add("MPIConstruction/P1Vector3", [&](benchmark::State& state) {
    measureConstruction<Math::SpatialVector<Real>, P1>(state, context, 3);
  });

  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
}
