/*
 * Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE or copy at https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file
 * @brief Paired local/subset/MPI AABB timings over all reference geometries.
 *
 * Mesh creation, geometric deformation, sharding and query generation are
 * excluded from timing. The synthetic partition assigns cells cyclically;
 * full overlap deliberately stresses filtering costs. Build includes owned
 * selection and first index construction. Queries compare the MPI wrapper
 * against the identical owned subset and the complete local shard. MPI
 * barriers/reductions surround each batch and are excluded from its timing.
 * Fixed iteration counts keep those benchmark-only collectives aligned.
 * Point has no nontrivial geometry order or distinct boundary query.
 */
#include <chrono>
#include <sstream>
#include <benchmark/benchmark.h>
#include <boost/mpi/collectives.hpp>
#include <Rodin/MPI/Location.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace
{
  using Clock = std::chrono::steady_clock;
  enum class Backend
  {
    FullShard,
    OwnedSubset,
    MPI
  };
  enum class Query
  {
    Build,
    Interior,
    Boundary,
    NearMiss
  };
  /// Small fixed batches amortize timing without dominating curved build runs.
  constexpr size_t QueryBatch = 16;
  /// Cheap queries/builds repeat for stable timing; expensive operations finish
  /// after one batch/build. This does not synchronize or constrain locator calls.
  [[maybe_unused]] constexpr double MinimumSampleSeconds = 0.005;
  /// Repeated manual samples; all ranks execute the same number of collectives.
  constexpr size_t TimingSamples = 5;
  /// Offset from a reference vertex, well outside numerical locator tolerance.
  constexpr Real MissOffset = Real(1e-4);
  /// Interior sample pulled toward the centroid to avoid partition boundaries.
  constexpr Real CentroidWeight = Real(0.75);
  /// Grid points per axis: two sizes expose scaling without oversized P4 runs.
  constexpr int64_t SmallResolution = 5;
  constexpr int64_t LargeResolution = 9;
  /// Distortion arguments use tenths; 0.2 is a moderate injective fixture shear.
  constexpr Real DistortionUnit = Real(0.1);
  constexpr int64_t DistortedFixture = 2;

  template <size_t K>
  void measure(benchmark::State& state, const Context::MPI& context, Polytope::Type type,
    Backend backend, Query category)
  {
    const auto& comm = context.getCommunicator();
    const size_t d = Polytope::Traits(type).getDimension();
    const size_t resolution = static_cast<size_t>(state.range(0));
    const bool overlap = state.range(1) != 0;
    const Real distortion = static_cast<Real>(state.range(2)) * DistortionUnit;
    auto parent = d == 0
      ? LocalMesh::Builder().initialize(1).nodes(1).vertex({Real(comm.rank())}).finalize()
      : LocalMesh::UniformGrid(type, Array<size_t>::Constant(d, resolution));
    if (d > 0)
    {
      Variational::RealH1Element<K> element(type);
      for (auto cell = parent.getCell(); cell; ++cell)
      {
        PointCloud nodes(d, element.getCount());
        for (size_t a = 0; a < element.getCount(); ++a)
        {
          Math::SpatialPoint x;
          cell->getTransformation().transform(x, element.getNode(a));
          // Global triangular shear (monotone in 1D) remains injective.
          x[d - 1] += distortion * x[0] * x[0];
          for (size_t j = 0; j < d; ++j)
            nodes(j, a) = x[j];
        }
        parent.setPolytopeTransformation({d, cell->getIndex()},
          new ParametricTransformation<Variational::RealH1Element<K>>(
            std::move(nodes), element));
      }
    }
    // Vertex owners are the smallest incident cell owner. This keeps lower-
    // dimensional ownership valid and measures the MPI wrapper's real selection
    // work, rather than marking every support vertex owned on every rank.
    std::vector<int> vertexOwner(parent.getVertexCount(), comm.size());
    std::vector<bool> partitionVertex(parent.getVertexCount(), false);
    std::vector<FlatSet<int>> holders(parent.getVertexCount());
    if (d > 0)
      for (auto cell = parent.getCell(); cell; ++cell)
      {
        const int owner = static_cast<int>(cell->getIndex() % comm.size());
        for (Index v : cell->getVertices())
        {
          vertexOwner[v] = std::min(vertexOwner[v], owner);
          partitionVertex[v] = partitionVertex[v] || owner == comm.rank();
          holders[v].insert(owner);
        }
      }
    Shard::Builder builder;
    builder.initialize(parent);
    for (auto cell = parent.getPolytope(d); cell; ++cell)
    {
      const int owner =
        d == 0 ? comm.rank() : static_cast<int>(cell->getIndex() % comm.size());
      const bool owned = owner == comm.rank();
      if (!owned && !overlap)
        continue;
      if (d > 0)
        for (Index v : cell->getVertices())
        {
          const auto vState = vertexOwner[v] == comm.rank() ? Shard::State::Owned
            : partitionVertex[v]                            ? Shard::State::Shared
                                                            : Shard::State::Ghost;
          const auto [local, inserted] = builder.include({0, v}, vState);
          if (!inserted)
            continue;
          if (vState != Shard::State::Owned)
            builder.setOwner(0, local, vertexOwner[v]);
          else
            for (int r = 0; r < comm.size(); ++r)
              if (r != comm.rank() && (overlap || holders[v].contains(r)))
                builder.halo(0, local, r);
        }
      const auto [local, inserted] = builder.include(
        {d, cell->getIndex()}, owned ? Shard::State::Owned : Shard::State::Ghost);
      if (!owned)
        builder.setOwner(d, local, owner);
      else if (overlap)
        for (int r = 0; r < comm.size(); ++r)
          if (r != comm.rank())
            builder.halo(d, local, r);
    }
    auto localShard = builder.finalize();
    if (d == 0)
    {
      // Each rank contributes one distinct point to the distributed point set.
      auto& ids = localShard.getPolytopeMap(0);
      ids.left[0] = comm.rank();
      ids.right.clear();
      ids.right.emplace(comm.rank(), 0);
      localShard.getHalo(0).clear();
    }
    auto mesh = MPIMesh::Builder(context).initialize(std::move(localShard)).finalize();
    const auto& shard = mesh.getShard();
    Location::AABB<LocalMesh>::Candidates candidates(d + 1);
    std::vector<Math::SpatialPoint> queries;
    for (auto cell = mesh.getPolytope(d); cell; ++cell)
    {
      if (!shard.isOwned(d, cell->getIndex()))
        continue;
      candidates[d].push_back(cell->getIndex());
      // Every backend repeats exactly the same query distribution, even when
      // adaptive batching gives cheap queries more repetitions.
      if (queries.size() >= QueryBatch)
        continue;
      const Polytope::Traits traits(type);
      Math::SpatialPoint rc = d == 0 ? Math::SpatialPoint(0) : traits.getVertex(0);
      if (d > 0 && category != Query::Boundary && category != Query::NearMiss)
        rc = CentroidWeight * traits.getCentroid() + (1 - CentroidWeight) * rc;
      Math::SpatialPoint x;
      cell->getTransformation().transform(x, rc);
      if (category == Query::NearMiss)
      {
        // All positive-dimensional fixture cells have nonnegative first
        // coordinates: probe just outside their common left domain boundary.
        x[0] = d == 0 ? x[0] + MissOffset : -MissOffset;
      }
      queries.push_back(x);
    }
    const bool anyEmpty =
      boost::mpi::all_reduce(comm, queries.empty(), std::logical_or<bool>());
    if (anyEmpty)
    {
      state.SkipWithError("Increase resolution so every rank owns a candidate");
      return;
    }
    Location::AABB<LocalMesh> full(shard), selected(shard, candidates);
    Location::AABB distributed(mesh);
    full.locate(d, queries.front());
    bool mismatch = false;
    for (const auto& x : queries)
    {
      const auto local = selected.locate(d, x);
      const auto lifted = distributed.locate(d, x);
      if (local.has_value() != lifted.has_value() ||
        (local && local->getPolytope().getIndex() != lifted->getPolytope().getIndex()) ||
        (category == Query::NearMiss && (local || full.locate(d, x))))
      {
        mismatch = true;
      }
    }
    if (boost::mpi::all_reduce(comm, mismatch, std::logical_or<bool>()))
    {
      state.SkipWithError("Local/MPI query mismatch or invalid miss fixture");
      return;
    }
    size_t checksum = 0;
    size_t queryOperations = 0;
#ifdef RODIN_AABB_WORKLOAD_DIAGNOSTICS
    Location::aabbWork = {};
#endif
    for (auto _ : state)
    {
      comm.barrier();
      const auto start = Clock::now();
      size_t operations = 1;
      if (category == Query::Build)
      {
        operations = 0;
        do
        {
          if (backend == Backend::MPI)
          {
            Location::AABB locator(mesh);
            checksum += locator.locate(d, queries.front()).has_value();
          }
          else if (backend == Backend::OwnedSubset)
          {
            Location::AABB<LocalMesh> locator(shard, candidates);
            checksum += locator.locate(d, queries.front()).has_value();
          }
          else
          {
            Location::AABB<LocalMesh> locator(shard);
            checksum += locator.locate(d, queries.front()).has_value();
          }
          ++operations;
        }
#ifdef RODIN_AABB_WORKLOAD_DIAGNOSTICS
        while (false);
#else
        while (std::chrono::duration<double>(Clock::now() - start).count() <
          MinimumSampleSeconds);
#endif
      }
      else
      {
        operations = 0;
        do
        {
          for (size_t q = 0; q < QueryBatch; ++q)
          {
            const auto& x = queries[(operations + q) % queries.size()];
            if (backend == Backend::MPI)
              checksum += distributed.locate(d, x).has_value();
            else if (backend == Backend::OwnedSubset)
              checksum += selected.locate(d, x).has_value();
            else
              checksum += full.locate(d, x).has_value();
          }
          operations += QueryBatch;
        }
#ifdef RODIN_AABB_WORKLOAD_DIAGNOSTICS
        while (false); // One batch suffices for counters; no diagnostic timing.
#else
        while (std::chrono::duration<double>(Clock::now() - start).count() <
          MinimumSampleSeconds);
#endif
      }
      const auto stop = Clock::now();
      benchmark::DoNotOptimize(checksum);
      queryOperations += operations;
      const double seconds =
        std::chrono::duration<double>(stop - start).count() / operations;
      state.SetIterationTime(
        boost::mpi::all_reduce(comm, seconds, boost::mpi::maximum<double>()));
    }
#ifdef RODIN_AABB_WORKLOAD_DIAGNOSTICS
    // Diagnostic builds count work; their timings are not performance evidence.
    const auto& work = Location::aabbWork;
    const double operations = queryOperations;
    state.counters["candidates_per_query"] = work.candidates / operations;
    state.counters["transforms_per_query"] = work.transforms / operations;
    state.counters["jacobians_per_query"] = work.jacobians / operations;
    state.counters["iterations_per_query"] = work.iterations / operations;
    state.counters["retries_per_query"] = work.retries / operations;
    state.counters["index_bytes"] = work.indexBytes;
#endif
    state.counters["queries"] = queryOperations;
    state.counters["stored"] = shard.getPolytopeCount(d);
    state.counters["owned"] = candidates[d].size();
    size_t candidateIds = backend == Backend::FullShard ? 0 : candidates[d].size();
    if (backend == Backend::MPI)
    {
      candidateIds = 0;
      for (size_t dimension = 0; dimension <= shard.getDimension(); ++dimension)
        for (Index i = 0; i < shard.getPolytopeCount(dimension); ++i)
          candidateIds += shard.isOwned(dimension, i);
    }
    state.counters["candidate_id_bytes"] = candidateIds * sizeof(Index);
    state.counters["ranks"] = comm.size();
  }

  template <size_t K>
  void registerOrder(const Context::MPI& context)
  {
    for (auto type : Polytope::Types)
    {
      if (type == Polytope::Type::Point && K > 1)
        continue;
      for (auto [backend, name] : {std::pair{Backend::FullShard, "FullShard"},
             {Backend::OwnedSubset, "OwnedSubset"}, {Backend::MPI, "MPI"}})
        for (auto [query, kind] :
          {std::pair{Query::Build, "Build"}, {Query::Interior, "Interior"},
            {Query::Boundary, "Boundary"}, {Query::NearMiss, "NearMiss"}})
        {
          std::ostringstream geometry;
          geometry << type;
          const auto label = std::string("MPIAABB/") + geometry.str() + "/P" +
            std::to_string(K) + "/" + name + "/" + kind;
          benchmark::RegisterBenchmark(label.c_str(),
            [&, type, backend, query](benchmark::State& state) {
              measure<K>(state, context, type, backend, query);
            })
            ->Args({SmallResolution, 0, 0})
            ->Args({SmallResolution, 1, 0})
            ->Args({SmallResolution, 0, DistortedFixture})
            ->Args({SmallResolution, 1, DistortedFixture})
            ->Args({LargeResolution, 0, DistortedFixture})
            ->Args({LargeResolution, 1, DistortedFixture})
            ->Iterations(TimingSamples)
            ->UseManualTime();
        }
    }
  }
}
int main(int argc, char** argv)
{
  boost::mpi::environment environment(argc, argv);
  boost::mpi::communicator world;
  Context::MPI context(environment, world);
  registerOrder<1>(context);
  registerOrder<2>(context);
  registerOrder<3>(context);
  registerOrder<4>(context);
  // Every rank writes its own result file; timings contain the slowest rank,
  // while diagnostic counters describe the rank that produced that file.
  std::vector<std::string> arguments;
  for (int i = 0; i < argc; ++i)
  {
    std::string arg = argv[i];
    if (world.rank() != 0 && arg.starts_with("--benchmark_out="))
      arg += ".rank" + std::to_string(world.rank());
    arguments.push_back(std::move(arg));
  }
  std::vector<char*> pointers;
  for (auto& arg : arguments)
    pointers.push_back(arg.data());
  pointers.push_back(nullptr);
  benchmark::Initialize(&argc, pointers.data());
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
}
