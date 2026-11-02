/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Matched tensor eikonal propagation with analytic accuracy counters.
 */

#include <benchmark/benchmark.h>
#include "Rodin/Eikonal/FIM.h"
#include "Rodin/Eikonal/FMM.h"
#include "Rodin/Variational/MatrixFunction.h"

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  /**
   * @brief Measures seeded propagation on matched metric-nonobtuse simplices.
   * @tparam Marching Whether to measure FMM instead of FIM.
   * @param state Dimensions, intervals, velocity ratio, and worker count.
   */
  template <bool Marching>
  void propagation(benchmark::State& state)
  {
    const size_t dimension = state.range(0);
    const size_t intervals = state.range(1);
    const Real ratio = state.range(2);
    const size_t threads = state.range(3);
    const auto geometry =
      dimension == 2 ? Polytope::Type::Triangle : Polytope::Type::Tetrahedron;
    auto mesh = dimension == 2
      ? LocalMesh::UniformGrid(geometry, {intervals + 1, intervals + 1})
      : LocalMesh::UniformGrid(geometry, {intervals + 1, intervals + 1, intervals + 1});
    mesh.scale(Real(1) / intervals);
    Math::SpatialMatrix<Real> transform =
      Math::SpatialMatrix<Real>::Identity(dimension, dimension);
    // An exactly specified orthogonal rotation followed by a velocity stretch.
    transform(0, 0) = Real(0.8) * ratio;
    transform(0, 1) = Real(-0.6);
    transform(1, 0) = Real(0.6) * ratio;
    transform(1, 1) = Real(0.8);
    Math::Vector<Real> exact(mesh.getVertexCount());
    IndexVector seeds;
    for (auto vertex = mesh.getVertex(); vertex; ++vertex)
    {
      const auto reference = vertex->getCoordinates();
      exact(vertex->getIndex()) = reference.norm();
      if (reference.squaredNorm() == 0)
        seeds.push_back(vertex->getIndex());
      mesh.setVertexCoordinates(vertex->getIndex(), transform * reference);
    }
    mesh.getConnectivity().compute(dimension, 0);
    mesh.getConnectivity().compute(0, dimension);
    mesh.getConnectivity().compute(0, 0);
    P1 space(mesh);
    GridFunction arrival(space), marchingReference(space), iterativeReference(space);
    MatrixFunction coefficient(transform * transform.transpose());
    Eikonal::FMM referenceMarching(marchingReference, coefficient);
    Eikonal::FIM referenceIterative(iterativeReference, coefficient);
    static constexpr Real LocalTolerance = Real(1e-6);
    static constexpr size_t IntegrationOrder = 8;
    referenceIterative.setLocalTolerance(LocalTolerance)
      .setIntegrationOrder(IntegrationOrder);
    Eikonal::FMM marching(arrival, coefficient);
    Eikonal::FIM iterative(arrival, coefficient);
    iterative.setLocalTolerance(LocalTolerance)
      .setIntegrationOrder(IntegrationOrder)
      .setThreadCount(threads);
    try
    {
      referenceMarching.seed(seeds).solve();
      referenceIterative.seed(seeds).solve();
      for (auto _ : state)
      {
        if constexpr (Marching)
          marching.seed(seeds).solve();
        else
          iterative.seed(seeds).solve();
        benchmark::DoNotOptimize(arrival.getData().data());
      }
      if constexpr (!Marching)
      {
        if (!(arrival.getData().array() == iterativeReference.getData().array()).all() ||
          iterative.getIterations() != referenceIterative.getIterations())
        {
          state.SkipWithError("FIM worker count changed the serial field or round count");
          return;
        }
      }
    }
    catch (const std::exception& error)
    {
      state.SkipWithError(error.what());
      return;
    }
    Real maximumError = 0, maximumDifference = 0;
    for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
    {
      const Index dof = space.getDOFs(0, vertex)[0];
      maximumError = std::max(maximumError, std::abs(arrival[dof] - exact(vertex)));
      maximumDifference =
        std::max(maximumDifference, std::abs(arrival[dof] - marchingReference[dof]));
    }
    state.counters["vertices"] = mesh.getVertexCount();
    state.counters["cells"] = mesh.getCellCount();
    state.counters["dimension"] = dimension;
    state.counters["velocity_ratio"] = ratio;
    state.counters["workers"] = threads;
    state.counters["rounds"] = Marching ? 0 : iterative.getIterations();
    state.counters["max_error"] = maximumError;
    state.counters["max_difference_fmm"] = maximumDifference;
    state.counters["search_resolution"] = Marching ? 0 : LocalTolerance;
    state.counters["path_order"] = Marching ? 0 : IntegrationOrder;
  }

  /**
   * @brief Registers matched isotropic and anisotropic cases in two and three dimensions.
   * @param cases Benchmark whose argument list is populated.
   * @param marching Whether only the serial FMM worker count applies.
   */
  void arguments(benchmark::internal::Benchmark* cases, bool marching)
  {
    for (size_t dimension : {2u, 3u})
    {
      for (size_t intervals :
        dimension == 2 ? std::vector<size_t>{16, 32} : std::vector<size_t>{4, 8})
      {
        for (size_t ratio : {1u, 2u})
        {
          cases->Args({static_cast<int64_t>(dimension), static_cast<int64_t>(intervals),
            static_cast<int64_t>(ratio), 1});
#ifdef RODIN_USE_OPENMP
          if (!marching)
          {
            for (size_t threads : {2u, 4u})
            {
              cases->Args(
                {static_cast<int64_t>(dimension), static_cast<int64_t>(intervals),
                  static_cast<int64_t>(ratio), static_cast<int64_t>(threads)});
            }
          }
#endif
        }
      }
    }
  }

  BENCHMARK_TEMPLATE(propagation, true)->Name("Eikonal/FMM")->Apply([](auto* cases) {
    arguments(cases, true);
  }) -> UseRealTime();
  BENCHMARK_TEMPLATE(propagation, false)->Name("Eikonal/FIM")->Apply([](auto* cases) {
    arguments(cases, false);
  }) -> UseRealTime();
}
