/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Isolates transformation and Jacobian evaluation for every geometry.
 */
#include <benchmark/benchmark.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational/H1.h>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace Rodin::Tests::Benchmarks
{
  template <size_t K>
  static void BM_ParametricEvaluation(benchmark::State& state, bool jacobian)
  {
    const auto geometry = static_cast<Polytope::Type>(state.range(0));
    const Polytope::Traits traits(geometry);
    // Compare native dimensions with 3D embeddings of lower-dimensional cells.
    const size_t physicalDimension = static_cast<size_t>(state.range(1));
    Variational::RealH1Element<K> fe(geometry);
    PointCloud nodes(physicalDimension, fe.getCount());
    for (size_t a = 0; a < fe.getCount(); ++a)
      for (size_t j = 0; j < physicalDimension; ++j)
        nodes(j, a) = j < traits.getDimension() ? fe.getNode(a)[j] : Real(1);
    ParametricTransformation transformation(std::move(nodes), fe);
    const Math::SpatialPoint reference = traits.getCentroid();
    Math::SpatialPoint x;
    Math::SpatialMatrix<Real> J;
    // Warm any shared finite-element basis tables before timing.
    transformation.transform(x, reference);
    transformation.jacobian(J, reference);
    for (auto _ : state)
    {
      benchmark::DoNotOptimize(&reference);
      if (jacobian)
      {
        transformation.jacobian(J, reference);
        benchmark::DoNotOptimize(&J);
      }
      else
      {
        transformation.transform(x, reference);
        benchmark::DoNotOptimize(&x);
      }
      benchmark::ClobberMemory();
    }
  }

  static void registerParametricBenchmarks()
  {
    using G = Polytope::Type;
    for (const auto& [geometry, label] :
      {std::pair{G::Point, "Point"}, {G::Segment, "Segment"}, {G::Triangle, "Triangle"},
        {G::Quadrilateral, "Quadrilateral"}, {G::Tetrahedron, "Tetrahedron"},
        {G::Hexahedron, "Hexahedron"}, {G::Pyramid, "Pyramid"}, {G::Wedge, "Wedge"}})
    {
      const std::string name = label;
      for (const size_t physicalDimension :
        {std::max(size_t(1), Polytope::Traits(geometry).getDimension()), size_t(3)})
      {
        // Three-dimensional cells need only one registration.
        const bool embedded = physicalDimension == 3;
        const std::string prefix = embedded ? "Parametric/" : "ParametricNative/";
        for (bool jacobian : {false, true})
        {
          const std::string operation = jacobian ? "Jacobian" : "Transform";
          benchmark::RegisterBenchmark((prefix + operation + "/P1/" + name).c_str(),
            &BM_ParametricEvaluation<1>, jacobian)
            ->Args({static_cast<int>(geometry), static_cast<int>(physicalDimension)});
          benchmark::RegisterBenchmark((prefix + operation + "/P2/" + name).c_str(),
            &BM_ParametricEvaluation<2>, jacobian)
            ->Args({static_cast<int>(geometry), static_cast<int>(physicalDimension)});
          benchmark::RegisterBenchmark((prefix + operation + "/P3/" + name).c_str(),
            &BM_ParametricEvaluation<3>, jacobian)
            ->Args({static_cast<int>(geometry), static_cast<int>(physicalDimension)});
          benchmark::RegisterBenchmark((prefix + operation + "/P4/" + name).c_str(),
            &BM_ParametricEvaluation<4>, jacobian)
            ->Args({static_cast<int>(geometry), static_cast<int>(physicalDimension)});
        }
        if (Polytope::Traits(geometry).getDimension() == 3)
          break;
      }
    }
  }

  const bool Registered = [] {
    registerParametricBenchmarks();
    return true;
  }();
}
