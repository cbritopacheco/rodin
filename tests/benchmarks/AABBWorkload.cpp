/*
 *          Copyright Carlos BRITO PACHECO 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @file
 * @brief Paired AABB workloads; see dev/profile_aabb.py for diagnostic counters.
 */

#include <Rodin/Location.h>
#include <Rodin/Variational.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Rodin;
using namespace Rodin::Geometry;

namespace
{
  using MeshType = Mesh<Context::Local>;
  using G = Polytope::Type;
  using Clock = std::chrono::steady_clock;
  /// At most 64 spatially distributed queries per class limit diagnostic cost.
  constexpr size_t QueryLimit = 64;
  /// Seven interleaved query groups cover every sample with alternating paired runs.
  /// Expensive groups run once; cheap groups repeat for at least BlockSeconds.
  [[maybe_unused]] constexpr size_t TimingBlocks = 7;
  /// Three paired constructions suffice for a median; expensive P4 builds can
  /// take seconds, so query groups and construction repetitions have separate budgets.
  [[maybe_unused]] constexpr size_t ConstructionRepetitions = 3;
  /// Each block lasts at least 3 ms; longer runs can scale this via the CLI.
  constexpr double BlockSeconds = 0.003;
  /// Face incidence uses roundoff-scale reference-coordinate comparisons.
  constexpr Real FaceTolerance = 1e-12;
  /// Nearby misses are displaced by this fraction of one cell width.
  constexpr Real MissCellFraction = 1e-4;
  /// Membership comparison is substantially tighter than fixture curvature.
  constexpr Real ComparisonTolerance = 1e-8;
#ifndef RODIN_AABB_WORKLOAD_DIAGNOSTICS
  volatile size_t checksum = 0;
#endif

  Math::SpatialPoint warp(Math::SpatialPoint x, size_t degree, Real distortion)
  {
    // A triangular shear preserves injectivity and volume in dimension >= 2.
    // In one dimension it is monotone for the nonnegative fixture coordinates.
    x[x.size() - 1] += distortion * std::pow(x[0], static_cast<Real>(degree));
    return x;
  }

  template <size_t K>
  void deform(MeshType& mesh, G type, Real distortion, size_t warpDegree)
  {
    Variational::RealH1Element<K> element(type);
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      PointCloud nodes(mesh.getDimension(), element.getCount());
      for (size_t a = 0; a < element.getCount(); ++a)
      {
        Math::SpatialPoint x;
        cell->getTransformation().transform(x, element.getNode(a));
        x = warp(x, warpDegree, distortion);
        for (size_t axis = 0; axis < mesh.getDimension(); ++axis)
          nodes(axis, a) = x[axis];
      }
      mesh.setPolytopeTransformation({mesh.getDimension(), cell->getIndex()},
        new ParametricTransformation<Variational::RealH1Element<K>>(
          std::move(nodes), element));
    }
  }

  struct Queries
  {
      std::vector<Math::SpatialPoint> interior;
      std::vector<Math::SpatialPoint> shared;
      std::vector<Math::SpatialPoint> miss;
  };

  Queries makeQueries(
    const MeshType& mesh, G type, size_t resolution, size_t degree, Real distortion)
  {
    Queries queries;
    if (type == G::Point)
    {
      queries.interior.push_back(mesh.getVertexCoordinates(0));
      auto outside = mesh.getVertexCoordinates(0);
      outside[0] += MissCellFraction;
      queries.miss.push_back(outside);
      return queries;
    }
    const Polytope::Traits traits(type);
    const auto& hs = traits.getHalfSpace();
    const size_t stride = std::max(size_t(1), mesh.getCellCount() / QueryLimit);
    for (auto cell = mesh.getCell(); cell; ++cell)
    {
      if (cell->getIndex() % stride != 0)
        continue;
      if (queries.interior.size() < QueryLimit)
      {
        Math::SpatialPoint x;
        const auto r =
          Real(0.75) * traits.getCentroid() + Real(0.25) * traits.getVertex(0);
        cell->getTransformation().transform(x, r);
        queries.interior.push_back(x);
      }
      if (queries.shared.size() >= QueryLimit)
        continue;
      for (Eigen::Index face = 0; face < hs.vector.size(); ++face)
      {
        Math::SpatialPoint r(mesh.getDimension());
        r.setZero();
        size_t vertices = 0;
        for (size_t v = 0; v < traits.getVertexCount(); ++v)
        {
          const auto& vertex = traits.getVertex(v);
          Real value = 0;
          for (size_t axis = 0; axis < mesh.getDimension(); ++axis)
            value += hs.matrix(face, axis) * vertex[axis];
          if (std::abs(value - hs.vector[face]) < FaceTolerance)
          {
            r += vertex;
            ++vertices;
          }
        }
        if (vertices == 0)
          throw std::runtime_error("Reference face has no vertices");
        r /= static_cast<Real>(vertices);
        Math::SpatialPoint x;
        cell->getTransformation().transform(x, r);
        // Undo the fixture shear to distinguish shared faces from domain faces.
        const Real first = x[0];
        if (mesh.getDimension() == 1)
          continue; // Segment shared vertices are generated separately below.
        auto original = x;
        original[original.size() - 1] -=
          distortion * std::pow(first, static_cast<Real>(degree));
        bool internal = true;
        for (size_t axis = 0; axis < mesh.getDimension(); ++axis)
          internal = internal && original[axis] > FaceTolerance &&
            original[axis] < Real(1) - FaceTolerance;
        if (internal)
        {
          queries.shared.push_back(x);
          break;
        }
      }
    }
    if (mesh.getDimension() == 1)
    {
      for (size_t v = 1; v + 1 < resolution && queries.shared.size() < QueryLimit; ++v)
      {
        Math::SpatialPoint x(1);
        x[0] = static_cast<Real>(v) / static_cast<Real>(resolution - 1);
        queries.shared.push_back(warp(x, degree, distortion));
      }
    }
    for (size_t i = 0; i < QueryLimit; ++i)
    {
      Math::SpatialPoint x(mesh.getDimension());
      x.setConstant(Real(0.5));
      // Stay away from corners: these misses probe the curved upper domain face.
      x[0] = Real(0.2) +
        Real(0.6) * (static_cast<Real>(i) + Real(0.5)) / static_cast<Real>(QueryLimit);
      x[x.size() - 1] = Real(1) + MissCellFraction / static_cast<Real>(resolution);
      queries.miss.push_back(warp(x, degree, distortion));
    }
    return queries;
  }

#ifndef RODIN_AABB_WORKLOAD_DIAGNOSTICS
  double median(std::vector<double> values)
  {
    std::sort(values.begin(), values.end());
    return values[values.size() / 2];
  }

  template <class Locator>
  size_t batch(const Locator& locator, const std::vector<Math::SpatialPoint>& queries)
  {
    size_t found = 0;
    for (const auto& x : queries)
    {
      const auto hit = locator.locate(x);
      found += hit.has_value();
      if (hit)
        found += hit->getPolytope().getIndex();
    }
    checksum = found;
    return found;
  }

#endif

  void run(G type, const char* name, size_t resolution, size_t degree, Real distortion,
    double seconds, size_t warpDegree)
  {
    const size_t dimension = Polytope::Traits(type).getDimension();
    Array<size_t> grid(dimension);
    grid.setConstant(resolution);
    MeshType mesh = type == G::Point
      ? MeshType::Builder().initialize(1).nodes(1).vertex({0.5}).finalize()
      : MeshType::UniformGrid(type, grid);
    Queries queries;
    if (warpDegree > degree)
      throw std::runtime_error("Warp power exceeds the geometry representation degree");
    if (type != G::Point)
    {
      mesh.scale(Real(1) / static_cast<Real>(resolution - 1));
      switch (degree)
      {
        case 1:
          deform<1>(mesh, type, distortion, warpDegree);
          break;
        case 2:
          deform<2>(mesh, type, distortion, warpDegree);
          break;
        case 4:
          deform<4>(mesh, type, distortion, warpDegree);
          break;
        default:
          throw std::runtime_error("Unsupported fixture degree");
      }
      queries = makeQueries(mesh, type, resolution, warpDegree, distortion);
    }
    else
      queries = makeQueries(mesh, type, resolution, warpDegree, distortion);
    Location::AABB off(mesh), on(mesh);
    off.setProjectionPruning(false);
    on.setProjectionPruning(true);
    // Warm both indices and the process-wide control conversion cache.
    off.locate(queries.interior.front());
    on.locate(queries.interior.front());
    const std::string prefix = std::string(name) + ',' + std::to_string(degree) + ',' +
      std::to_string(resolution) + ',' + std::to_string(mesh.getCellCount()) + ',' +
      std::to_string(distortion) + ',' + std::to_string(warpDegree);

    for (const auto& [label, values] : {std::pair{"interior", &queries.interior},
           std::pair{"shared", &queries.shared}, std::pair{"miss", &queries.miss}})
    {
      if (values->empty())
        continue;
      for (const auto& x : *values)
      {
        const auto a = off.locate(x), b = on.locate(x);
        if (a.has_value() != b.has_value() ||
          a.has_value() == (std::string(label) == "miss"))
          throw std::runtime_error(prefix + ": membership mismatch in " + label);
        if (a &&
          (a->getPolytope().getIndex() != b->getPolytope().getIndex() ||
            (a->getReferenceCoordinates() - b->getReferenceCoordinates()).norm() >
              ComparisonTolerance))
          throw std::runtime_error(prefix + ": returned point mismatch");
      }
#ifdef RODIN_AABB_WORKLOAD_DIAGNOSTICS
      for (bool enabled : {false, true})
      {
        auto& locator = enabled ? on : off;
        for (size_t i = 0; i < values->size(); ++i)
        {
          Location::aabbWork = {};
          locator.locate((*values)[i]);
          const auto& c = Location::aabbWork;
          std::cout << prefix << ',' << label << ',' << enabled << ',' << i << ','
                    << c.candidates << ',' << c.transforms << ',' << c.jacobians << ','
                    << c.iterations << ',' << c.retries << ',' << c.boxCandidates << ','
                    << c.projectionRejected << ',' << c.indexBytes << '\n';
        }
      }
#else
      double weightedTime[2] = {0, 0};
      for (size_t block = 0; block < TimingBlocks; ++block)
      {
        std::vector<Math::SpatialPoint> group;
        for (size_t i = block; i < values->size(); i += TimingBlocks)
          group.push_back((*values)[i]);
        if (group.empty())
          continue;
        for (size_t turn = 0; turn < 2; ++turn)
        {
          const size_t enabled = (block + turn) % 2;
          const auto& locator = enabled ? on : off;
          const auto start = Clock::now();
          const auto cpuStart = std::clock();
          size_t batches = 0;
          double elapsed;
          do
          {
            batch(locator, group);
            ++batches;
            elapsed = std::chrono::duration<double>(Clock::now() - start).count();
          } while (elapsed < seconds);
          const double cpuSeconds = static_cast<double>(std::clock() - cpuStart) /
            static_cast<double>(CLOCKS_PER_SEC);
          weightedTime[enabled] += cpuSeconds * 1e9 / static_cast<double>(batches);
        }
      }
      for (size_t enabled = 0; enabled < 2; ++enabled)
        std::cout << prefix << ',' << label << ',' << enabled << ',' << values->size()
                  << ',' << weightedTime[enabled] / static_cast<double>(values->size())
                  << '\n';
#endif
    }
#ifndef RODIN_AABB_WORKLOAD_DIAGNOSTICS
    std::vector<double> builds[2];
    for (size_t block = 0; block < ConstructionRepetitions; ++block)
      for (size_t turn = 0; turn < 2; ++turn)
      {
        const size_t enabled = (block + turn) % 2;
        const auto cpuStart = std::clock();
        {
          Location::AABB locator(mesh);
          locator.setProjectionPruning(enabled);
          batch(locator, {queries.interior.front()});
        }
        builds[enabled].push_back(static_cast<double>(std::clock() - cpuStart) * 1e9 /
          static_cast<double>(CLOCKS_PER_SEC));
      }
    for (size_t enabled = 0; enabled < 2; ++enabled)
      std::cout << prefix << ",build," << enabled << ",1," << median(builds[enabled])
                << '\n';
#endif
  }
}

int main(int argc, char** argv)
{
  try
  {
    const double seconds = argc > 1 ? std::stod(argv[1]) : BlockSeconds;
    const std::string filter = argc > 2 ? argv[2] : "";
    const std::string scenario = argc > 3 ? argv[3] : "";
    // Optional power holds the physical map fixed while changing FE order.
    const size_t warpDegree = argc > 4 ? std::stoul(argv[4]) : 0;
    if (!std::isfinite(seconds) || !(seconds > 0))
      throw std::runtime_error("Timing block duration must be positive");
    std::cout << std::setprecision(12);
#ifdef RODIN_AABB_WORKLOAD_DIAGNOSTICS
    std::cout
      << "geometry,degree,resolution,cells,distortion,warp_degree,query,pruning,sample,"
         "candidates,transforms,jacobians,iterations,retries,box_candidates,"
         "projection_rejected,index_bytes\n";
#else
    std::cout << "geometry,degree,resolution,cells,distortion,warp_degree,query,pruning,"
                 "samples,ns\n";
#endif
    for (const auto& [type, name] :
      {std::pair{G::Point, "Point"}, std::pair{G::Segment, "Segment"},
        std::pair{G::Triangle, "Triangle"}, std::pair{G::Quadrilateral, "Quadrilateral"},
        std::pair{G::Tetrahedron, "Tetrahedron"}, std::pair{G::Hexahedron, "Hexahedron"},
        std::pair{G::Pyramid, "Pyramid"}, std::pair{G::Wedge, "Wedge"}})
    {
      if (!filter.empty() && filter != name)
        continue;
      if (type == G::Point)
      {
        run(type, name, 1, 1, 0, seconds, 1);
        continue;
      }
      const size_t dim = Polytope::Traits(type).getDimension();
      // Resolution counts grid vertices per axis, not cells per axis.
      const std::vector<size_t> sizes = dim == 1 ? std::vector<size_t>{16, 64, 256}
        : dim == 2                               ? std::vector<size_t>{4, 8, 16}
                                                 : std::vector<size_t>{3, 5, 8};
      for (size_t degree : {1, 2, 4})
        for (Real distortion : {Real(0), Real(1), Real(4)})
          for (size_t size : sizes)
          {
            const std::string key = std::to_string(degree) + "/" +
              std::to_string(static_cast<int>(distortion)) + "/" + std::to_string(size);
            if (scenario.empty() || scenario == key)
              run(type, name, size, degree, distortion, seconds,
                warpDegree == 0 ? degree : warpDegree);
          }
    }
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
