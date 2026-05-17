/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include <benchmark/benchmark.h>

#include <Rodin/Adaptation.h>
#include <Rodin/Assembly/Default.h>
#include <Rodin/Geometry.h>
#include <Rodin/Geometry/LevelSetDiscretizerTriangles.h>
#include <Rodin/Geometry/TriangleMeshOptimizer.h>
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Solver/SparseLU.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Adaptation::TargetMatrixOptimization;
using namespace Rodin::Geometry;
using namespace Rodin::Solver;
using namespace Rodin::Variational;

namespace Rodin::Tests::Benchmarks
{
  namespace
  {
    constexpr Attribute Interface = 30;
    constexpr Attribute Boundary = 40;
    constexpr Attribute Negative = 1;
    constexpr Attribute Positive = 2;
    constexpr Real Pi = 3.14159265358979323846;

    enum class ShapeCase : int
    {
      CircleOrbit = 0,
      SquareOrbit = 1,
      TriangleOrbit = 2,
      CircleEnterLeave = 3
    };

    enum class MetricCase : int
    {
      SquaredDistance = 0,
      AreaDistortion = 1,
      ShapeSize = 2,
      ShapeDistortion = 3
    };

    struct MeshStats
    {
      Real minQuality = 0;
      Index inverted = 0;
      Real coverage = 0;
      Real signedArea = 0;
    };

    struct StageCounters
    {
      Index cutCells = 0;
      Index finalCells = 0;
      Index interfaceEdges = 0;
      Index uncutCells = 0;
      Index snappedCrossings = 0;
      Index pathologicalCuts = 0;
      Index splits = 0;
      Index collapses = 0;
      Index swaps = 0;
      Real cutMinQuality = 0;
      Real finalMinQuality = 0;
      Real fitRms = 0;
      Real fitMax = 0;
      Real coverage = 0;
      Real signedArea = 0;
      Real maxInterfaceDeviation = 0;
      Index inverted = 0;
      Index tmopFailures = 0;
    };

    Math::SpatialPoint orbitCenter(Real t)
    {
      return Math::SpatialPoint{
        Real(0.5) + Real(0.14) * std::cos(Real(2) * Pi * t),
        Real(0.5) + Real(0.11) * std::sin(Real(2) * Pi * t)
      };
    }

    Math::SpatialPoint enterLeaveCenter(Real t)
    {
      return Math::SpatialPoint{
        Real(-0.20) + Real(1.40) * t,
        Real(0.50) + Real(0.08) * std::sin(Real(2) * Pi * t)
      };
    }

    Real circleValue(
        const Math::SpatialPoint& x,
        const Math::SpatialPoint& c,
        Real radius)
    {
      return std::hypot(x[0] - c[0], x[1] - c[1]) - radius;
    }

    Math::SpatialPoint circleGradient(
        const Math::SpatialPoint& x,
        const Math::SpatialPoint& c)
    {
      const Real dx = x[0] - c[0];
      const Real dy = x[1] - c[1];
      const Real r = std::hypot(dx, dy);
      if (r <= Real(1e-14))
        return Math::SpatialPoint{0, 0};
      return Math::SpatialPoint{dx / r, dy / r};
    }

    Real squareValue(const Math::SpatialPoint& x, Real t)
    {
      const auto c = orbitCenter(t);
      const Real half = Real(0.16);
      const Real ax = std::abs(x[0] - c[0]) - half;
      const Real ay = std::abs(x[1] - c[1]) - half;
      const Real ox = std::max(ax, Real(0));
      const Real oy = std::max(ay, Real(0));
      const Real outside = std::hypot(ox, oy);
      const Real inside = std::min(std::max(ax, ay), Real(0));
      return outside + inside;
    }

    Math::SpatialPoint squareGradient(const Math::SpatialPoint& x, Real t)
    {
      const auto c = orbitCenter(t);
      const Real half = Real(0.16);
      const Real dx = x[0] - c[0];
      const Real dy = x[1] - c[1];
      const Real ax = std::abs(dx) - half;
      const Real ay = std::abs(dy) - half;
      const Real ox = std::max(ax, Real(0));
      const Real oy = std::max(ay, Real(0));
      const Real n = std::hypot(ox, oy);
      const Real sx = dx >= Real(0) ? Real(1) : Real(-1);
      const Real sy = dy >= Real(0) ? Real(1) : Real(-1);
      if (n > Real(1e-14))
        return Math::SpatialPoint{sx * ox / n, sy * oy / n};
      if (ax >= ay)
        return Math::SpatialPoint{sx, 0};
      return Math::SpatialPoint{0, sy};
    }

    std::array<Math::SpatialPoint, 3> triangleVertices(Real t)
    {
      const auto c = orbitCenter(t);
      return {{
        Math::SpatialPoint{c[0], c[1] + Real(0.20)},
        Math::SpatialPoint{c[0] - Real(0.18), c[1] - Real(0.12)},
        Math::SpatialPoint{c[0] + Real(0.18), c[1] - Real(0.12)}
      }};
    }

    Real triangleValue(const Math::SpatialPoint& x, Real t)
    {
      const auto vertices = triangleVertices(t);
      Real phi = -std::numeric_limits<Real>::infinity();
      for (std::uint8_t i = 0; i < 3; ++i)
      {
        const auto& a = vertices[i];
        const auto& b = vertices[(i + 1) % 3];
        const Real ex = b[0] - a[0];
        const Real ey = b[1] - a[1];
        const Real len = std::hypot(ex, ey);
        const Real signedEdge =
          -(ex * (x[1] - a[1]) - ey * (x[0] - a[0])) / len;
        phi = std::max(phi, signedEdge);
      }
      return phi;
    }

    Math::SpatialPoint triangleGradient(const Math::SpatialPoint& x, Real t)
    {
      const auto vertices = triangleVertices(t);
      Real phi = -std::numeric_limits<Real>::infinity();
      Math::SpatialPoint gradient{1, 0};
      for (std::uint8_t i = 0; i < 3; ++i)
      {
        const auto& a = vertices[i];
        const auto& b = vertices[(i + 1) % 3];
        const Real ex = b[0] - a[0];
        const Real ey = b[1] - a[1];
        const Real len = std::hypot(ex, ey);
        const Real signedEdge =
          -(ex * (x[1] - a[1]) - ey * (x[0] - a[0])) / len;
        if (signedEdge > phi)
        {
          phi = signedEdge;
          gradient = Math::SpatialPoint{ey / len, -ex / len};
        }
      }
      return gradient;
    }

    Real levelSetValue(ShapeCase shape, const Math::SpatialPoint& x, Real t)
    {
      switch (shape)
      {
        case ShapeCase::CircleOrbit:
          return circleValue(x, orbitCenter(t), Real(0.18));
        case ShapeCase::SquareOrbit:
          return squareValue(x, t);
        case ShapeCase::TriangleOrbit:
          return triangleValue(x, t);
        case ShapeCase::CircleEnterLeave:
          return circleValue(x, enterLeaveCenter(t), Real(0.18));
      }
      return 1;
    }

    Math::SpatialPoint levelSetGradient(
        ShapeCase shape,
        const Math::SpatialPoint& x,
        Real t)
    {
      switch (shape)
      {
        case ShapeCase::CircleOrbit:
          return circleGradient(x, orbitCenter(t));
        case ShapeCase::SquareOrbit:
          return squareGradient(x, t);
        case ShapeCase::TriangleOrbit:
          return triangleGradient(x, t);
        case ShapeCase::CircleEnterLeave:
          return circleGradient(x, enterLeaveCenter(t));
      }
      return Math::SpatialPoint{0, 0};
    }

    Real boxBoundaryValue(const Math::SpatialPoint& x)
    {
      const std::array<Real, 4> d = {{
        std::abs(x[0]), std::abs(x[0] - Real(1)),
        std::abs(x[1]), std::abs(x[1] - Real(1))
      }};
      const auto side =
        static_cast<size_t>(std::min_element(d.begin(), d.end()) - d.begin());
      if (side == 0)
        return x[0];
      if (side == 1)
        return x[0] - Real(1);
      if (side == 2)
        return x[1];
      return x[1] - Real(1);
    }

    Math::SpatialPoint boxBoundaryGradient(const Math::SpatialPoint& x)
    {
      const std::array<Real, 4> d = {{
        std::abs(x[0]), std::abs(x[0] - Real(1)),
        std::abs(x[1]), std::abs(x[1] - Real(1))
      }};
      const auto side =
        static_cast<size_t>(std::min_element(d.begin(), d.end()) - d.begin());
      if (side < 2)
        return Math::SpatialPoint{1, 0};
      return Math::SpatialPoint{0, 1};
    }

    Real signedTriangleArea(
        const Math::SpatialPoint& a,
        const Math::SpatialPoint& b,
        const Math::SpatialPoint& c)
    {
      return Real(0.5)
        * ((b[0] - a[0]) * (c[1] - a[1])
         - (b[1] - a[1]) * (c[0] - a[0]));
    }

    Real triangleQuality(
        const Math::SpatialPoint& a,
        const Math::SpatialPoint& b,
        const Math::SpatialPoint& c)
    {
      const Real area = std::abs(signedTriangleArea(a, b, c));
      const Real den =
        (b - a).squaredNorm() + (c - b).squaredNorm()
      + (a - c).squaredNorm();
      if (den <= Real(0))
        return Real(0);
      return Real(4) * std::sqrt(Real(3)) * area / den;
    }

    MeshStats meshStats(const LocalMesh& mesh)
    {
      MeshStats stats;
      stats.minQuality = std::numeric_limits<Real>::infinity();
      const auto& conn = mesh.getConnectivity();
      for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
      {
        const auto& cell = conn.getPolytope(2, c);
        const auto x0 = mesh.getVertexCoordinates(cell(0));
        const auto x1 = mesh.getVertexCoordinates(cell(1));
        const auto x2 = mesh.getVertexCoordinates(cell(2));
        const Real area = signedTriangleArea(x0, x1, x2);
        if (area <= Real(0))
          stats.inverted++;
        stats.coverage += std::abs(area);
        stats.signedArea += area;
        stats.minQuality =
          std::min(stats.minQuality, triangleQuality(x0, x1, x2));
      }
      if (!std::isfinite(stats.minQuality))
        stats.minQuality = 0;
      return stats;
    }

    void annotateBoundary(LocalMesh& mesh)
    {
      const auto& conn = mesh.getConnectivity();
      for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
      {
        const auto& edge = conn.getPolytope(1, e);
        const auto x0 = mesh.getVertexCoordinates(edge(0));
        const auto x1 = mesh.getVertexCoordinates(edge(1));
        const Real eps = Real(1e-12);
        const bool on0 =
          std::abs(x0[0]) <= eps || std::abs(x0[0] - Real(1)) <= eps
       || std::abs(x0[1]) <= eps || std::abs(x0[1] - Real(1)) <= eps;
        const bool on1 =
          std::abs(x1[0]) <= eps || std::abs(x1[0] - Real(1)) <= eps
       || std::abs(x1[1]) <= eps || std::abs(x1[1] - Real(1)) <= eps;
        if (on0 && on1)
          mesh.setAttribute({ 1, e }, Boundary);
      }
    }

    LocalMesh makeBackground(size_t resolution)
    {
      auto mesh =
        LocalMesh::UniformGrid(Polytope::Type::Triangle, {resolution, resolution});
      mesh.scale(Real(1) / static_cast<Real>(resolution - 1));
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      return mesh;
    }

    LevelSetDiscretizerTrianglesResult cutShape(
        LocalMesh& background,
        ShapeCase shape,
        Real t,
        Real snapTolerance,
        Real minCutQuality)
    {
      background.getConnectivity().compute(2, 1);
      background.getConnectivity().compute(1, 0);

      P1<Real, LocalMesh> phiSpace(background);
      GridFunction phi(phiSpace);
      for (Index v = 0; v < background.getVertexCount(); ++v)
        phi[v] = levelSetValue(shape, background.getVertexCoordinates(v), t);

      auto cut = LevelSetDiscretizerTriangles(phi)
        .setSignTolerance(1e-12)
        .setInterfaceAttribute(Interface)
        .setNegativeCellAttribute(Negative)
        .setPositiveCellAttribute(Positive)
        .setCrossingSnapTolerance(snapTolerance)
        .setMinCutQuality(minCutQuality)
        .discretize();
      annotateBoundary(cut.mesh);
      return cut;
    }

    std::pair<Real, Real> interfaceFit(
        const LocalMesh& mesh,
        ShapeCase shape,
        Real t)
    {
      const auto& conn = mesh.getConnectivity();
      std::vector<char> onInterface(mesh.getVertexCount(), 0);
      for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
      {
        const auto attr = mesh.getAttribute(1, e);
        if (!attr || *attr != Interface)
          continue;
        const auto& edge = conn.getPolytope(1, e);
        onInterface[edge(0)] = 1;
        onInterface[edge(1)] = 1;
      }
      Real max = 0;
      Real sq = 0;
      Index count = 0;
      for (Index v = 0; v < static_cast<Index>(mesh.getVertexCount()); ++v)
      {
        if (!onInterface[v])
          continue;
        const Real d = std::abs(
            levelSetValue(shape, mesh.getVertexCoordinates(v), t));
        max = std::max(max, d);
        sq += d * d;
        count++;
      }
      return {
        count > 0 ? std::sqrt(sq / static_cast<Real>(count)) : Real(0),
        max
      };
    }

    template <class Metric>
    bool solveTMOPWithMetric(
        LocalMesh& mesh,
        ShapeCase shape,
        Real t,
        Metric metric,
        Real qualityWeight,
        Index maxIterations)
    {
      P1 space(mesh, mesh.getSpaceDimension());
      GridFunction u(space);
      u.getData().setZero();
      TrialFunction du(space);
      TestFunction v(space);

      auto phiValue = [shape, t](const Math::SpatialPoint& x)
      {
        return levelSetValue(shape, x, t);
      };
      auto phiGradient = [shape, t](const Math::SpatialPoint& x)
      {
        return levelSetGradient(shape, x, t);
      };

      QualityTerm quality(
          metric,
          OrientedEquilateralSameAreaTargetJacobian(mesh),
          qualityWeight);
      DeviationTerm deviation(1.0);
      AnalyticLevelSetFitTerm fit(
          phiValue, phiGradient, Optional<Attribute>(Interface), 1.0);
      AnalyticLevelSetFitTerm boundaryFit(
          boxBoundaryValue, boxBoundaryGradient, Optional<Attribute>(Boundary), 1.0);

      Variational::Problem problem(du, v);
      problem =
          quality.tangent(u, du, v)
        + deviation.tangent(du, v)
        + fit.tangent(u, du, v)
        + boundaryFit.tangent(u, du, v)
        + quality.residual(u, v)
        + deviation.residual(u, v)
        + fit.residual(u, v)
        + boundaryFit.residual(u, v);

      SparseLU solver{problem};
      NewtonSolver newton(solver);
      newton.setMaxIterations(maxIterations)
            .setDampingFactor(1.0)
            .setAbsoluteTolerance(1e-10)
            .setRelativeTolerance(1e-8);
      try
      {
        newton.solve(u);
      }
      catch (...)
      {
        return false;
      }

      const auto& data = u.getData();
      for (Eigen::Index i = 0; i < data.size(); ++i)
      {
        if (!std::isfinite(data(i)))
          return false;
      }

      const Index nv = mesh.getVertexCount();
      std::vector<Math::SpatialPoint> originalVertices;
      originalVertices.reserve(static_cast<size_t>(nv));
      for (Index vtx = 0; vtx < nv; ++vtx)
        originalVertices.push_back(mesh.getVertexCoordinates(vtx));

      for (Index vtx = 0; vtx < nv; ++vtx)
      {
        auto x = mesh.getVertexCoordinates(vtx);
        x[0] += data(vtx);
        x[1] += data(vtx + nv);
        mesh.setVertexCoordinates(vtx, x);
      }
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);

      const auto stats = meshStats(mesh);
      const bool valid =
        stats.inverted == 0
        && std::isfinite(stats.coverage)
        && std::isfinite(stats.signedArea)
        && stats.coverage > Real(0.8)
        && stats.coverage < Real(1.2);
      if (!valid)
      {
        for (Index vtx = 0; vtx < nv; ++vtx)
          mesh.setVertexCoordinates(vtx, originalVertices[vtx]);
        mesh.getConnectivity().compute(2, 1);
        mesh.getConnectivity().compute(1, 0);
        return false;
      }
      return true;
    }

    bool solveTMOP(
        LocalMesh& mesh,
        ShapeCase shape,
        Real t,
        MetricCase metricCase = MetricCase::SquaredDistance)
    {
      switch (metricCase)
      {
        case MetricCase::SquaredDistance:
          return solveTMOPWithMetric(
              mesh, shape, t, SquaredDistanceMetric{}, 0.05, 5);
        case MetricCase::AreaDistortion:
          return solveTMOPWithMetric(
              mesh, shape, t, AreaDistortionMetric{}, 0.05, 5);
        case MetricCase::ShapeSize:
          return solveTMOPWithMetric(
              mesh, shape, t, ShapeSizeMetric{}, 0.03, 5);
        case MetricCase::ShapeDistortion:
          // The pure shape barrier has a scale nullspace and still lacks
          // line-searched step rejection in NewtonSolver; keep it in the
          // benchmark matrix, but with a conservative one-step diagnostic.
          return solveTMOPWithMetric(
              mesh, shape, t, ShapeDistortionMetric{}, 0.02, 1);
      }
      return false;
    }

    bool solveSquaredTMOP(
        LocalMesh& mesh,
        ShapeCase shape,
        Real t,
        Real qualityWeight)
    {
      return solveTMOPWithMetric(
          mesh, shape, t, SquaredDistanceMetric{}, qualityWeight, 5);
    }

    TriangleMeshOptimizerReport coarsen(LocalMesh& mesh, size_t resolution)
    {
      std::vector<char> frozen(mesh.getVertexCount(), 0);
      const auto& conn = mesh.getConnectivity();
      for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
      {
        const auto attr = mesh.getAttribute(1, e);
        if (!attr || *attr != Interface)
          continue;
        const auto& edge = conn.getPolytope(1, e);
        frozen[edge(0)] = 1;
        frozen[edge(1)] = 1;
      }

      return TriangleMeshOptimizer()
        .setHMin(Real(0.35) / static_cast<Real>(resolution - 1))
        .setHMax(std::numeric_limits<Real>::infinity())
        .setMinQuality(0.1)
        .setMaxIterations(4)
        .setProtectedVertices(std::move(frozen))
        .optimize(mesh);
    }

    StageCounters setCounters(
        benchmark::State& st,
        const LevelSetDiscretizerTrianglesResult& cut,
        const LocalMesh& finalMesh,
        const TriangleMeshOptimizerReport* optReport,
        ShapeCase shape,
        bool tmopSucceeded,
        Real t)
    {
      StageCounters counters;
      const auto stats = meshStats(finalMesh);
      const auto fit = interfaceFit(finalMesh, shape, t);
      counters.cutCells = cut.mesh.getCellCount();
      counters.finalCells = finalMesh.getCellCount();
      counters.interfaceEdges = cut.report.interfaceEdgeProvenance.size();
      counters.uncutCells = cut.report.uncutCellCount;
      counters.snappedCrossings = cut.report.snappedCrossingCount;
      counters.pathologicalCuts = cut.report.pathologicalCutCount;
      counters.cutMinQuality = cut.report.minOutputCellQuality;
      counters.finalMinQuality = stats.minQuality;
      counters.fitRms = fit.first;
      counters.fitMax = fit.second;
      counters.coverage = stats.coverage;
      counters.signedArea = stats.signedArea;
      counters.maxInterfaceDeviation = cut.report.maxInterfaceDeviation;
      counters.inverted = stats.inverted;
      counters.tmopFailures = tmopSucceeded ? 0 : 1;
      if (optReport)
      {
        counters.splits = optReport->splits;
        counters.collapses = optReport->collapses;
        counters.swaps = optReport->swaps;
      }

      st.counters["cut_cells"] = benchmark::Counter(counters.cutCells);
      st.counters["final_cells"] = benchmark::Counter(counters.finalCells);
      st.counters["interface_edges"] = benchmark::Counter(counters.interfaceEdges);
      st.counters["uncut_cells"] = benchmark::Counter(counters.uncutCells);
      st.counters["snapped"] = benchmark::Counter(counters.snappedCrossings);
      st.counters["pathological"] = benchmark::Counter(counters.pathologicalCuts);
      st.counters["qmin_cut"] = benchmark::Counter(counters.cutMinQuality);
      st.counters["qmin_final"] = benchmark::Counter(counters.finalMinQuality);
      st.counters["fit_rms"] = benchmark::Counter(counters.fitRms);
      st.counters["fit_max"] = benchmark::Counter(counters.fitMax);
      st.counters["coverage"] = benchmark::Counter(counters.coverage);
      st.counters["signed_area"] = benchmark::Counter(counters.signedArea);
      st.counters["max_iface_dev"] =
        benchmark::Counter(counters.maxInterfaceDeviation);
      st.counters["inverted"] = benchmark::Counter(counters.inverted);
      st.counters["tmop_failed"] = benchmark::Counter(counters.tmopFailures);
      st.counters["splits"] = benchmark::Counter(counters.splits);
      st.counters["collapses"] = benchmark::Counter(counters.collapses);
      st.counters["swaps"] = benchmark::Counter(counters.swaps);
      return counters;
    }

    void demoteInterface(LocalMesh& mesh)
    {
      auto& conn = mesh.getConnectivity();
      for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
      {
        const auto attr = mesh.getAttribute(1, e);
        if (attr && *attr == Interface)
          mesh.setAttribute({1, e}, Attribute(0));
      }
    }
  }

  static void BM_LevelSetPipeline_CutExact(benchmark::State& st)
  {
    const auto resolution = static_cast<size_t>(st.range(0));
    const auto shape = static_cast<ShapeCase>(st.range(1));
    StageCounters counters;
    for (auto _ : st)
    {
      auto background = makeBackground(resolution);
      auto cut = cutShape(background, shape, Real(0.25), Real(0), Real(0));
      benchmark::DoNotOptimize(cut.mesh.getCellCount());
      counters =
        setCounters(st, cut, cut.mesh, nullptr, shape, true, Real(0.25));
    }
    st.SetItemsProcessed(st.iterations() * counters.cutCells);
  }

  static void BM_LevelSetPipeline_CutQualityAware(benchmark::State& st)
  {
    const auto resolution = static_cast<size_t>(st.range(0));
    const auto shape = static_cast<ShapeCase>(st.range(1));
    StageCounters counters;
    for (auto _ : st)
    {
      auto background = makeBackground(resolution);
      auto cut = cutShape(background, shape, Real(0.25), Real(0.05), Real(0.2));
      benchmark::DoNotOptimize(cut.mesh.getCellCount());
      counters =
        setCounters(st, cut, cut.mesh, nullptr, shape, true, Real(0.25));
    }
    st.SetItemsProcessed(st.iterations() * counters.cutCells);
  }

  static void BM_LevelSetPipeline_TMOPMetricFreshCut(benchmark::State& st)
  {
    const auto resolution = static_cast<size_t>(st.range(0));
    const auto shape = static_cast<ShapeCase>(st.range(1));
    const auto metric = static_cast<MetricCase>(st.range(2));
    StageCounters counters;
    for (auto _ : st)
    {
      auto background = makeBackground(resolution);
      auto cut = cutShape(background, shape, Real(0.25), Real(0.05), Real(0.2));
      auto optimized = cut.mesh;
      const bool tmopOk = solveTMOP(optimized, shape, Real(0.25), metric);
      benchmark::DoNotOptimize(optimized.getCellCount());
      counters =
        setCounters(st, cut, optimized, nullptr, shape, tmopOk, Real(0.25));
    }
    st.SetItemsProcessed(st.iterations() * counters.cutCells);
  }

  static void BM_LevelSetPipeline_CoarsenAfterTMOP(benchmark::State& st)
  {
    const auto resolution = static_cast<size_t>(st.range(0));
    const auto shape = static_cast<ShapeCase>(st.range(1));
    StageCounters counters;
    for (auto _ : st)
    {
      auto background = makeBackground(resolution);
      auto cut = cutShape(background, shape, Real(0.25), Real(0.05), Real(0.2));
      auto optimized = cut.mesh;
      const bool tmopOk =
        solveTMOP(optimized, shape, Real(0.25), MetricCase::SquaredDistance);
      const auto report = coarsen(optimized, resolution);
      benchmark::DoNotOptimize(optimized.getCellCount());
      counters =
        setCounters(st, cut, optimized, &report, shape, tmopOk, Real(0.25));
    }
    st.SetItemsProcessed(st.iterations() * counters.finalCells);
  }

  static void BM_LevelSetPipeline_TMOPSquaredWeightFreshCut(
      benchmark::State& st)
  {
    const auto resolution = static_cast<size_t>(st.range(0));
    const auto shape = static_cast<ShapeCase>(st.range(1));
    const Real qualityWeight =
      static_cast<Real>(st.range(2)) / static_cast<Real>(1000);
    StageCounters counters;
    for (auto _ : st)
    {
      auto background = makeBackground(resolution);
      auto cut = cutShape(background, shape, Real(0.25), Real(0.05), Real(0.2));
      auto optimized = cut.mesh;
      const bool tmopOk =
        solveSquaredTMOP(optimized, shape, Real(0.25), qualityWeight);
      benchmark::DoNotOptimize(optimized.getCellCount());
      counters =
        setCounters(st, cut, optimized, nullptr, shape, tmopOk, Real(0.25));
    }
    st.counters["quality_weight"] = benchmark::Counter(qualityWeight);
    st.SetItemsProcessed(st.iterations() * counters.cutCells);
  }

  static void BM_LevelSetPipeline_CarryForwardOrbit(benchmark::State& st)
  {
    const auto resolution = static_cast<size_t>(st.range(0));
    const auto steps = static_cast<Index>(st.range(1));
    const auto shape = static_cast<ShapeCase>(st.range(2));
    const auto metric = static_cast<MetricCase>(st.range(3));
    StageCounters counters;
    for (auto _ : st)
    {
      auto background = makeBackground(resolution);
      Index tmopFailures = 0;
      for (Index step = 0; step < steps; ++step)
      {
        const Real t = static_cast<Real>(step)
          / static_cast<Real>(std::max<Index>(1, steps - 1));
        auto cut = cutShape(background, shape, t, Real(0.05), Real(0.2));
        auto optimized = cut.mesh;
        const bool tmopOk = solveTMOP(optimized, shape, t, metric);
        if (!tmopOk)
          tmopFailures++;
        const auto report = coarsen(optimized, resolution);
        optimized.getConnectivity().compute(2, 1);
        optimized.getConnectivity().compute(1, 0);
        counters = setCounters(st, cut, optimized, &report, shape, tmopOk, t);
        demoteInterface(optimized);
        background = std::move(optimized);
      }
      counters.tmopFailures = tmopFailures;
      benchmark::DoNotOptimize(background.getCellCount());
    }
    st.counters["orbit_steps"] = benchmark::Counter(steps);
    st.counters["tmop_failed"] = benchmark::Counter(counters.tmopFailures);
    st.SetItemsProcessed(st.iterations() * counters.finalCells);
  }

  static void BM_LevelSetPipeline_CarryForwardSquaredWeight(
      benchmark::State& st)
  {
    const auto resolution = static_cast<size_t>(st.range(0));
    const auto steps = static_cast<Index>(st.range(1));
    const auto shape = static_cast<ShapeCase>(st.range(2));
    const Real qualityWeight =
      static_cast<Real>(st.range(3)) / static_cast<Real>(1000);
    StageCounters counters;
    for (auto _ : st)
    {
      auto background = makeBackground(resolution);
      Index tmopFailures = 0;
      for (Index step = 0; step < steps; ++step)
      {
        const Real t = static_cast<Real>(step)
          / static_cast<Real>(std::max<Index>(1, steps - 1));
        auto cut = cutShape(background, shape, t, Real(0.05), Real(0.2));
        auto optimized = cut.mesh;
        const bool tmopOk =
          solveSquaredTMOP(optimized, shape, t, qualityWeight);
        if (!tmopOk)
          tmopFailures++;
        const auto report = coarsen(optimized, resolution);
        optimized.getConnectivity().compute(2, 1);
        optimized.getConnectivity().compute(1, 0);
        counters = setCounters(st, cut, optimized, &report, shape, tmopOk, t);
        demoteInterface(optimized);
        background = std::move(optimized);
      }
      counters.tmopFailures = tmopFailures;
      benchmark::DoNotOptimize(background.getCellCount());
    }
    st.counters["orbit_steps"] = benchmark::Counter(steps);
    st.counters["quality_weight"] = benchmark::Counter(qualityWeight);
    st.counters["tmop_failed"] = benchmark::Counter(counters.tmopFailures);
    st.SetItemsProcessed(st.iterations() * counters.finalCells);
  }

#define RODIN_CUT_EXACT_BENCH(SHAPE_LABEL, SHAPE_ENUM) \
  BENCHMARK(BM_LevelSetPipeline_CutExact) \
    ->Name("LevelSetPipeline/CutExact/" SHAPE_LABEL "/24") \
    ->Args({24, static_cast<int>(ShapeCase::SHAPE_ENUM)}) \
    ->Unit(benchmark::kMillisecond);

#define RODIN_CUT_QUALITY_BENCH(SHAPE_LABEL, SHAPE_ENUM, RESOLUTION) \
  BENCHMARK(BM_LevelSetPipeline_CutQualityAware) \
    ->Name("LevelSetPipeline/CutQualityAware/" SHAPE_LABEL "/" #RESOLUTION) \
    ->Args({RESOLUTION, static_cast<int>(ShapeCase::SHAPE_ENUM)}) \
    ->Unit(benchmark::kMillisecond);

#define RODIN_TMOP_METRIC_BENCH(SHAPE_LABEL, SHAPE_ENUM, METRIC_LABEL, METRIC_ENUM) \
  BENCHMARK(BM_LevelSetPipeline_TMOPMetricFreshCut) \
    ->Name("LevelSetPipeline/TMOPFresh/" SHAPE_LABEL "/" METRIC_LABEL "/24") \
    ->Args({24, static_cast<int>(ShapeCase::SHAPE_ENUM), \
      static_cast<int>(MetricCase::METRIC_ENUM)}) \
    ->Unit(benchmark::kMillisecond);

#define RODIN_COARSEN_BENCH(SHAPE_LABEL, SHAPE_ENUM) \
  BENCHMARK(BM_LevelSetPipeline_CoarsenAfterTMOP) \
    ->Name("LevelSetPipeline/CoarsenAfterTMOP/" SHAPE_LABEL "/24") \
    ->Args({24, static_cast<int>(ShapeCase::SHAPE_ENUM)}) \
    ->Unit(benchmark::kMillisecond);

#define RODIN_ORBIT_BENCH(SHAPE_LABEL, SHAPE_ENUM, METRIC_LABEL, METRIC_ENUM) \
  BENCHMARK(BM_LevelSetPipeline_CarryForwardOrbit) \
    ->Name("LevelSetPipeline/CarryForward/" SHAPE_LABEL "/" METRIC_LABEL "/24x8") \
    ->Args({24, 8, static_cast<int>(ShapeCase::SHAPE_ENUM), \
      static_cast<int>(MetricCase::METRIC_ENUM)}) \
    ->Unit(benchmark::kMillisecond);

#define RODIN_WEIGHT_FRESH_BENCH(SHAPE_LABEL, SHAPE_ENUM, WEIGHT_MILLI) \
  BENCHMARK(BM_LevelSetPipeline_TMOPSquaredWeightFreshCut) \
    ->Name("LevelSetPipeline/TMOPWeightFresh/" SHAPE_LABEL "/SquaredDistance/w" #WEIGHT_MILLI "/24") \
    ->Args({24, static_cast<int>(ShapeCase::SHAPE_ENUM), WEIGHT_MILLI}) \
    ->Unit(benchmark::kMillisecond);

#define RODIN_WEIGHT_ORBIT_BENCH(SHAPE_LABEL, SHAPE_ENUM, WEIGHT_MILLI) \
  BENCHMARK(BM_LevelSetPipeline_CarryForwardSquaredWeight) \
    ->Name("LevelSetPipeline/CarryForwardWeight/" SHAPE_LABEL "/SquaredDistance/w" #WEIGHT_MILLI "/24x8") \
    ->Args({24, 8, static_cast<int>(ShapeCase::SHAPE_ENUM), WEIGHT_MILLI}) \
    ->Unit(benchmark::kMillisecond);

  RODIN_CUT_EXACT_BENCH("Circle", CircleOrbit)
  RODIN_CUT_EXACT_BENCH("Square", SquareOrbit)
  RODIN_CUT_EXACT_BENCH("Triangle", TriangleOrbit)
  RODIN_CUT_EXACT_BENCH("CircleEnterLeave", CircleEnterLeave)

  RODIN_CUT_QUALITY_BENCH("Circle", CircleOrbit, 24)
  RODIN_CUT_QUALITY_BENCH("Circle", CircleOrbit, 48)
  RODIN_CUT_QUALITY_BENCH("Square", SquareOrbit, 24)
  RODIN_CUT_QUALITY_BENCH("Triangle", TriangleOrbit, 24)
  RODIN_CUT_QUALITY_BENCH("CircleEnterLeave", CircleEnterLeave, 24)

  RODIN_TMOP_METRIC_BENCH("Circle", CircleOrbit, "SquaredDistance", SquaredDistance)
  RODIN_TMOP_METRIC_BENCH("Circle", CircleOrbit, "AreaDistortion", AreaDistortion)
  RODIN_TMOP_METRIC_BENCH("Circle", CircleOrbit, "ShapeSize", ShapeSize)
  RODIN_TMOP_METRIC_BENCH("Circle", CircleOrbit, "ShapeDistortion", ShapeDistortion)
  RODIN_TMOP_METRIC_BENCH("Square", SquareOrbit, "SquaredDistance", SquaredDistance)
  RODIN_TMOP_METRIC_BENCH("Square", SquareOrbit, "AreaDistortion", AreaDistortion)
  RODIN_TMOP_METRIC_BENCH("Square", SquareOrbit, "ShapeSize", ShapeSize)
  RODIN_TMOP_METRIC_BENCH("Square", SquareOrbit, "ShapeDistortion", ShapeDistortion)
  RODIN_TMOP_METRIC_BENCH("Triangle", TriangleOrbit, "SquaredDistance", SquaredDistance)
  RODIN_TMOP_METRIC_BENCH("Triangle", TriangleOrbit, "AreaDistortion", AreaDistortion)
  RODIN_TMOP_METRIC_BENCH("Triangle", TriangleOrbit, "ShapeSize", ShapeSize)
  RODIN_TMOP_METRIC_BENCH("Triangle", TriangleOrbit, "ShapeDistortion", ShapeDistortion)
  RODIN_TMOP_METRIC_BENCH("CircleEnterLeave", CircleEnterLeave, "SquaredDistance", SquaredDistance)
  RODIN_TMOP_METRIC_BENCH("CircleEnterLeave", CircleEnterLeave, "AreaDistortion", AreaDistortion)
  RODIN_TMOP_METRIC_BENCH("CircleEnterLeave", CircleEnterLeave, "ShapeSize", ShapeSize)
  RODIN_TMOP_METRIC_BENCH("CircleEnterLeave", CircleEnterLeave, "ShapeDistortion", ShapeDistortion)

  RODIN_WEIGHT_FRESH_BENCH("Circle", CircleOrbit, 50)
  RODIN_WEIGHT_FRESH_BENCH("Circle", CircleOrbit, 100)
  RODIN_WEIGHT_FRESH_BENCH("Circle", CircleOrbit, 200)
  RODIN_WEIGHT_FRESH_BENCH("Circle", CircleOrbit, 500)

  RODIN_COARSEN_BENCH("Circle", CircleOrbit)
  RODIN_COARSEN_BENCH("Square", SquareOrbit)
  RODIN_COARSEN_BENCH("Triangle", TriangleOrbit)
  RODIN_COARSEN_BENCH("CircleEnterLeave", CircleEnterLeave)

  RODIN_ORBIT_BENCH("Circle", CircleOrbit, "SquaredDistance", SquaredDistance)
  RODIN_ORBIT_BENCH("Circle", CircleOrbit, "AreaDistortion", AreaDistortion)
  RODIN_ORBIT_BENCH("Circle", CircleOrbit, "ShapeSize", ShapeSize)
  RODIN_ORBIT_BENCH("Circle", CircleOrbit, "ShapeDistortion", ShapeDistortion)
  RODIN_ORBIT_BENCH("Square", SquareOrbit, "SquaredDistance", SquaredDistance)
  RODIN_ORBIT_BENCH("Triangle", TriangleOrbit, "SquaredDistance", SquaredDistance)
  RODIN_ORBIT_BENCH("CircleEnterLeave", CircleEnterLeave, "SquaredDistance", SquaredDistance)

  RODIN_WEIGHT_ORBIT_BENCH("Circle", CircleOrbit, 50)
  RODIN_WEIGHT_ORBIT_BENCH("Circle", CircleOrbit, 100)
  RODIN_WEIGHT_ORBIT_BENCH("Circle", CircleOrbit, 200)

#undef RODIN_CUT_EXACT_BENCH
#undef RODIN_CUT_QUALITY_BENCH
#undef RODIN_TMOP_METRIC_BENCH
#undef RODIN_COARSEN_BENCH
#undef RODIN_ORBIT_BENCH
#undef RODIN_WEIGHT_FRESH_BENCH
#undef RODIN_WEIGHT_ORBIT_BENCH
}
