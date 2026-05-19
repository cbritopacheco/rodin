/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <array>
#include <chrono>
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

#include <type_traits>
#include <Eigen/SparseLU>

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
      CircleEnterLeave = 3,
      CircleFigureEight = 4,
      SquareFigureEight = 5,
      TriangleFigureEight = 6,
      WavyCircleOrbit = 7,
      WavyCircleFigureEight = 8,
      WavyCircleEnterLeave = 9
    };

    enum class MetricCase : int
    {
      SquaredDistance = 0,
      AreaDistortion = 1,
      ShapeSize = 2,
      ShapeDistortion = 3
    };

    enum class StageCase : int
    {
      TMOPOnly = 0,
      OptimizeOnly = 1,
      TMOPThenOptimize = 2,     // diagnostic: topology changes after TMOP
      OptimizeThenTMOP = 3,     // production: topology first, TMOP last
      OptimizeNoFeature = 4,
      OptimizeNoInterior = 5
    };

    constexpr StageCase ProductionStageCase = StageCase::OptimizeThenTMOP;

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
      Index smooths = 0;
      Index featureSmooths = 0;
      Real cutMinQuality = 0;
      Real finalMinQuality = 0;
      Real fitRms = 0;
      Real fitMax = 0;
      Real curvedFitRms = 0;
      Real curvedFitMax = 0;
      Real curvedMinQuality = 0;
      Real curvedMinJacobian = 0;
      Index curvedInvalidSamples = 0;
      bool hasP2Diagnostics = false;
      Real coverage = 0;
      Real signedArea = 0;
      Real maxInterfaceDeviation = 0;
      Real cutSeconds = 0;
      Real optimizeSeconds = 0;
      Real tmopSeconds = 0;
      Real tmopAssemblySeconds = 0;
      Real tmopSolveSeconds = 0;
      Real tmopMeritSeconds = 0;
      Index tmopIterations = 0;
      Index inverted = 0;
      Index tmopFailures = 0;
    };

    struct StageAverages
    {
      Index samples = 0;
      Real cutCells = 0;
      Real finalCells = 0;
      Real interfaceEdges = 0;
      Real uncutCells = 0;
      Real snappedCrossings = 0;
      Real pathologicalCuts = 0;
      Real cutMinQuality = 0;
      Real finalMinQuality = 0;
      Real minFinalQuality = std::numeric_limits<Real>::infinity();
      Real maxFinalQuality = -std::numeric_limits<Real>::infinity();
      Real finalStepQuality = 0;
      Real fitRms = 0;
      Real fitMax = 0;
      Real bestFitRms = std::numeric_limits<Real>::infinity();
      Real worstFitRms = 0;
      Real finalFitRms = 0;
      Real bestFitMax = std::numeric_limits<Real>::infinity();
      Real worstFitMax = 0;
      Real finalFitMax = 0;
      Real curvedFitRms = 0;
      Real curvedFitMax = 0;
      Real curvedMinQuality = 0;
      Real minCurvedQuality = std::numeric_limits<Real>::infinity();
      Real maxCurvedQuality = -std::numeric_limits<Real>::infinity();
      Real finalCurvedQuality = 0;
      Real curvedMinJacobian = std::numeric_limits<Real>::infinity();
      Real bestCurvedFitRms = std::numeric_limits<Real>::infinity();
      Real worstCurvedFitRms = 0;
      Real finalCurvedFitRms = 0;
      Real bestCurvedFitMax = std::numeric_limits<Real>::infinity();
      Real worstCurvedFitMax = 0;
      Real finalCurvedFitMax = 0;
      Real bestCurvedMinJacobian = -std::numeric_limits<Real>::infinity();
      Real finalCurvedMinJacobian = 0;
      Real curvedInvalidSamples = 0;
      Real bestCurvedInvalidSamples = std::numeric_limits<Real>::infinity();
      Real worstCurvedInvalidSamples = 0;
      Real finalCurvedInvalidSamples = 0;
      Real coverage = 0;
      Real signedArea = 0;
      Real finalCoverage = 0;
      Real bestCoverage = -std::numeric_limits<Real>::infinity();
      Real worstCoverage = std::numeric_limits<Real>::infinity();
      Real maxInterfaceDeviation = 0;
      Real cutSeconds = 0;
      Real optimizeSeconds = 0;
      Real tmopSeconds = 0;
      Real tmopAssemblySeconds = 0;
      Real tmopSolveSeconds = 0;
      Real tmopMeritSeconds = 0;
      Real tmopIterations = 0;
      Real inverted = 0;
      Real tmopFailures = 0;
      bool hasP2Diagnostics = false;

      void add(const StageCounters& counters)
      {
        ++samples;
        cutCells += static_cast<Real>(counters.cutCells);
        finalCells += static_cast<Real>(counters.finalCells);
        interfaceEdges += static_cast<Real>(counters.interfaceEdges);
        uncutCells += static_cast<Real>(counters.uncutCells);
        snappedCrossings += static_cast<Real>(counters.snappedCrossings);
        pathologicalCuts += static_cast<Real>(counters.pathologicalCuts);
        cutMinQuality += counters.cutMinQuality;
        finalMinQuality += counters.finalMinQuality;
        minFinalQuality =
          std::min(minFinalQuality, counters.finalMinQuality);
        maxFinalQuality =
          std::max(maxFinalQuality, counters.finalMinQuality);
        finalStepQuality = counters.finalMinQuality;
        fitRms += counters.fitRms;
        fitMax += counters.fitMax;
        bestFitRms = std::min(bestFitRms, counters.fitRms);
        worstFitRms = std::max(worstFitRms, counters.fitRms);
        finalFitRms = counters.fitRms;
        bestFitMax = std::min(bestFitMax, counters.fitMax);
        worstFitMax = std::max(worstFitMax, counters.fitMax);
        finalFitMax = counters.fitMax;
        curvedFitRms += counters.curvedFitRms;
        curvedFitMax += counters.curvedFitMax;
        curvedMinQuality += counters.curvedMinQuality;
        minCurvedQuality =
          std::min(minCurvedQuality, counters.curvedMinQuality);
        maxCurvedQuality =
          std::max(maxCurvedQuality, counters.curvedMinQuality);
        finalCurvedQuality = counters.curvedMinQuality;
        bestCurvedFitRms =
          std::min(bestCurvedFitRms, counters.curvedFitRms);
        worstCurvedFitRms =
          std::max(worstCurvedFitRms, counters.curvedFitRms);
        finalCurvedFitRms = counters.curvedFitRms;
        bestCurvedFitMax =
          std::min(bestCurvedFitMax, counters.curvedFitMax);
        worstCurvedFitMax =
          std::max(worstCurvedFitMax, counters.curvedFitMax);
        finalCurvedFitMax = counters.curvedFitMax;
        if (std::isfinite(counters.curvedMinJacobian))
        {
          curvedMinJacobian =
            std::min(curvedMinJacobian, counters.curvedMinJacobian);
          bestCurvedMinJacobian =
            std::max(bestCurvedMinJacobian, counters.curvedMinJacobian);
        }
        finalCurvedMinJacobian = counters.curvedMinJacobian;
        curvedInvalidSamples +=
          static_cast<Real>(counters.curvedInvalidSamples);
        bestCurvedInvalidSamples = std::min(
            bestCurvedInvalidSamples,
            static_cast<Real>(counters.curvedInvalidSamples));
        worstCurvedInvalidSamples = std::max(
            worstCurvedInvalidSamples,
            static_cast<Real>(counters.curvedInvalidSamples));
        finalCurvedInvalidSamples =
          static_cast<Real>(counters.curvedInvalidSamples);
        hasP2Diagnostics = hasP2Diagnostics || counters.hasP2Diagnostics;
        coverage += counters.coverage;
        bestCoverage = std::max(bestCoverage, counters.coverage);
        worstCoverage = std::min(worstCoverage, counters.coverage);
        finalCoverage = counters.coverage;
        signedArea += counters.signedArea;
        maxInterfaceDeviation += counters.maxInterfaceDeviation;
        cutSeconds += counters.cutSeconds;
        optimizeSeconds += counters.optimizeSeconds;
        tmopSeconds += counters.tmopSeconds;
        tmopAssemblySeconds += counters.tmopAssemblySeconds;
        tmopSolveSeconds += counters.tmopSolveSeconds;
        tmopMeritSeconds += counters.tmopMeritSeconds;
        tmopIterations += static_cast<Real>(counters.tmopIterations);
        inverted += static_cast<Real>(counters.inverted);
        tmopFailures += static_cast<Real>(counters.tmopFailures);
      }

      void publish(benchmark::State& st) const
      {
        if (samples == 0)
          return;
        const Real inv = Real(1) / static_cast<Real>(samples);
        st.counters["avg_cut_cells"] = benchmark::Counter(cutCells * inv);
        st.counters["avg_final_cells"] =
          benchmark::Counter(finalCells * inv);
        st.counters["avg_interface_edges"] =
          benchmark::Counter(interfaceEdges * inv);
        st.counters["avg_uncut_cells"] =
          benchmark::Counter(uncutCells * inv);
        st.counters["avg_snapped"] =
          benchmark::Counter(snappedCrossings * inv);
        st.counters["avg_pathological"] =
          benchmark::Counter(pathologicalCuts * inv);
        st.counters["avg_qmin_cut"] =
          benchmark::Counter(cutMinQuality * inv);
        st.counters["avg_qmin_final"] =
          benchmark::Counter(finalMinQuality * inv);
        st.counters["final_qmin_final"] =
          benchmark::Counter(finalStepQuality);
        st.counters["best_qmin_final"] =
          benchmark::Counter(std::isfinite(maxFinalQuality) ? maxFinalQuality : Real(0));
        st.counters["worst_qmin_final"] =
          benchmark::Counter(std::isfinite(minFinalQuality) ? minFinalQuality : Real(0));
        st.counters["min_qmin_final"] =
          benchmark::Counter(minFinalQuality);
        st.counters["avg_fit_rms"] = benchmark::Counter(fitRms * inv);
        st.counters["final_fit_rms"] = benchmark::Counter(finalFitRms);
        st.counters["best_fit_rms"] =
          benchmark::Counter(std::isfinite(bestFitRms) ? bestFitRms : Real(0));
        st.counters["worst_fit_rms"] = benchmark::Counter(worstFitRms);
        st.counters["avg_fit_max"] = benchmark::Counter(fitMax * inv);
        st.counters["final_fit_max"] = benchmark::Counter(finalFitMax);
        st.counters["best_fit_max"] =
          benchmark::Counter(std::isfinite(bestFitMax) ? bestFitMax : Real(0));
        st.counters["worst_fit_max"] = benchmark::Counter(worstFitMax);
        st.counters["avg_curved_fit_rms"] =
          benchmark::Counter(curvedFitRms * inv);
        st.counters["avg_curved_qmin"] =
          benchmark::Counter(curvedMinQuality * inv);
        st.counters["final_curved_qmin"] =
          benchmark::Counter(finalCurvedQuality);
        st.counters["best_curved_qmin"] =
          benchmark::Counter(std::isfinite(maxCurvedQuality) ? maxCurvedQuality : Real(0));
        st.counters["worst_curved_qmin"] =
          benchmark::Counter(std::isfinite(minCurvedQuality) ? minCurvedQuality : Real(0));
        st.counters["final_curved_fit_rms"] =
          benchmark::Counter(finalCurvedFitRms);
        st.counters["best_curved_fit_rms"] =
          benchmark::Counter(std::isfinite(bestCurvedFitRms) ? bestCurvedFitRms : Real(0));
        st.counters["worst_curved_fit_rms"] =
          benchmark::Counter(worstCurvedFitRms);
        st.counters["avg_curved_fit_max"] =
          benchmark::Counter(curvedFitMax * inv);
        st.counters["final_curved_fit_max"] =
          benchmark::Counter(finalCurvedFitMax);
        st.counters["best_curved_fit_max"] =
          benchmark::Counter(std::isfinite(bestCurvedFitMax) ? bestCurvedFitMax : Real(0));
        st.counters["worst_curved_fit_max"] =
          benchmark::Counter(worstCurvedFitMax);
        st.counters["min_curved_jac"] =
          benchmark::Counter(std::isfinite(curvedMinJacobian) ? curvedMinJacobian : Real(0));
        st.counters["final_curved_min_jac"] =
          benchmark::Counter(finalCurvedMinJacobian);
        st.counters["best_curved_min_jac"] =
          benchmark::Counter(std::isfinite(bestCurvedMinJacobian) ? bestCurvedMinJacobian : Real(0));
        st.counters["worst_curved_min_jac"] =
          benchmark::Counter(std::isfinite(curvedMinJacobian) ? curvedMinJacobian : Real(0));
        st.counters["avg_curved_invalid"] =
          benchmark::Counter(curvedInvalidSamples * inv);
        st.counters["final_curved_invalid"] =
          benchmark::Counter(finalCurvedInvalidSamples);
        st.counters["best_curved_invalid"] =
          benchmark::Counter(std::isfinite(bestCurvedInvalidSamples) ? bestCurvedInvalidSamples : Real(0));
        st.counters["worst_curved_invalid"] =
          benchmark::Counter(worstCurvedInvalidSamples);
        if (hasP2Diagnostics)
        {
          st.counters["avg_p2_fit_rms"] =
            benchmark::Counter(curvedFitRms * inv);
          st.counters["avg_p2_qmin"] =
            benchmark::Counter(curvedMinQuality * inv);
          st.counters["final_p2_qmin"] =
            benchmark::Counter(finalCurvedQuality);
          st.counters["best_p2_qmin"] =
            benchmark::Counter(std::isfinite(maxCurvedQuality) ? maxCurvedQuality : Real(0));
          st.counters["worst_p2_qmin"] =
            benchmark::Counter(std::isfinite(minCurvedQuality) ? minCurvedQuality : Real(0));
          st.counters["final_p2_fit_rms"] =
            benchmark::Counter(finalCurvedFitRms);
          st.counters["best_p2_fit_rms"] =
            benchmark::Counter(std::isfinite(bestCurvedFitRms) ? bestCurvedFitRms : Real(0));
          st.counters["worst_p2_fit_rms"] =
            benchmark::Counter(worstCurvedFitRms);
          st.counters["avg_p2_fit_max"] =
            benchmark::Counter(curvedFitMax * inv);
          st.counters["final_p2_fit_max"] =
            benchmark::Counter(finalCurvedFitMax);
          st.counters["best_p2_fit_max"] =
            benchmark::Counter(std::isfinite(bestCurvedFitMax) ? bestCurvedFitMax : Real(0));
          st.counters["worst_p2_fit_max"] =
            benchmark::Counter(worstCurvedFitMax);
          st.counters["final_p2_min_jac"] =
            benchmark::Counter(finalCurvedMinJacobian);
          st.counters["best_p2_min_jac"] =
            benchmark::Counter(std::isfinite(bestCurvedMinJacobian) ? bestCurvedMinJacobian : Real(0));
          st.counters["worst_p2_min_jac"] =
            benchmark::Counter(std::isfinite(curvedMinJacobian) ? curvedMinJacobian : Real(0));
          st.counters["final_p2_invalid"] =
            benchmark::Counter(finalCurvedInvalidSamples);
          st.counters["best_p2_invalid"] =
            benchmark::Counter(std::isfinite(bestCurvedInvalidSamples) ? bestCurvedInvalidSamples : Real(0));
          st.counters["worst_p2_invalid"] =
            benchmark::Counter(worstCurvedInvalidSamples);
        }
        st.counters["avg_coverage"] = benchmark::Counter(coverage * inv);
        st.counters["final_coverage"] = benchmark::Counter(finalCoverage);
        st.counters["best_coverage"] =
          benchmark::Counter(std::isfinite(bestCoverage) ? bestCoverage : Real(0));
        st.counters["worst_coverage"] =
          benchmark::Counter(std::isfinite(worstCoverage) ? worstCoverage : Real(0));
        st.counters["avg_signed_area"] =
          benchmark::Counter(signedArea * inv);
        st.counters["avg_max_iface_dev"] =
          benchmark::Counter(maxInterfaceDeviation * inv);
        st.counters["avg_cut_ms"] =
          benchmark::Counter(Real(1000) * cutSeconds * inv);
        st.counters["avg_optimize_ms"] =
          benchmark::Counter(Real(1000) * optimizeSeconds * inv);
        st.counters["avg_tmop_ms"] =
          benchmark::Counter(Real(1000) * tmopSeconds * inv);
        st.counters["avg_tmop_assembly_ms"] =
          benchmark::Counter(Real(1000) * tmopAssemblySeconds * inv);
        st.counters["avg_tmop_solve_ms"] =
          benchmark::Counter(Real(1000) * tmopSolveSeconds * inv);
        st.counters["avg_tmop_merit_ms"] =
          benchmark::Counter(Real(1000) * tmopMeritSeconds * inv);
        st.counters["avg_tmop_iterations"] =
          benchmark::Counter(tmopIterations * inv);
        st.counters["avg_inverted"] = benchmark::Counter(inverted * inv);
        st.counters["avg_tmop_failed"] =
          benchmark::Counter(tmopFailures * inv);
      }
    };

    using BenchClock = std::chrono::steady_clock;

    Real elapsedSeconds(BenchClock::time_point start)
    {
      return std::chrono::duration<Real>(BenchClock::now() - start).count();
    }

    struct TMOPSolveStats
    {
      Real assemblySeconds = 0;
      Real solveSeconds = 0;
      Real meritSeconds = 0;
      Index iterations = 0;
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

    Math::SpatialPoint figureEightCenter(Real t)
    {
      return Math::SpatialPoint{
        Real(0.5) + Real(0.22) * std::sin(Real(2) * Pi * t),
        Real(0.5) + Real(0.15) * std::sin(Real(4) * Pi * t + Real(0.35))
      };
    }

    Math::SpatialPoint shapeCenter(ShapeCase shape, Real t)
    {
      switch (shape)
      {
        case ShapeCase::CircleOrbit:
        case ShapeCase::SquareOrbit:
        case ShapeCase::TriangleOrbit:
        case ShapeCase::WavyCircleOrbit:
          return orbitCenter(t);
        case ShapeCase::CircleEnterLeave:
        case ShapeCase::WavyCircleEnterLeave:
          return enterLeaveCenter(t);
        case ShapeCase::CircleFigureEight:
        case ShapeCase::SquareFigureEight:
        case ShapeCase::TriangleFigureEight:
        case ShapeCase::WavyCircleFigureEight:
          return figureEightCenter(t);
      }
      return orbitCenter(t);
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

    Real squareValue(
        const Math::SpatialPoint& x,
        const Math::SpatialPoint& c)
    {
      const Real half = Real(0.16);
      const Real ax = std::abs(x[0] - c[0]) - half;
      const Real ay = std::abs(x[1] - c[1]) - half;
      const Real ox = std::max(ax, Real(0));
      const Real oy = std::max(ay, Real(0));
      const Real outside = std::hypot(ox, oy);
      const Real inside = std::min(std::max(ax, ay), Real(0));
      return outside + inside;
    }

    Math::SpatialPoint squareGradient(
        const Math::SpatialPoint& x,
        const Math::SpatialPoint& c)
    {
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

    std::array<Math::SpatialPoint, 3> triangleVertices(
        const Math::SpatialPoint& c)
    {
      return {{
        Math::SpatialPoint{c[0], c[1] + Real(0.20)},
        Math::SpatialPoint{c[0] - Real(0.18), c[1] - Real(0.12)},
        Math::SpatialPoint{c[0] + Real(0.18), c[1] - Real(0.12)}
      }};
    }

    Real triangleValue(
        const Math::SpatialPoint& x,
        const Math::SpatialPoint& c)
    {
      const auto vertices = triangleVertices(c);
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

    Math::SpatialPoint triangleGradient(
        const Math::SpatialPoint& x,
        const Math::SpatialPoint& c)
    {
      const auto vertices = triangleVertices(c);
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

    Real wavyCircleRadius(Real theta, Real t)
    {
      return Real(0.17)
        + Real(0.035) * std::sin(Real(5) * theta + Real(0.65) * std::sin(Real(2) * Pi * t));
    }

    Real wavyCircleRadiusDerivative(Real theta, Real t)
    {
      return Real(0.175) * std::cos(
          Real(5) * theta + Real(0.65) * std::sin(Real(2) * Pi * t));
    }

    Real wavyCircleValue(
        const Math::SpatialPoint& x,
        const Math::SpatialPoint& c,
        Real t)
    {
      const Real dx = x[0] - c[0];
      const Real dy = x[1] - c[1];
      const Real theta = std::atan2(dy, dx);
      return std::hypot(dx, dy) - wavyCircleRadius(theta, t);
    }

    Math::SpatialPoint wavyCircleGradient(
        const Math::SpatialPoint& x,
        const Math::SpatialPoint& c,
        Real t)
    {
      const Real dx = x[0] - c[0];
      const Real dy = x[1] - c[1];
      const Real rho2 = dx * dx + dy * dy;
      const Real rho = std::sqrt(rho2);
      if (rho <= Real(1e-14))
        return Math::SpatialPoint{1, 0};
      const Real theta = std::atan2(dy, dx);
      const Real dr = wavyCircleRadiusDerivative(theta, t);
      return Math::SpatialPoint{
        dx / rho + dr * dy / rho2,
        dy / rho - dr * dx / rho2
      };
    }

    Real levelSetValue(ShapeCase shape, const Math::SpatialPoint& x, Real t)
    {
      const auto c = shapeCenter(shape, t);
      switch (shape)
      {
        case ShapeCase::CircleOrbit:
        case ShapeCase::CircleFigureEight:
          return circleValue(x, c, Real(0.18));
        case ShapeCase::SquareOrbit:
        case ShapeCase::SquareFigureEight:
          return squareValue(x, c);
        case ShapeCase::TriangleOrbit:
        case ShapeCase::TriangleFigureEight:
          return triangleValue(x, c);
        case ShapeCase::CircleEnterLeave:
          return circleValue(x, c, Real(0.18));
        case ShapeCase::WavyCircleOrbit:
        case ShapeCase::WavyCircleFigureEight:
        case ShapeCase::WavyCircleEnterLeave:
          return wavyCircleValue(x, c, t);
      }
      return 1;
    }

    Math::SpatialPoint levelSetGradient(
        ShapeCase shape,
        const Math::SpatialPoint& x,
        Real t)
    {
      const auto c = shapeCenter(shape, t);
      switch (shape)
      {
        case ShapeCase::CircleOrbit:
        case ShapeCase::CircleFigureEight:
          return circleGradient(x, c);
        case ShapeCase::SquareOrbit:
        case ShapeCase::SquareFigureEight:
          return squareGradient(x, c);
        case ShapeCase::TriangleOrbit:
        case ShapeCase::TriangleFigureEight:
          return triangleGradient(x, c);
        case ShapeCase::CircleEnterLeave:
          return circleGradient(x, c);
        case ShapeCase::WavyCircleOrbit:
        case ShapeCase::WavyCircleFigureEight:
        case ShapeCase::WavyCircleEnterLeave:
          return wavyCircleGradient(x, c, t);
      }
      return Math::SpatialPoint{0, 0};
    }

    Math::SpatialPoint projectToLevelSet(
        ShapeCase shape,
        Real t,
        const Math::SpatialPoint& x,
        int iterations = 4)
    {
      Math::SpatialPoint p = x;
      for (int i = 0; i < iterations; ++i)
      {
        const Real f = levelSetValue(shape, p, t);
        const auto g = levelSetGradient(shape, p, t);
        const Real gg = g[0] * g[0] + g[1] * g[1];
        if (gg < Real(1e-30))
          break;
        p = Math::SpatialPoint{ p[0] - f * g[0] / gg,
                                p[1] - f * g[1] / gg };
      }
      return p;
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

    struct CurvedInterfaceStats
    {
      Real fitRms = 0;
      Real fitMax = 0;
      Real minQuality = std::numeric_limits<Real>::infinity();
      Real minJacobian = std::numeric_limits<Real>::infinity();
      Index invalidJacobianSamples = 0;
    };

    Math::SpatialPoint triangleReferenceVertex(std::uint8_t local)
    {
      if (local == 1)
        return Math::SpatialPoint{1, 0};
      if (local == 2)
        return Math::SpatialPoint{0, 1};
      return Math::SpatialPoint{0, 0};
    }

    Math::SpatialPoint triangleReferenceEdgePoint(
        std::uint8_t localA,
        std::uint8_t localB,
        Real s)
    {
      return (Real(1) - s) * triangleReferenceVertex(localA)
           + s * triangleReferenceVertex(localB);
    }

    bool referencePointOnTriangleEdge(
        const Math::SpatialPoint& rc,
        std::uint8_t localA,
        std::uint8_t localB)
    {
      const Real eps = Real(1e-12);
      if ((localA == 0 && localB == 1) || (localA == 1 && localB == 0))
        return std::abs(rc[1]) <= eps;
      if ((localA == 1 && localB == 2) || (localA == 2 && localB == 1))
        return std::abs(rc[0] + rc[1] - Real(1)) <= eps;
      return std::abs(rc[0]) <= eps;
    }

    bool referencePointIsTriangleVertex(const Math::SpatialPoint& rc)
    {
      const Real eps = Real(1e-12);
      return (std::abs(rc[0]) <= eps && std::abs(rc[1]) <= eps)
          || (std::abs(rc[0] - Real(1)) <= eps && std::abs(rc[1]) <= eps)
          || (std::abs(rc[0]) <= eps && std::abs(rc[1] - Real(1)) <= eps);
    }

    Optional<std::array<std::uint8_t, 2>> localEdgeForMeshEdge(
        const Geometry::Polytope::Key& cell,
        const Geometry::Polytope::Key& edge)
    {
      for (std::uint8_t a = 0; a < 3; ++a)
        for (std::uint8_t b = a + 1; b < 3; ++b)
        {
          const bool same =
            (cell(a) == edge(0) && cell(b) == edge(1))
         || (cell(a) == edge(1) && cell(b) == edge(0));
          if (same)
            return std::array<std::uint8_t, 2>{{a, b}};
        }
      return {};
    }

    std::vector<std::array<std::uint8_t, 2>> interfaceLocalEdges(
        const LocalMesh& mesh,
        Index cellIndex)
    {
      std::vector<std::array<std::uint8_t, 2>> edges;
      const auto& conn = mesh.getConnectivity();
      if (conn.getGeometry(2, cellIndex) != Polytope::Type::Triangle)
        return edges;

      const auto& cell = conn.getPolytope(2, cellIndex);
      for (Index edgeIndex : conn.getIncidence({2, 1}, cellIndex))
      {
        const auto attr = mesh.getAttribute(1, edgeIndex);
        if (!attr || *attr != Interface)
          continue;
        if (const auto edge = localEdgeForMeshEdge(
              cell, conn.getPolytope(1, edgeIndex)))
          edges.push_back(*edge);
      }
      return edges;
    }

    void installQuadraticInterfaceGeometry(
        LocalMesh& mesh,
        ShapeCase shape,
        Real t)
    {
      if (mesh.getConnectivity().getIncidence(2, 1).size() == 0)
        return;

      Variational::RealH1Element<2> fe(Polytope::Type::Triangle);
      Variational::RealH1Element<2> edgeFe(Polytope::Type::Segment);
      const auto& conn = mesh.getConnectivity();
      for (Index edgeIndex = 0;
           edgeIndex < static_cast<Index>(conn.getCount(1));
           ++edgeIndex)
      {
        const auto attr = mesh.getAttribute(1, edgeIndex);
        const bool isInterface = attr && *attr == Interface;

        const auto& edge = conn.getPolytope(1, edgeIndex);
        const auto x0 = mesh.getVertexCoordinates(edge(0));
        const auto x1 = mesh.getVertexCoordinates(edge(1));
        PointCloud pm(mesh.getSpaceDimension(), edgeFe.getCount());
        for (size_t local = 0; local < edgeFe.getCount(); ++local)
        {
          const Real s = edgeFe.getNode(local)[0];
          Math::SpatialPoint x = (Real(1) - s) * x0 + s * x1;
          if (isInterface
              && s > Real(1e-12)
              && s < Real(1) - Real(1e-12))
            x = projectToLevelSet(shape, t, x, 8);
          for (size_t d = 0; d < static_cast<size_t>(x.size()); ++d)
            pm(d, local) = x[d];
        }
        mesh.setPolytopeTransformation(
            {1, edgeIndex},
            new ParametricTransformation<Variational::RealH1Element<2>>(
                pm, Variational::RealH1Element<2>(edgeFe)));
      }

      for (Index cellIndex = 0;
           cellIndex < static_cast<Index>(mesh.getCellCount());
           ++cellIndex)
      {
        if (conn.getGeometry(2, cellIndex) != Polytope::Type::Triangle)
          continue;

        const auto featureEdges = interfaceLocalEdges(mesh, cellIndex);
        const auto& cell = conn.getPolytope(2, cellIndex);
        const auto x0 = mesh.getVertexCoordinates(cell(0));
        const auto x1 = mesh.getVertexCoordinates(cell(1));
        const auto x2 = mesh.getVertexCoordinates(cell(2));

        PointCloud pm(mesh.getSpaceDimension(), fe.getCount());
        for (size_t local = 0; local < fe.getCount(); ++local)
        {
          const auto& rc = fe.getNode(local);
          Math::SpatialPoint x =
            (Real(1) - rc[0] - rc[1]) * x0 + rc[0] * x1 + rc[1] * x2;

          const bool edgeInterior = !referencePointIsTriangleVertex(rc);
          bool onInterfaceEdge = false;
          for (const auto& e : featureEdges)
            onInterfaceEdge = onInterfaceEdge
              || referencePointOnTriangleEdge(rc, e[0], e[1]);
          if (edgeInterior && onInterfaceEdge)
            x = projectToLevelSet(shape, t, x, 6);

          for (size_t d = 0; d < static_cast<size_t>(x.size()); ++d)
            pm(d, local) = x[d];
        }

        auto* trans =
          new ParametricTransformation<Variational::RealH1Element<2>>(
              pm, Variational::RealH1Element<2>(fe));
        bool valid = true;
        for (const Math::SpatialPoint& rc : {
              Math::SpatialPoint{Real(1) / Real(3), Real(1) / Real(3)},
              Math::SpatialPoint{Real(0.5), Real(0.25)},
              Math::SpatialPoint{Real(0.25), Real(0.5)},
              Math::SpatialPoint{Real(0.25), Real(0.25)}})
        {
          Math::SpatialMatrix<Real> J;
          trans->jacobian(J, rc);
          valid = valid && J.rows() == 2 && J.cols() == 2
            && std::isfinite(J.determinant()) && J.determinant() > Real(1e-14);
        }
        if (valid)
          mesh.setPolytopeTransformation({2, cellIndex}, trans);
        else
          delete trans;
      }
    }

    CurvedInterfaceStats curvedInterfaceStats(
        const LocalMesh& mesh,
        ShapeCase shape,
        Real t)
    {
      CurvedInterfaceStats stats;
      const auto& conn = mesh.getConnectivity();
      Real sq = 0;
      Index fitSamples = 0;

      for (Index cellIndex = 0;
           cellIndex < static_cast<Index>(mesh.getCellCount());
           ++cellIndex)
      {
        if (conn.getGeometry(2, cellIndex) != Polytope::Type::Triangle)
          continue;

        const auto cellIt = mesh.getPolytope(2, cellIndex);
        const auto& trans = cellIt->getTransformation();
        for (const Math::SpatialPoint& rc : {
              Math::SpatialPoint{Real(1) / Real(3), Real(1) / Real(3)},
              Math::SpatialPoint{Real(0.5), Real(0.25)},
              Math::SpatialPoint{Real(0.25), Real(0.5)},
              Math::SpatialPoint{Real(0.25), Real(0.25)}})
        {
          Math::SpatialMatrix<Real> J;
          trans.jacobian(J, rc);
          if (J.rows() == 2 && J.cols() == 2)
          {
            const Real det = J.determinant();
            stats.minJacobian = std::min(stats.minJacobian, det);
            if (!(det > Real(0)) || !std::isfinite(det))
              stats.invalidJacobianSamples++;
          }
        }

        for (const auto& e : interfaceLocalEdges(mesh, cellIndex))
        {
          for (Real s : {Real(0.25), Real(0.5), Real(0.75)})
          {
            Math::SpatialPoint x;
            trans.transform(x, triangleReferenceEdgePoint(e[0], e[1], s));
            const Real phi = std::abs(levelSetValue(shape, x, t));
            stats.fitMax = std::max(stats.fitMax, phi);
            sq += phi * phi;
            fitSamples++;
          }
        }
      }

      if (fitSamples > 0)
        stats.fitRms = std::sqrt(sq / static_cast<Real>(fitSamples));
      for (Index cellIndex = 0;
           cellIndex < static_cast<Index>(mesh.getCellCount());
           ++cellIndex)
      {
        if (conn.getGeometry(2, cellIndex) != Polytope::Type::Triangle)
          continue;
        const auto cellIt = mesh.getPolytope(2, cellIndex);
        const auto& trans = cellIt->getTransformation();
        for (Index i = 0; i < 3; ++i)
        {
          for (Index j = 0; j < 3 - i; ++j)
          {
            auto sample = [&](Index a, Index b)
            {
              Math::SpatialPoint x;
              trans.transform(
                  x,
                  Math::SpatialPoint{
                    static_cast<Real>(a) / Real(3),
                    static_cast<Real>(b) / Real(3)});
              return x;
            };
            const auto x00 = sample(i, j);
            const auto x10 = sample(i + 1, j);
            const auto x01 = sample(i, j + 1);
            stats.minQuality =
              std::min(stats.minQuality, triangleQuality(x00, x10, x01));
            if (i + j + 1 < 3)
            {
              const auto x11 = sample(i + 1, j + 1);
              stats.minQuality =
                std::min(stats.minQuality, triangleQuality(x10, x11, x01));
            }
          }
        }
      }
      if (!std::isfinite(stats.minQuality))
        stats.minQuality = 0;
      if (!std::isfinite(stats.minJacobian))
        stats.minJacobian = 0;
      return stats;
    }

    template <class FES, class Data>
    void applyParametricDisplacement(
        LocalMesh& mesh,
        const Variational::GridFunction<FES, Data>& u)
    {
      const auto& fes = u.getFiniteElementSpace();
      const auto& data = u.getData();
      const size_t sdim = mesh.getSpaceDimension();
      if (sdim != 2 || fes.getVectorDimension() != 2)
        return;

      const auto& conn = mesh.getConnectivity();
      std::vector<PointCloud> edgePointClouds(conn.getCount(1));
      std::vector<PointCloud> cellPointClouds(mesh.getCellCount());
      std::vector<char> hasCellPointCloud(mesh.getCellCount(), 0);
      for (Index edgeIndex = 0;
           edgeIndex < static_cast<Index>(conn.getCount(1));
           ++edgeIndex)
      {
        const auto edgeIterator = mesh.getPolytope(1, edgeIndex);
        const auto& edge = *edgeIterator;
        const auto& fe =
          fes.getFiniteElement(edge.getDimension(), edge.getIndex());
        const size_t localSize = fe.getCount();
        PointCloud pm(sdim, localSize / sdim);
        for (size_t node = 0; node < localSize / sdim; ++node)
        {
          const auto& rc = fe.getNode(node * sdim);
          const Geometry::Point point(edge, rc);
          auto x = point.getPhysicalCoordinates();
          for (size_t c = 0; c < sdim; ++c)
          {
            const size_t local = node * sdim + c;
            const Index global = fes.getGlobalIndex(
                {edge.getDimension(), edge.getIndex()},
                static_cast<Index>(local));
            x[c] += data(global);
          }
          for (size_t c = 0; c < sdim; ++c)
            pm(c, node) = x[c];
        }
        edgePointClouds[static_cast<size_t>(edgeIndex)] = std::move(pm);
      }

      for (Index cellIndex = 0;
           cellIndex < static_cast<Index>(mesh.getCellCount());
           ++cellIndex)
      {
        if (conn.getGeometry(2, cellIndex) != Polytope::Type::Triangle)
          continue;
        const auto cellIterator = mesh.getPolytope(2, cellIndex);
        const auto& cell = *cellIterator;
        const auto& fe =
          fes.getFiniteElement(cell.getDimension(), cell.getIndex());
        const size_t localSize = fe.getCount();
        PointCloud pm(sdim, localSize / sdim);
        for (size_t node = 0; node < localSize / sdim; ++node)
        {
          const auto& rc = fe.getNode(node * sdim);
          const Geometry::Point point(cell, rc);
          auto x = point.getPhysicalCoordinates();
          for (size_t c = 0; c < sdim; ++c)
          {
            const size_t local = node * sdim + c;
            const Index global = fes.getGlobalIndex(
                {cell.getDimension(), cell.getIndex()},
                static_cast<Index>(local));
            x[c] += data(global);
          }
          for (size_t c = 0; c < sdim; ++c)
            pm(c, node) = x[c];
        }
        cellPointClouds[static_cast<size_t>(cellIndex)] = std::move(pm);
        hasCellPointCloud[static_cast<size_t>(cellIndex)] = 1;
      }

      for (Index vtx = 0; vtx < static_cast<Index>(mesh.getVertexCount());
           ++vtx)
      {
        auto x = mesh.getVertexCoordinates(vtx);
        for (size_t c = 0; c < sdim; ++c)
        {
          const Index global = fes.getGlobalIndex(
              {0, vtx}, static_cast<Index>(c));
          x[c] += data(global);
        }
        mesh.setVertexCoordinates(vtx, x);
      }

      // setVertexCoordinates() flushes all cached/custom transformations.
      // Reinstall the displaced P2 maps afterwards so the benchmark measures
      // the actual curved geometry, including edge-midpoint motion.
      for (Index edgeIndex = 0;
           edgeIndex < static_cast<Index>(conn.getCount(1));
           ++edgeIndex)
      {
        mesh.setPolytopeTransformation(
            {1, edgeIndex},
            new ParametricTransformation<Variational::RealH1Element<2>>(
                edgePointClouds[static_cast<size_t>(edgeIndex)],
                Variational::RealH1Element<2>(Polytope::Type::Segment)));
      }
      for (Index cellIndex = 0;
           cellIndex < static_cast<Index>(mesh.getCellCount());
           ++cellIndex)
      {
        if (!hasCellPointCloud[static_cast<size_t>(cellIndex)])
          continue;
        mesh.setPolytopeTransformation(
            {2, cellIndex},
            new ParametricTransformation<Variational::RealH1Element<2>>(
                cellPointClouds[static_cast<size_t>(cellIndex)],
                Variational::RealH1Element<2>(Polytope::Type::Triangle)));
      }
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
        Index maxIterations,
        TMOPSolveStats* tmopStats = nullptr)
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
          IdealElementTargetJacobian(mesh),
          qualityWeight);
      DeviationTerm deviation(1.0);
      AnalyticLevelSetFitTerm fit(
          phiValue, phiGradient, Optional<Attribute>(Interface), 1.0);
      AnalyticLevelSetFitTerm boundaryFit(
          boxBoundaryValue, boxBoundaryGradient, Optional<Attribute>(Boundary), 1.0);

      // Local Armijo/validity line-searched Newton. Unlike NewtonSolver
      // (constant damping only) this backtracks the step until the barrier
      // quality energy is finite and does not increase, so a barrier metric
      // (ShapeDistortion) cannot cross det<=0 and is actually usable.
      auto makeTangent = [&]()
      {
        return quality.tangent(u, du, v) + deviation.tangent(du, v)
             + fit.tangent(u, du, v) + boundaryFit.tangent(u, du, v);
      };
      auto makeResidual = [&]()
      {
        return quality.residual(u, v) + deviation.residual(u, v)
             + fit.residual(u, v) + boundaryFit.residual(u, v);
      };
      auto meritEnergy = [&]()
      {
        return quality.energy(u)
          + deviation.energy(u)
          + fit.energy(u)
          + boundaryFit.energy(u);
      };

      try
      {
        for (Index it = 0; it < maxIterations; ++it)
        {
          if (tmopStats)
            tmopStats->iterations++;
          auto start = BenchClock::now();
          LinearForm R(v);
          R = makeResidual();
          R.assemble();
          if (tmopStats)
            tmopStats->assemblySeconds += elapsedSeconds(start);
          const auto r = R.getVector();
          if (!r.allFinite())
            return false;
          if (r.norm() <= Real(1e-10))
            break;

          start = BenchClock::now();
          BilinearForm J(du, v);
          J = makeTangent();
          J.assemble();
          if (tmopStats)
            tmopStats->assemblySeconds += elapsedSeconds(start);

          start = BenchClock::now();
          Eigen::SparseLU<std::decay_t<decltype(J.getOperator())>> lu;
          lu.compute(J.getOperator());
          if (lu.info() != Eigen::Success)
            return false;
          const Math::Vector<Real> rhs = -r;          // Newton: J dx = -R
          const Math::Vector<Real> dx = lu.solve(rhs);
          if (tmopStats)
            tmopStats->solveSeconds += elapsedSeconds(start);
          if (lu.info() != Eigen::Success || !dx.allFinite())
            return false;

          start = BenchClock::now();
          const Real e0 = meritEnergy();
          if (tmopStats)
            tmopStats->meritSeconds += elapsedSeconds(start);
          const Math::Vector<Real> u0 = u.getData();
          Real alpha = Real(1);
          bool accepted = false;
          for (int ls = 0; ls < 30; ++ls)
          {
            u.getData() = u0 + alpha * dx;
            start = BenchClock::now();
            const Real e = meritEnergy();
            if (tmopStats)
              tmopStats->meritSeconds += elapsedSeconds(start);
            if (std::isfinite(e) && e <= e0 * (Real(1) + Real(1e-12)))
            {
              accepted = true;
              break;
            }
            alpha *= Real(0.5);
          }
          if (!accepted)
          {
            u.getData() = u0;                          // keep best so far
            break;
          }
          if (alpha * dx.norm() <= Real(1e-10))
            break;
        }
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
        MetricCase metricCase = MetricCase::SquaredDistance,
        TMOPSolveStats* stats = nullptr)
    {
      switch (metricCase)
      {
        case MetricCase::SquaredDistance:
          return solveTMOPWithMetric(
              mesh, shape, t, SquaredDistanceMetric{}, 0.05, 5, stats);
        case MetricCase::AreaDistortion:
          return solveTMOPWithMetric(
              mesh, shape, t, AreaDistortionMetric{}, 0.05, 5, stats);
        case MetricCase::ShapeSize:
          return solveTMOPWithMetric(
              mesh, shape, t, ShapeSizeMetric{}, 0.03, 5, stats);
        case MetricCase::ShapeDistortion:
          // The pure shape barrier has a scale nullspace and still lacks
          // line-searched step rejection in NewtonSolver; keep it in the
          // benchmark matrix, but with a conservative one-step diagnostic.
          return solveTMOPWithMetric(
              mesh, shape, t, ShapeDistortionMetric{}, 0.02, 1, stats);
      }
      return false;
    }

    template <class Metric>
    bool solveCurvedTMOPWithMetric(
        LocalMesh& mesh,
        ShapeCase shape,
        Real t,
        Metric metric,
        Real qualityWeight,
        Index maxIterations,
        TMOPSolveStats* tmopStats = nullptr)
    {
      VectorH1<2, LocalMesh> space(std::integral_constant<size_t, 2>{}, mesh, 2);
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

      CurvedQualityTargetJacobian target(mesh, 0.75);
      QualityTerm quality(
          metric,
          target,
          qualityWeight);
      quality.setQuadratureOrder(4);
      DeviationTerm deviation(0.25);
      AnalyticLevelSetFitTerm fit(
          phiValue, phiGradient, Optional<Attribute>(Interface), 2.0);
      AnalyticLevelSetFitTerm boundaryFit(
          boxBoundaryValue, boxBoundaryGradient, Optional<Attribute>(Boundary), 1.0);

      auto makeTangent = [&]()
      {
        return quality.tangent(u, du, v) + deviation.tangent(du, v)
             + fit.tangent(u, du, v) + boundaryFit.tangent(u, du, v);
      };
      auto makeResidual = [&]()
      {
        return quality.residual(u, v) + deviation.residual(u, v)
             + fit.residual(u, v) + boundaryFit.residual(u, v);
      };
      auto meritEnergy = [&]()
      {
        return quality.energy(u)
          + deviation.energy(u)
          + fit.energy(u)
          + boundaryFit.energy(u);
      };

      try
      {
        for (Index it = 0; it < maxIterations; ++it)
        {
          if (tmopStats)
            tmopStats->iterations++;
          auto start = BenchClock::now();
          LinearForm R(v);
          R = makeResidual();
          R.assemble();
          if (tmopStats)
            tmopStats->assemblySeconds += elapsedSeconds(start);
          const auto r = R.getVector();
          if (!r.allFinite())
            return false;
          if (r.norm() <= Real(1e-10))
            break;

          start = BenchClock::now();
          BilinearForm J(du, v);
          J = makeTangent();
          J.assemble();
          if (tmopStats)
            tmopStats->assemblySeconds += elapsedSeconds(start);

          start = BenchClock::now();
          Eigen::SparseLU<std::decay_t<decltype(J.getOperator())>> lu;
          lu.compute(J.getOperator());
          if (lu.info() != Eigen::Success)
            return false;
          const Math::Vector<Real> dx = lu.solve(-r);
          if (tmopStats)
            tmopStats->solveSeconds += elapsedSeconds(start);
          if (lu.info() != Eigen::Success || !dx.allFinite())
            return false;

          start = BenchClock::now();
          const Real e0 = meritEnergy();
          if (tmopStats)
            tmopStats->meritSeconds += elapsedSeconds(start);
          const Math::Vector<Real> u0 = u.getData();
          Real alpha = Real(1);
          bool accepted = false;
          for (int ls = 0; ls < 30; ++ls)
          {
            u.getData() = u0 + alpha * dx;
            start = BenchClock::now();
            const Real e = meritEnergy();
            if (tmopStats)
              tmopStats->meritSeconds += elapsedSeconds(start);
            if (std::isfinite(e) && e <= e0 * (Real(1) + Real(1e-12)))
            {
              accepted = true;
              break;
            }
            alpha *= Real(0.5);
          }
          if (!accepted)
          {
            u.getData() = u0;
            break;
          }
          if (alpha * dx.norm() <= Real(1e-10))
            break;
        }
      }
      catch (...)
      {
        return false;
      }

      if (!u.getData().allFinite())
        return false;
      LocalMesh beforeApply(mesh);
      applyParametricDisplacement(mesh, u);
      const auto curved = curvedInterfaceStats(mesh, shape, t);
      if (curved.invalidJacobianSamples != 0
          || !(curved.minJacobian > Real(0))
          || !std::isfinite(curved.minJacobian))
      {
        mesh = std::move(beforeApply);
        return false;
      }
      return true;
    }

    bool solveCurvedTMOP(
        LocalMesh& mesh,
        ShapeCase shape,
        Real t,
        MetricCase metricCase = MetricCase::ShapeSize,
        TMOPSolveStats* stats = nullptr)
    {
      switch (metricCase)
      {
        case MetricCase::SquaredDistance:
          return solveCurvedTMOPWithMetric(
              mesh, shape, t, SquaredDistanceMetric{}, 0.04, 4, stats);
        case MetricCase::AreaDistortion:
          return solveCurvedTMOPWithMetric(
              mesh, shape, t, AreaDistortionMetric{}, 0.04, 4, stats);
        case MetricCase::ShapeSize:
          return solveCurvedTMOPWithMetric(
              mesh, shape, t, ShapeSizeMetric{}, 0.03, 4, stats);
        case MetricCase::ShapeDistortion:
          return solveCurvedTMOPWithMetric(
              mesh, shape, t, ShapeDistortionMetric{}, 0.02, 1, stats);
      }
      return false;
    }

    bool solveSquaredTMOP(
        LocalMesh& mesh,
        ShapeCase shape,
        Real t,
        Real qualityWeight,
        TMOPSolveStats* stats = nullptr)
    {
      return solveTMOPWithMetric(
          mesh, shape, t, SquaredDistanceMetric{}, qualityWeight, 5, stats);
    }

    TriangleMeshOptimizerParameters defaultOptimizerParameters(
        size_t resolution)
    {
      const Real h = Real(1) / static_cast<Real>(resolution - 1);
      return TriangleMeshOptimizerParameters::levelSetCarryForward(h);
    }

    TriangleMeshOptimizerParameters optimizerParametersForStage(
        size_t resolution,
        StageCase stage)
    {
      auto parameters = defaultOptimizerParameters(resolution);
      if (stage == StageCase::OptimizeNoFeature)
        parameters.featureSmoothingPasses = 0;
      if (stage == StageCase::OptimizeNoInterior)
        parameters.smoothingPasses = 0;
      return parameters;
    }

    TriangleMeshOptimizerReport coarsen(
        LocalMesh& mesh,
        size_t resolution,
        ShapeCase shape,
        Real t,
        const TriangleMeshOptimizerParameters& parameters)
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

      // Snap a candidate point back onto phi=0 via two Newton steps on the
      // analytic level set: x <- x - phi(x) grad(x) / |grad|^2. Keeps the
      // fit exact while interface vertices are tangentially smoothed.
      auto projectToInterface =
        [shape, t](const Math::SpatialPoint& x) -> Math::SpatialPoint
      {
        return projectToLevelSet(shape, t, x, 2);
      };

      return TriangleMeshOptimizer()
        .setParameters(parameters)
        .setFeatureProjection(projectToInterface)
        .setProtectedVertices(std::move(frozen))
        .optimize(mesh);
    }

    TriangleMeshOptimizerReport coarsen(
        LocalMesh& mesh,
        size_t resolution,
        ShapeCase shape,
        Real t)
    {
      return coarsen(mesh, resolution, shape, t,
          defaultOptimizerParameters(resolution));
    }

    StageCounters collectStageCounters(
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
      const auto curved = curvedInterfaceStats(finalMesh, shape, t);
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
      counters.curvedFitRms = curved.fitRms;
      counters.curvedFitMax = curved.fitMax;
      counters.curvedMinQuality = curved.minQuality;
      counters.curvedMinJacobian = curved.minJacobian;
      counters.curvedInvalidSamples = curved.invalidJacobianSamples;
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
        counters.smooths = optReport->smooths;
        counters.featureSmooths = optReport->featureSmooths;
      }

      return counters;
    }

    void publishStageCounters(
        benchmark::State& st,
        const StageCounters& counters)
    {
      st.counters["cut_cells"] = benchmark::Counter(counters.cutCells);
      st.counters["final_cells"] = benchmark::Counter(counters.finalCells);
      st.counters["interface_edges"] = benchmark::Counter(counters.interfaceEdges);
      st.counters["uncut_cells"] = benchmark::Counter(counters.uncutCells);
      st.counters["snapped"] = benchmark::Counter(counters.snappedCrossings);
      st.counters["pathological"] = benchmark::Counter(counters.pathologicalCuts);
      st.counters["qmin_cut"] = benchmark::Counter(counters.cutMinQuality);
      st.counters["qmin_final"] = benchmark::Counter(counters.finalMinQuality);
      st.counters["avg_qmin_final"] =
        benchmark::Counter(counters.finalMinQuality);
      st.counters["final_qmin_final"] =
        benchmark::Counter(counters.finalMinQuality);
      st.counters["best_qmin_final"] =
        benchmark::Counter(counters.finalMinQuality);
      st.counters["worst_qmin_final"] =
        benchmark::Counter(counters.finalMinQuality);
      st.counters["fit_rms"] = benchmark::Counter(counters.fitRms);
      st.counters["avg_fit_rms"] = benchmark::Counter(counters.fitRms);
      st.counters["final_fit_rms"] = benchmark::Counter(counters.fitRms);
      st.counters["best_fit_rms"] = benchmark::Counter(counters.fitRms);
      st.counters["worst_fit_rms"] = benchmark::Counter(counters.fitRms);
      st.counters["fit_max"] = benchmark::Counter(counters.fitMax);
      st.counters["avg_fit_max"] = benchmark::Counter(counters.fitMax);
      st.counters["final_fit_max"] = benchmark::Counter(counters.fitMax);
      st.counters["best_fit_max"] = benchmark::Counter(counters.fitMax);
      st.counters["worst_fit_max"] = benchmark::Counter(counters.fitMax);
      st.counters["curved_fit_rms"] =
        benchmark::Counter(counters.curvedFitRms);
      st.counters["curved_qmin"] =
        benchmark::Counter(counters.curvedMinQuality);
      st.counters["avg_curved_qmin"] =
        benchmark::Counter(counters.curvedMinQuality);
      st.counters["final_curved_qmin"] =
        benchmark::Counter(counters.curvedMinQuality);
      st.counters["best_curved_qmin"] =
        benchmark::Counter(counters.curvedMinQuality);
      st.counters["worst_curved_qmin"] =
        benchmark::Counter(counters.curvedMinQuality);
      st.counters["avg_curved_fit_rms"] =
        benchmark::Counter(counters.curvedFitRms);
      st.counters["final_curved_fit_rms"] =
        benchmark::Counter(counters.curvedFitRms);
      st.counters["best_curved_fit_rms"] =
        benchmark::Counter(counters.curvedFitRms);
      st.counters["worst_curved_fit_rms"] =
        benchmark::Counter(counters.curvedFitRms);
      st.counters["curved_fit_max"] =
        benchmark::Counter(counters.curvedFitMax);
      st.counters["avg_curved_fit_max"] =
        benchmark::Counter(counters.curvedFitMax);
      st.counters["final_curved_fit_max"] =
        benchmark::Counter(counters.curvedFitMax);
      st.counters["best_curved_fit_max"] =
        benchmark::Counter(counters.curvedFitMax);
      st.counters["worst_curved_fit_max"] =
        benchmark::Counter(counters.curvedFitMax);
      st.counters["curved_min_jac"] =
        benchmark::Counter(counters.curvedMinJacobian);
      st.counters["final_curved_min_jac"] =
        benchmark::Counter(counters.curvedMinJacobian);
      st.counters["best_curved_min_jac"] =
        benchmark::Counter(counters.curvedMinJacobian);
      st.counters["worst_curved_min_jac"] =
        benchmark::Counter(counters.curvedMinJacobian);
      st.counters["curved_invalid"] =
        benchmark::Counter(counters.curvedInvalidSamples);
      st.counters["avg_curved_invalid"] =
        benchmark::Counter(counters.curvedInvalidSamples);
      st.counters["final_curved_invalid"] =
        benchmark::Counter(counters.curvedInvalidSamples);
      st.counters["best_curved_invalid"] =
        benchmark::Counter(counters.curvedInvalidSamples);
      st.counters["worst_curved_invalid"] =
        benchmark::Counter(counters.curvedInvalidSamples);
      if (counters.hasP2Diagnostics)
      {
        st.counters["p2_fit_rms"] =
          benchmark::Counter(counters.curvedFitRms);
        st.counters["p2_qmin"] =
          benchmark::Counter(counters.curvedMinQuality);
        st.counters["avg_p2_qmin"] =
          benchmark::Counter(counters.curvedMinQuality);
        st.counters["final_p2_qmin"] =
          benchmark::Counter(counters.curvedMinQuality);
        st.counters["best_p2_qmin"] =
          benchmark::Counter(counters.curvedMinQuality);
        st.counters["worst_p2_qmin"] =
          benchmark::Counter(counters.curvedMinQuality);
        st.counters["avg_p2_fit_rms"] =
          benchmark::Counter(counters.curvedFitRms);
        st.counters["final_p2_fit_rms"] =
          benchmark::Counter(counters.curvedFitRms);
        st.counters["best_p2_fit_rms"] =
          benchmark::Counter(counters.curvedFitRms);
        st.counters["worst_p2_fit_rms"] =
          benchmark::Counter(counters.curvedFitRms);
        st.counters["p2_fit_max"] =
          benchmark::Counter(counters.curvedFitMax);
        st.counters["avg_p2_fit_max"] =
          benchmark::Counter(counters.curvedFitMax);
        st.counters["final_p2_fit_max"] =
          benchmark::Counter(counters.curvedFitMax);
        st.counters["best_p2_fit_max"] =
          benchmark::Counter(counters.curvedFitMax);
        st.counters["worst_p2_fit_max"] =
          benchmark::Counter(counters.curvedFitMax);
        st.counters["p2_min_jac"] =
          benchmark::Counter(counters.curvedMinJacobian);
        st.counters["final_p2_min_jac"] =
          benchmark::Counter(counters.curvedMinJacobian);
        st.counters["best_p2_min_jac"] =
          benchmark::Counter(counters.curvedMinJacobian);
        st.counters["worst_p2_min_jac"] =
          benchmark::Counter(counters.curvedMinJacobian);
        st.counters["p2_invalid"] =
          benchmark::Counter(counters.curvedInvalidSamples);
        st.counters["avg_p2_invalid"] =
          benchmark::Counter(counters.curvedInvalidSamples);
        st.counters["final_p2_invalid"] =
          benchmark::Counter(counters.curvedInvalidSamples);
        st.counters["best_p2_invalid"] =
          benchmark::Counter(counters.curvedInvalidSamples);
        st.counters["worst_p2_invalid"] =
          benchmark::Counter(counters.curvedInvalidSamples);
      }
      st.counters["coverage"] = benchmark::Counter(counters.coverage);
      st.counters["avg_coverage"] = benchmark::Counter(counters.coverage);
      st.counters["final_coverage"] = benchmark::Counter(counters.coverage);
      st.counters["best_coverage"] = benchmark::Counter(counters.coverage);
      st.counters["worst_coverage"] = benchmark::Counter(counters.coverage);
      st.counters["signed_area"] = benchmark::Counter(counters.signedArea);
      st.counters["max_iface_dev"] =
        benchmark::Counter(counters.maxInterfaceDeviation);
      st.counters["cut_ms"] =
        benchmark::Counter(Real(1000) * counters.cutSeconds);
      st.counters["optimize_ms"] =
        benchmark::Counter(Real(1000) * counters.optimizeSeconds);
      st.counters["tmop_ms"] =
        benchmark::Counter(Real(1000) * counters.tmopSeconds);
      st.counters["tmop_assembly_ms"] =
        benchmark::Counter(Real(1000) * counters.tmopAssemblySeconds);
      st.counters["tmop_solve_ms"] =
        benchmark::Counter(Real(1000) * counters.tmopSolveSeconds);
      st.counters["tmop_merit_ms"] =
        benchmark::Counter(Real(1000) * counters.tmopMeritSeconds);
      st.counters["tmop_iterations"] =
        benchmark::Counter(counters.tmopIterations);
      st.counters["inverted"] = benchmark::Counter(counters.inverted);
      st.counters["tmop_failed"] = benchmark::Counter(counters.tmopFailures);
      st.counters["splits"] = benchmark::Counter(counters.splits);
      st.counters["collapses"] = benchmark::Counter(counters.collapses);
      st.counters["swaps"] = benchmark::Counter(counters.swaps);
      st.counters["smooths"] = benchmark::Counter(counters.smooths);
      st.counters["feature_smooths"] =
        benchmark::Counter(counters.featureSmooths);
    }

    StageCounters setCounters(
        benchmark::State& st,
        const LevelSetDiscretizerTrianglesResult& cut,
        const LocalMesh& finalMesh,
        const TriangleMeshOptimizerReport* optReport,
        ShapeCase shape,
        bool tmopSucceeded,
        Real t,
        Real cutSeconds = 0,
        Real optimizeSeconds = 0,
        Real tmopSeconds = 0,
        const TMOPSolveStats* tmopStats = nullptr,
        bool hasP2Diagnostics = false)
    {
      auto counters = collectStageCounters(
          cut, finalMesh, optReport, shape, tmopSucceeded, t);
      counters.cutSeconds = cutSeconds;
      counters.optimizeSeconds = optimizeSeconds;
      counters.tmopSeconds = tmopSeconds;
      counters.hasP2Diagnostics = hasP2Diagnostics;
      if (tmopStats)
      {
        counters.tmopAssemblySeconds = tmopStats->assemblySeconds;
        counters.tmopSolveSeconds = tmopStats->solveSeconds;
        counters.tmopMeritSeconds = tmopStats->meritSeconds;
        counters.tmopIterations = tmopStats->iterations;
      }
      publishStageCounters(st, counters);
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
      const auto cutStart = BenchClock::now();
      auto cut = cutShape(background, shape, Real(0.25), Real(0), Real(0));
      const Real cutSeconds = elapsedSeconds(cutStart);
      benchmark::DoNotOptimize(cut.mesh.getCellCount());
      counters =
        setCounters(st, cut, cut.mesh, nullptr, shape, true, Real(0.25),
            cutSeconds);
    }
    st.counters["resolution"] = benchmark::Counter(resolution);
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
      const auto cutStart = BenchClock::now();
      auto cut = cutShape(background, shape, Real(0.25), Real(0.05), Real(0.2));
      const Real cutSeconds = elapsedSeconds(cutStart);
      benchmark::DoNotOptimize(cut.mesh.getCellCount());
      counters =
        setCounters(st, cut, cut.mesh, nullptr, shape, true, Real(0.25),
            cutSeconds);
    }
    st.counters["resolution"] = benchmark::Counter(resolution);
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
      const auto cutStart = BenchClock::now();
      auto cut = cutShape(background, shape, Real(0.25), Real(0.05), Real(0.2));
      const Real cutSeconds = elapsedSeconds(cutStart);
      auto optimized = cut.mesh;
      TMOPSolveStats tmopStats;
      const auto tmopStart = BenchClock::now();
      const bool tmopOk =
        solveTMOP(optimized, shape, Real(0.25), metric, &tmopStats);
      const Real tmopSeconds = elapsedSeconds(tmopStart);
      benchmark::DoNotOptimize(optimized.getCellCount());
      counters =
        setCounters(st, cut, optimized, nullptr, shape, tmopOk, Real(0.25),
            cutSeconds, 0, tmopSeconds, &tmopStats);
    }
    st.counters["resolution"] = benchmark::Counter(resolution);
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
      const auto cutStart = BenchClock::now();
      auto cut = cutShape(background, shape, Real(0.25), Real(0.05), Real(0.2));
      const Real cutSeconds = elapsedSeconds(cutStart);
      auto optimized = cut.mesh;
      TMOPSolveStats tmopStats;
      const auto tmopStart = BenchClock::now();
      const bool tmopOk =
        solveTMOP(optimized, shape, Real(0.25),
            MetricCase::SquaredDistance, &tmopStats);
      const Real tmopSeconds = elapsedSeconds(tmopStart);
      const auto optStart = BenchClock::now();
      const auto report = coarsen(optimized, resolution, shape, Real(0.25));
      const Real optimizeSeconds = elapsedSeconds(optStart);
      benchmark::DoNotOptimize(optimized.getCellCount());
      counters =
        setCounters(st, cut, optimized, &report, shape, tmopOk, Real(0.25),
            cutSeconds, optimizeSeconds, tmopSeconds, &tmopStats);
    }
    st.counters["resolution"] = benchmark::Counter(resolution);
    st.SetItemsProcessed(st.iterations() * counters.finalCells);
  }

  static void BM_LevelSetPipeline_StageFreshCut(benchmark::State& st)
  {
    const auto resolution = static_cast<size_t>(st.range(0));
    const auto shape = static_cast<ShapeCase>(st.range(1));
    const auto metric = static_cast<MetricCase>(st.range(2));
    const auto stage = static_cast<StageCase>(st.range(3));
    StageCounters counters;
    StageAverages averages;
    for (auto _ : st)
    {
      averages = StageAverages{};
      auto background = makeBackground(resolution);
      const auto cutStart = BenchClock::now();
      auto cut = cutShape(background, shape, Real(0.25), Real(0.05), Real(0.2));
      const Real cutSeconds = elapsedSeconds(cutStart);
      auto optimized = cut.mesh;
      bool tmopOk = true;
      TriangleMeshOptimizerReport report;
      TMOPSolveStats tmopStats;
      Real optimizeSeconds = 0;
      Real tmopSeconds = 0;
      const auto config = optimizerParametersForStage(resolution, stage);
      switch (stage)
      {
        case StageCase::TMOPOnly:
        {
          const auto start = BenchClock::now();
          tmopOk = solveTMOP(
              optimized, shape, Real(0.25), metric, &tmopStats);
          tmopSeconds += elapsedSeconds(start);
          break;
        }
        case StageCase::OptimizeOnly:
        {
          const auto start = BenchClock::now();
          report = coarsen(optimized, resolution, shape, Real(0.25), config);
          optimizeSeconds += elapsedSeconds(start);
          break;
        }
        case StageCase::TMOPThenOptimize:
        case StageCase::OptimizeNoFeature:
        case StageCase::OptimizeNoInterior:
        {
          auto start = BenchClock::now();
          tmopOk = solveTMOP(
              optimized, shape, Real(0.25), metric, &tmopStats);
          tmopSeconds += elapsedSeconds(start);
          start = BenchClock::now();
          report = coarsen(optimized, resolution, shape, Real(0.25), config);
          optimizeSeconds += elapsedSeconds(start);
          break;
        }
        case StageCase::OptimizeThenTMOP:
        {
          auto start = BenchClock::now();
          report = coarsen(optimized, resolution, shape, Real(0.25), config);
          optimizeSeconds += elapsedSeconds(start);
          start = BenchClock::now();
          tmopOk = solveTMOP(
              optimized, shape, Real(0.25), metric, &tmopStats);
          tmopSeconds += elapsedSeconds(start);
          break;
        }
      }
      benchmark::DoNotOptimize(optimized.getCellCount());
      counters =
        setCounters(st, cut, optimized, &report, shape, tmopOk, Real(0.25),
            cutSeconds, optimizeSeconds, tmopSeconds, &tmopStats);
      averages.add(counters);
    }
    averages.publish(st);
    st.counters["resolution"] = benchmark::Counter(resolution);
    st.counters["stage_case"] = benchmark::Counter(static_cast<int>(stage));
    st.counters["production_order"] =
      benchmark::Counter(stage == ProductionStageCase ? 1 : 0);
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
      const auto cutStart = BenchClock::now();
      auto cut = cutShape(background, shape, Real(0.25), Real(0.05), Real(0.2));
      const Real cutSeconds = elapsedSeconds(cutStart);
      auto optimized = cut.mesh;
      TMOPSolveStats tmopStats;
      const auto tmopStart = BenchClock::now();
      const bool tmopOk =
        solveSquaredTMOP(
            optimized, shape, Real(0.25), qualityWeight, &tmopStats);
      const Real tmopSeconds = elapsedSeconds(tmopStart);
      benchmark::DoNotOptimize(optimized.getCellCount());
      counters =
        setCounters(st, cut, optimized, nullptr, shape, tmopOk, Real(0.25),
            cutSeconds, 0, tmopSeconds, &tmopStats);
    }
    st.counters["quality_weight"] = benchmark::Counter(qualityWeight);
    st.counters["resolution"] = benchmark::Counter(resolution);
    st.SetItemsProcessed(st.iterations() * counters.cutCells);
  }

  static void BM_LevelSetPipeline_CarryForwardOrbit(benchmark::State& st)
  {
    const auto resolution = static_cast<size_t>(st.range(0));
    const auto steps = static_cast<Index>(st.range(1));
    const auto shape = static_cast<ShapeCase>(st.range(2));
    const auto metric = static_cast<MetricCase>(st.range(3));
    StageCounters counters;
    StageAverages averages;
    for (auto _ : st)
    {
      averages = StageAverages{};
      auto background = makeBackground(resolution);
      Index tmopFailures = 0;
      for (Index step = 0; step < steps; ++step)
      {
        const Real t = static_cast<Real>(step)
          / static_cast<Real>(std::max<Index>(1, steps - 1));
        const auto cutStart = BenchClock::now();
        auto cut = cutShape(background, shape, t, Real(0.05), Real(0.2));
        const Real cutSeconds = elapsedSeconds(cutStart);
        auto optimized = cut.mesh;
        const auto optStart = BenchClock::now();
        const auto report = coarsen(optimized, resolution, shape, t);
        const Real optimizeSeconds = elapsedSeconds(optStart);
        TMOPSolveStats tmopStats;
        const auto tmopStart = BenchClock::now();
        const bool tmopOk =
          solveTMOP(optimized, shape, t, metric, &tmopStats);
        const Real tmopSeconds = elapsedSeconds(tmopStart);
        if (!tmopOk)
          tmopFailures++;
        optimized.getConnectivity().compute(2, 1);
        optimized.getConnectivity().compute(1, 0);
        counters = setCounters(st, cut, optimized, &report, shape, tmopOk, t,
            cutSeconds, optimizeSeconds, tmopSeconds, &tmopStats);
        averages.add(counters);
        demoteInterface(optimized);
        background = std::move(optimized);
      }
      counters.tmopFailures = tmopFailures;
      benchmark::DoNotOptimize(background.getCellCount());
    }
    averages.publish(st);
    st.counters["resolution"] = benchmark::Counter(resolution);
    st.counters["orbit_steps"] = benchmark::Counter(steps);
    st.counters["production_order"] = benchmark::Counter(1);
    st.counters["tmop_failed"] = benchmark::Counter(counters.tmopFailures);
    st.SetItemsProcessed(st.iterations() * counters.finalCells);
  }

  static void BM_LevelSetPipeline_CurvedFreshCut(benchmark::State& st)
  {
    const auto resolution = static_cast<size_t>(st.range(0));
    const auto shape = static_cast<ShapeCase>(st.range(1));
    const auto metric = static_cast<MetricCase>(st.range(2));
    StageCounters counters;
    StageAverages averages;
    for (auto _ : st)
    {
      averages = StageAverages{};
      auto background = makeBackground(resolution);
      const auto cutStart = BenchClock::now();
      auto cut = cutShape(background, shape, Real(0.25), Real(0.05), Real(0.2));
      const Real cutSeconds = elapsedSeconds(cutStart);
      auto optimized = cut.mesh;

      const auto optStart = BenchClock::now();
      const auto report = coarsen(optimized, resolution, shape, Real(0.25));
      const Real optimizeSeconds = elapsedSeconds(optStart);

      optimized.getConnectivity().compute(2, 1);
      optimized.getConnectivity().compute(1, 0);
      auto curved = optimized;
      installQuadraticInterfaceGeometry(curved, shape, Real(0.25));

      TMOPSolveStats tmopStats;
      const auto tmopStart = BenchClock::now();
      const bool tmopOk =
        solveCurvedTMOP(curved, shape, Real(0.25), metric, &tmopStats);
      const Real tmopSeconds = elapsedSeconds(tmopStart);

      counters = setCounters(st, cut, curved, &report, shape, tmopOk,
          Real(0.25), cutSeconds, optimizeSeconds, tmopSeconds, &tmopStats,
          true);
      averages.add(counters);
      benchmark::DoNotOptimize(curved.getCellCount());
    }
    averages.publish(st);
    st.counters["resolution"] = benchmark::Counter(resolution);
    st.counters["curved_geometry"] = benchmark::Counter(1);
    st.counters["production_order"] = benchmark::Counter(1);
    st.SetItemsProcessed(st.iterations() * counters.finalCells);
  }

  static void BM_LevelSetPipeline_CurvedCarryForwardOrbit(
      benchmark::State& st)
  {
    const auto resolution = static_cast<size_t>(st.range(0));
    const auto steps = static_cast<Index>(st.range(1));
    const auto shape = static_cast<ShapeCase>(st.range(2));
    const auto metric = static_cast<MetricCase>(st.range(3));
    StageCounters counters;
    StageAverages averages;
    for (auto _ : st)
    {
      averages = StageAverages{};
      auto background = makeBackground(resolution);
      Index tmopFailures = 0;
      for (Index step = 0; step < steps; ++step)
      {
        const Real t = static_cast<Real>(step)
          / static_cast<Real>(std::max<Index>(1, steps - 1));
        const auto cutStart = BenchClock::now();
        auto cut = cutShape(background, shape, t, Real(0.05), Real(0.2));
        const Real cutSeconds = elapsedSeconds(cutStart);
        auto optimized = cut.mesh;

        const auto optStart = BenchClock::now();
        const auto report = coarsen(optimized, resolution, shape, t);
        const Real optimizeSeconds = elapsedSeconds(optStart);

        optimized.getConnectivity().compute(2, 1);
        optimized.getConnectivity().compute(1, 0);
        auto curved = optimized;
        installQuadraticInterfaceGeometry(curved, shape, t);

        TMOPSolveStats tmopStats;
        const auto tmopStart = BenchClock::now();
        const bool tmopOk =
          solveCurvedTMOP(curved, shape, t, metric, &tmopStats);
        const Real tmopSeconds = elapsedSeconds(tmopStart);
        if (!tmopOk)
          tmopFailures++;

        counters = setCounters(st, cut, curved, &report, shape, tmopOk, t,
            cutSeconds, optimizeSeconds, tmopSeconds, &tmopStats, true);
        averages.add(counters);

        // Carry the endpoint motion produced by the curved TMOP solve into
        // the next linear cut. The cutter is still P1/topological, so clear
        // the P2 transformations after preserving the updated vertices.
        demoteInterface(curved);
        curved.flush();
        background = std::move(curved);
      }
      counters.tmopFailures = tmopFailures;
      benchmark::DoNotOptimize(background.getCellCount());
    }
    averages.publish(st);
    st.counters["resolution"] = benchmark::Counter(resolution);
    st.counters["orbit_steps"] = benchmark::Counter(steps);
    st.counters["curved_geometry"] = benchmark::Counter(1);
    st.counters["production_order"] = benchmark::Counter(1);
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
    StageAverages averages;
    for (auto _ : st)
    {
      averages = StageAverages{};
      auto background = makeBackground(resolution);
      Index tmopFailures = 0;
      for (Index step = 0; step < steps; ++step)
      {
        const Real t = static_cast<Real>(step)
          / static_cast<Real>(std::max<Index>(1, steps - 1));
        const auto cutStart = BenchClock::now();
        auto cut = cutShape(background, shape, t, Real(0.05), Real(0.2));
        const Real cutSeconds = elapsedSeconds(cutStart);
        auto optimized = cut.mesh;
        const auto optStart = BenchClock::now();
        const auto report = coarsen(optimized, resolution, shape, t);
        const Real optimizeSeconds = elapsedSeconds(optStart);
        TMOPSolveStats tmopStats;
        const auto tmopStart = BenchClock::now();
        const bool tmopOk =
          solveSquaredTMOP(
              optimized, shape, t, qualityWeight, &tmopStats);
        const Real tmopSeconds = elapsedSeconds(tmopStart);
        if (!tmopOk)
          tmopFailures++;
        optimized.getConnectivity().compute(2, 1);
        optimized.getConnectivity().compute(1, 0);
        counters = setCounters(st, cut, optimized, &report, shape, tmopOk, t,
            cutSeconds, optimizeSeconds, tmopSeconds, &tmopStats);
        averages.add(counters);
        demoteInterface(optimized);
        background = std::move(optimized);
      }
      counters.tmopFailures = tmopFailures;
      benchmark::DoNotOptimize(background.getCellCount());
    }
    averages.publish(st);
    st.counters["resolution"] = benchmark::Counter(resolution);
    st.counters["orbit_steps"] = benchmark::Counter(steps);
    st.counters["quality_weight"] = benchmark::Counter(qualityWeight);
    st.counters["production_order"] = benchmark::Counter(1);
    st.counters["tmop_failed"] = benchmark::Counter(counters.tmopFailures);
    st.SetItemsProcessed(st.iterations() * counters.finalCells);
  }

  static void BM_LevelSetPipeline_CarryForwardStage(benchmark::State& st)
  {
    const auto resolution = static_cast<size_t>(st.range(0));
    const auto steps = static_cast<Index>(st.range(1));
    const auto shape = static_cast<ShapeCase>(st.range(2));
    const auto metric = static_cast<MetricCase>(st.range(3));
    const auto stage = static_cast<StageCase>(st.range(4));
    StageCounters counters;
    StageAverages averages;
    for (auto _ : st)
    {
      averages = StageAverages{};
      auto background = makeBackground(resolution);
      Index tmopFailures = 0;
      TriangleMeshOptimizerReport accumulated;
      for (Index step = 0; step < steps; ++step)
      {
        const Real t = static_cast<Real>(step)
          / static_cast<Real>(std::max<Index>(1, steps - 1));
        const auto cutStart = BenchClock::now();
        auto cut = cutShape(background, shape, t, Real(0.05), Real(0.2));
        const Real cutSeconds = elapsedSeconds(cutStart);
        auto optimized = cut.mesh;
        bool tmopOk = true;
        TriangleMeshOptimizerReport report;
        TMOPSolveStats tmopStats;
        Real optimizeSeconds = 0;
        Real tmopSeconds = 0;
        const auto config = optimizerParametersForStage(resolution, stage);
        switch (stage)
        {
          case StageCase::TMOPOnly:
          {
            const auto start = BenchClock::now();
            tmopOk = solveTMOP(optimized, shape, t, metric, &tmopStats);
            tmopSeconds += elapsedSeconds(start);
            break;
          }
          case StageCase::OptimizeOnly:
          {
            const auto start = BenchClock::now();
            report = coarsen(optimized, resolution, shape, t, config);
            optimizeSeconds += elapsedSeconds(start);
            break;
          }
          case StageCase::TMOPThenOptimize:
          case StageCase::OptimizeNoFeature:
          case StageCase::OptimizeNoInterior:
          {
            auto start = BenchClock::now();
            tmopOk = solveTMOP(optimized, shape, t, metric, &tmopStats);
            tmopSeconds += elapsedSeconds(start);
            start = BenchClock::now();
            report = coarsen(optimized, resolution, shape, t, config);
            optimizeSeconds += elapsedSeconds(start);
            break;
          }
          case StageCase::OptimizeThenTMOP:
          {
            auto start = BenchClock::now();
            report = coarsen(optimized, resolution, shape, t, config);
            optimizeSeconds += elapsedSeconds(start);
            start = BenchClock::now();
            tmopOk = solveTMOP(optimized, shape, t, metric, &tmopStats);
            tmopSeconds += elapsedSeconds(start);
            break;
          }
        }
        if (!tmopOk)
          tmopFailures++;
        accumulated.splits += report.splits;
        accumulated.collapses += report.collapses;
        accumulated.swaps += report.swaps;
        accumulated.smooths += report.smooths;
        accumulated.featureSmooths += report.featureSmooths;
        optimized.getConnectivity().compute(2, 1);
        optimized.getConnectivity().compute(1, 0);
        counters = setCounters(st, cut, optimized, &accumulated, shape, tmopOk, t,
            cutSeconds, optimizeSeconds, tmopSeconds, &tmopStats);
        auto stepCounters =
          collectStageCounters(cut, optimized, &report, shape, tmopOk, t);
        stepCounters.cutSeconds = cutSeconds;
        stepCounters.optimizeSeconds = optimizeSeconds;
        stepCounters.tmopSeconds = tmopSeconds;
        stepCounters.tmopAssemblySeconds = tmopStats.assemblySeconds;
        stepCounters.tmopSolveSeconds = tmopStats.solveSeconds;
        stepCounters.tmopMeritSeconds = tmopStats.meritSeconds;
        stepCounters.tmopIterations = tmopStats.iterations;
        averages.add(stepCounters);
        demoteInterface(optimized);
        background = std::move(optimized);
      }
      counters.tmopFailures = tmopFailures;
      benchmark::DoNotOptimize(background.getCellCount());
    }
    averages.publish(st);
    st.counters["resolution"] = benchmark::Counter(resolution);
    st.counters["orbit_steps"] = benchmark::Counter(steps);
    st.counters["stage_case"] = benchmark::Counter(static_cast<int>(stage));
    st.counters["production_order"] =
      benchmark::Counter(stage == ProductionStageCase ? 1 : 0);
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

#define RODIN_RESOLUTION_ORBIT_BENCH(SHAPE_LABEL, SHAPE_ENUM, METRIC_LABEL, METRIC_ENUM, RESOLUTION, STEPS) \
  BENCHMARK(BM_LevelSetPipeline_CarryForwardOrbit) \
    ->Name("LevelSetPipeline/Resolution/" SHAPE_LABEL "/" METRIC_LABEL "/" #RESOLUTION "x" #STEPS) \
    ->Args({RESOLUTION, STEPS, static_cast<int>(ShapeCase::SHAPE_ENUM), \
      static_cast<int>(MetricCase::METRIC_ENUM)}) \
    ->Unit(benchmark::kMillisecond);

#define RODIN_CURVED_FRESH_BENCH(SHAPE_LABEL, SHAPE_ENUM, METRIC_LABEL, METRIC_ENUM, RESOLUTION) \
  BENCHMARK(BM_LevelSetPipeline_CurvedFreshCut) \
    ->Name("LevelSetPipeline/CurvedFresh/" SHAPE_LABEL "/" METRIC_LABEL "/" #RESOLUTION) \
    ->Args({RESOLUTION, static_cast<int>(ShapeCase::SHAPE_ENUM), \
      static_cast<int>(MetricCase::METRIC_ENUM)}) \
    ->Unit(benchmark::kMillisecond);

#define RODIN_CURVED_ORBIT_BENCH(SHAPE_LABEL, SHAPE_ENUM, METRIC_LABEL, METRIC_ENUM, RESOLUTION, STEPS) \
  BENCHMARK(BM_LevelSetPipeline_CurvedCarryForwardOrbit) \
    ->Name("LevelSetPipeline/CurvedCarryForward/" SHAPE_LABEL "/" METRIC_LABEL "/" #RESOLUTION "x" #STEPS) \
    ->Args({RESOLUTION, STEPS, static_cast<int>(ShapeCase::SHAPE_ENUM), \
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

#define RODIN_STAGE_FRESH_BENCH(SHAPE_LABEL, SHAPE_ENUM, METRIC_LABEL, METRIC_ENUM, STAGE_LABEL, STAGE_ENUM) \
  BENCHMARK(BM_LevelSetPipeline_StageFreshCut) \
    ->Name("LevelSetPipeline/StageFresh/" SHAPE_LABEL "/" METRIC_LABEL "/" STAGE_LABEL "/24") \
    ->Args({24, static_cast<int>(ShapeCase::SHAPE_ENUM), \
      static_cast<int>(MetricCase::METRIC_ENUM), \
      static_cast<int>(StageCase::STAGE_ENUM)}) \
    ->Unit(benchmark::kMillisecond);

#define RODIN_STAGE_ORBIT_BENCH(SHAPE_LABEL, SHAPE_ENUM, METRIC_LABEL, METRIC_ENUM, STAGE_LABEL, STAGE_ENUM, STEPS) \
  BENCHMARK(BM_LevelSetPipeline_CarryForwardStage) \
    ->Name("LevelSetPipeline/CarryForwardStage/" SHAPE_LABEL "/" METRIC_LABEL "/" STAGE_LABEL "/24x" #STEPS) \
    ->Args({24, STEPS, static_cast<int>(ShapeCase::SHAPE_ENUM), \
      static_cast<int>(MetricCase::METRIC_ENUM), \
      static_cast<int>(StageCase::STAGE_ENUM)}) \
    ->Unit(benchmark::kMillisecond);

  RODIN_CUT_EXACT_BENCH("Circle", CircleOrbit)
  RODIN_CUT_EXACT_BENCH("Square", SquareOrbit)
  RODIN_CUT_EXACT_BENCH("Triangle", TriangleOrbit)
  RODIN_CUT_EXACT_BENCH("CircleEnterLeave", CircleEnterLeave)
  RODIN_CUT_EXACT_BENCH("CircleFigureEight", CircleFigureEight)
  RODIN_CUT_EXACT_BENCH("SquareFigureEight", SquareFigureEight)
  RODIN_CUT_EXACT_BENCH("TriangleFigureEight", TriangleFigureEight)
  RODIN_CUT_EXACT_BENCH("WavyCircle", WavyCircleOrbit)
  RODIN_CUT_EXACT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight)
  RODIN_CUT_EXACT_BENCH("WavyCircleEnterLeave", WavyCircleEnterLeave)

  RODIN_CUT_QUALITY_BENCH("Circle", CircleOrbit, 24)
  RODIN_CUT_QUALITY_BENCH("Circle", CircleOrbit, 48)
  RODIN_CUT_QUALITY_BENCH("Circle", CircleOrbit, 10)
  RODIN_CUT_QUALITY_BENCH("Circle", CircleOrbit, 12)
  RODIN_CUT_QUALITY_BENCH("Circle", CircleOrbit, 16)
  RODIN_CUT_QUALITY_BENCH("Circle", CircleOrbit, 32)
  RODIN_CUT_QUALITY_BENCH("Square", SquareOrbit, 24)
  RODIN_CUT_QUALITY_BENCH("Triangle", TriangleOrbit, 24)
  RODIN_CUT_QUALITY_BENCH("CircleEnterLeave", CircleEnterLeave, 24)
  RODIN_CUT_QUALITY_BENCH("CircleFigureEight", CircleFigureEight, 24)
  RODIN_CUT_QUALITY_BENCH("SquareFigureEight", SquareFigureEight, 24)
  RODIN_CUT_QUALITY_BENCH("TriangleFigureEight", TriangleFigureEight, 24)
  RODIN_CUT_QUALITY_BENCH("WavyCircle", WavyCircleOrbit, 24)
  RODIN_CUT_QUALITY_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, 24)
  RODIN_CUT_QUALITY_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, 10)
  RODIN_CUT_QUALITY_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, 12)
  RODIN_CUT_QUALITY_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, 16)
  RODIN_CUT_QUALITY_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, 32)
  RODIN_CUT_QUALITY_BENCH("WavyCircleEnterLeave", WavyCircleEnterLeave, 24)

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
  RODIN_TMOP_METRIC_BENCH("CircleFigureEight", CircleFigureEight, "SquaredDistance", SquaredDistance)
  RODIN_TMOP_METRIC_BENCH("SquareFigureEight", SquareFigureEight, "SquaredDistance", SquaredDistance)
  RODIN_TMOP_METRIC_BENCH("TriangleFigureEight", TriangleFigureEight, "SquaredDistance", SquaredDistance)
  RODIN_TMOP_METRIC_BENCH("WavyCircle", WavyCircleOrbit, "SquaredDistance", SquaredDistance)
  RODIN_TMOP_METRIC_BENCH("WavyCircle", WavyCircleOrbit, "AreaDistortion", AreaDistortion)
  RODIN_TMOP_METRIC_BENCH("WavyCircle", WavyCircleOrbit, "ShapeSize", ShapeSize)
  RODIN_TMOP_METRIC_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "SquaredDistance", SquaredDistance)
  RODIN_TMOP_METRIC_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "AreaDistortion", AreaDistortion)
  RODIN_TMOP_METRIC_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "ShapeSize", ShapeSize)

  RODIN_WEIGHT_FRESH_BENCH("Circle", CircleOrbit, 50)
  RODIN_WEIGHT_FRESH_BENCH("Circle", CircleOrbit, 100)
  RODIN_WEIGHT_FRESH_BENCH("Circle", CircleOrbit, 200)
  RODIN_WEIGHT_FRESH_BENCH("Circle", CircleOrbit, 500)

  RODIN_COARSEN_BENCH("Circle", CircleOrbit)
  RODIN_COARSEN_BENCH("Square", SquareOrbit)
  RODIN_COARSEN_BENCH("Triangle", TriangleOrbit)
  RODIN_COARSEN_BENCH("CircleEnterLeave", CircleEnterLeave)
  RODIN_COARSEN_BENCH("CircleFigureEight", CircleFigureEight)
  RODIN_COARSEN_BENCH("SquareFigureEight", SquareFigureEight)
  RODIN_COARSEN_BENCH("TriangleFigureEight", TriangleFigureEight)
  RODIN_COARSEN_BENCH("WavyCircle", WavyCircleOrbit)
  RODIN_COARSEN_BENCH("WavyCircleFigureEight", WavyCircleFigureEight)

  RODIN_STAGE_FRESH_BENCH("Circle", CircleOrbit, "SquaredDistance", SquaredDistance, "TMOPOnly", TMOPOnly)
  RODIN_STAGE_FRESH_BENCH("Circle", CircleOrbit, "SquaredDistance", SquaredDistance, "OptimizeOnly", OptimizeOnly)
  RODIN_STAGE_FRESH_BENCH("Circle", CircleOrbit, "SquaredDistance", SquaredDistance, "TMOPThenOptimize", TMOPThenOptimize)
  RODIN_STAGE_FRESH_BENCH("Circle", CircleOrbit, "SquaredDistance", SquaredDistance, "OptimizeThenTMOP", OptimizeThenTMOP)
  RODIN_STAGE_FRESH_BENCH("Circle", CircleOrbit, "SquaredDistance", SquaredDistance, "OptimizeNoFeature", OptimizeNoFeature)
  RODIN_STAGE_FRESH_BENCH("Circle", CircleOrbit, "SquaredDistance", SquaredDistance, "OptimizeNoInterior", OptimizeNoInterior)
  RODIN_STAGE_FRESH_BENCH("Circle", CircleOrbit, "AreaDistortion", AreaDistortion, "TMOPThenOptimize", TMOPThenOptimize)
  RODIN_STAGE_FRESH_BENCH("Circle", CircleOrbit, "ShapeSize", ShapeSize, "TMOPThenOptimize", TMOPThenOptimize)
  RODIN_STAGE_FRESH_BENCH("Square", SquareOrbit, "SquaredDistance", SquaredDistance, "TMOPThenOptimize", TMOPThenOptimize)
  RODIN_STAGE_FRESH_BENCH("Square", SquareOrbit, "AreaDistortion", AreaDistortion, "TMOPThenOptimize", TMOPThenOptimize)
  RODIN_STAGE_FRESH_BENCH("Square", SquareOrbit, "ShapeSize", ShapeSize, "TMOPThenOptimize", TMOPThenOptimize)
  RODIN_STAGE_FRESH_BENCH("Square", SquareOrbit, "SquaredDistance", SquaredDistance, "OptimizeOnly", OptimizeOnly)
  RODIN_STAGE_FRESH_BENCH("Triangle", TriangleOrbit, "SquaredDistance", SquaredDistance, "TMOPThenOptimize", TMOPThenOptimize)
  RODIN_STAGE_FRESH_BENCH("Triangle", TriangleOrbit, "AreaDistortion", AreaDistortion, "TMOPThenOptimize", TMOPThenOptimize)
  RODIN_STAGE_FRESH_BENCH("Triangle", TriangleOrbit, "ShapeSize", ShapeSize, "TMOPThenOptimize", TMOPThenOptimize)
  RODIN_STAGE_FRESH_BENCH("Triangle", TriangleOrbit, "SquaredDistance", SquaredDistance, "OptimizeOnly", OptimizeOnly)
  RODIN_STAGE_FRESH_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "SquaredDistance", SquaredDistance, "TMOPOnly", TMOPOnly)
  RODIN_STAGE_FRESH_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "SquaredDistance", SquaredDistance, "OptimizeOnly", OptimizeOnly)
  RODIN_STAGE_FRESH_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "SquaredDistance", SquaredDistance, "TMOPThenOptimize", TMOPThenOptimize)
  RODIN_STAGE_FRESH_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "SquaredDistance", SquaredDistance, "OptimizeThenTMOP", OptimizeThenTMOP)
  RODIN_STAGE_FRESH_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "SquaredDistance", SquaredDistance, "OptimizeNoFeature", OptimizeNoFeature)
  RODIN_STAGE_FRESH_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "SquaredDistance", SquaredDistance, "OptimizeNoInterior", OptimizeNoInterior)
  RODIN_STAGE_FRESH_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "AreaDistortion", AreaDistortion, "TMOPThenOptimize", TMOPThenOptimize)
  RODIN_STAGE_FRESH_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "ShapeSize", ShapeSize, "TMOPThenOptimize", TMOPThenOptimize)

  RODIN_ORBIT_BENCH("Circle", CircleOrbit, "SquaredDistance", SquaredDistance)
  BENCHMARK(BM_LevelSetPipeline_CarryForwardOrbit)
    ->Name("LevelSetPipeline/CarryForward/Circle/SquaredDistance/24x20")
    ->Args({24, 20, static_cast<int>(ShapeCase::CircleOrbit),
      static_cast<int>(MetricCase::SquaredDistance)})
    ->Unit(benchmark::kMillisecond);
  RODIN_ORBIT_BENCH("Circle", CircleOrbit, "AreaDistortion", AreaDistortion)
  BENCHMARK(BM_LevelSetPipeline_CarryForwardOrbit)
    ->Name("LevelSetPipeline/CarryForward/Circle/AreaDistortion/24x20")
    ->Args({24, 20, static_cast<int>(ShapeCase::CircleOrbit),
      static_cast<int>(MetricCase::AreaDistortion)})
    ->Unit(benchmark::kMillisecond);
  RODIN_ORBIT_BENCH("Circle", CircleOrbit, "ShapeSize", ShapeSize)
  BENCHMARK(BM_LevelSetPipeline_CarryForwardOrbit)
    ->Name("LevelSetPipeline/CarryForward/Circle/ShapeSize/24x20")
    ->Args({24, 20, static_cast<int>(ShapeCase::CircleOrbit),
      static_cast<int>(MetricCase::ShapeSize)})
    ->Unit(benchmark::kMillisecond);
  RODIN_ORBIT_BENCH("Circle", CircleOrbit, "ShapeDistortion", ShapeDistortion)
  RODIN_ORBIT_BENCH("Square", SquareOrbit, "SquaredDistance", SquaredDistance)
  RODIN_ORBIT_BENCH("Square", SquareOrbit, "AreaDistortion", AreaDistortion)
  RODIN_ORBIT_BENCH("Square", SquareOrbit, "ShapeSize", ShapeSize)
  RODIN_ORBIT_BENCH("Triangle", TriangleOrbit, "SquaredDistance", SquaredDistance)
  RODIN_ORBIT_BENCH("Triangle", TriangleOrbit, "AreaDistortion", AreaDistortion)
  RODIN_ORBIT_BENCH("Triangle", TriangleOrbit, "ShapeSize", ShapeSize)
  RODIN_ORBIT_BENCH("CircleEnterLeave", CircleEnterLeave, "SquaredDistance", SquaredDistance)
  RODIN_ORBIT_BENCH("CircleFigureEight", CircleFigureEight, "SquaredDistance", SquaredDistance)
  RODIN_ORBIT_BENCH("SquareFigureEight", SquareFigureEight, "SquaredDistance", SquaredDistance)
  RODIN_ORBIT_BENCH("TriangleFigureEight", TriangleFigureEight, "SquaredDistance", SquaredDistance)
  RODIN_ORBIT_BENCH("WavyCircle", WavyCircleOrbit, "SquaredDistance", SquaredDistance)
  RODIN_ORBIT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "SquaredDistance", SquaredDistance)
  RODIN_ORBIT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "AreaDistortion", AreaDistortion)
  RODIN_ORBIT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "ShapeSize", ShapeSize)
  RODIN_ORBIT_BENCH("WavyCircleEnterLeave", WavyCircleEnterLeave, "SquaredDistance", SquaredDistance)

  RODIN_RESOLUTION_ORBIT_BENCH("Circle", CircleOrbit, "SquaredDistance", SquaredDistance, 10, 8)
  RODIN_RESOLUTION_ORBIT_BENCH("Circle", CircleOrbit, "SquaredDistance", SquaredDistance, 12, 8)
  RODIN_RESOLUTION_ORBIT_BENCH("Circle", CircleOrbit, "SquaredDistance", SquaredDistance, 16, 8)
  RODIN_RESOLUTION_ORBIT_BENCH("Circle", CircleOrbit, "SquaredDistance", SquaredDistance, 24, 8)
  RODIN_RESOLUTION_ORBIT_BENCH("Circle", CircleOrbit, "SquaredDistance", SquaredDistance, 32, 8)
  RODIN_RESOLUTION_ORBIT_BENCH("Circle", CircleOrbit, "ShapeSize", ShapeSize, 10, 8)
  RODIN_RESOLUTION_ORBIT_BENCH("Circle", CircleOrbit, "ShapeSize", ShapeSize, 12, 8)
  RODIN_RESOLUTION_ORBIT_BENCH("Circle", CircleOrbit, "ShapeSize", ShapeSize, 16, 8)
  RODIN_RESOLUTION_ORBIT_BENCH("Circle", CircleOrbit, "ShapeSize", ShapeSize, 24, 8)
  RODIN_RESOLUTION_ORBIT_BENCH("Circle", CircleOrbit, "ShapeSize", ShapeSize, 32, 8)
  RODIN_RESOLUTION_ORBIT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "SquaredDistance", SquaredDistance, 10, 8)
  RODIN_RESOLUTION_ORBIT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "SquaredDistance", SquaredDistance, 12, 8)
  RODIN_RESOLUTION_ORBIT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "SquaredDistance", SquaredDistance, 16, 8)
  RODIN_RESOLUTION_ORBIT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "SquaredDistance", SquaredDistance, 24, 8)
  RODIN_RESOLUTION_ORBIT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "SquaredDistance", SquaredDistance, 32, 8)
  RODIN_RESOLUTION_ORBIT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "ShapeSize", ShapeSize, 10, 8)
  RODIN_RESOLUTION_ORBIT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "ShapeSize", ShapeSize, 12, 8)
  RODIN_RESOLUTION_ORBIT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "ShapeSize", ShapeSize, 16, 8)
  RODIN_RESOLUTION_ORBIT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "ShapeSize", ShapeSize, 24, 8)
  RODIN_RESOLUTION_ORBIT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "ShapeSize", ShapeSize, 32, 8)

  RODIN_CURVED_FRESH_BENCH("Circle", CircleOrbit, "ShapeSize", ShapeSize, 10)
  RODIN_CURVED_FRESH_BENCH("Circle", CircleOrbit, "ShapeSize", ShapeSize, 16)
  RODIN_CURVED_FRESH_BENCH("Circle", CircleOrbit, "ShapeSize", ShapeSize, 24)
  RODIN_CURVED_FRESH_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "ShapeSize", ShapeSize, 10)
  RODIN_CURVED_FRESH_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "ShapeSize", ShapeSize, 16)
  RODIN_CURVED_FRESH_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "ShapeSize", ShapeSize, 24)
  RODIN_CURVED_FRESH_BENCH("Circle", CircleOrbit, "SquaredDistance", SquaredDistance, 10)
  RODIN_CURVED_FRESH_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "SquaredDistance", SquaredDistance, 10)

  RODIN_CURVED_ORBIT_BENCH("Circle", CircleOrbit, "ShapeSize", ShapeSize, 10, 8)
  RODIN_CURVED_ORBIT_BENCH("Circle", CircleOrbit, "ShapeSize", ShapeSize, 16, 8)
  RODIN_CURVED_ORBIT_BENCH("Circle", CircleOrbit, "ShapeSize", ShapeSize, 24, 8)
  RODIN_CURVED_ORBIT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "ShapeSize", ShapeSize, 10, 8)
  RODIN_CURVED_ORBIT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "ShapeSize", ShapeSize, 16, 8)
  RODIN_CURVED_ORBIT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "ShapeSize", ShapeSize, 24, 8)

  RODIN_STAGE_ORBIT_BENCH("Circle", CircleOrbit, "SquaredDistance", SquaredDistance, "TMOPOnly", TMOPOnly, 8)
  RODIN_STAGE_ORBIT_BENCH("Circle", CircleOrbit, "SquaredDistance", SquaredDistance, "OptimizeOnly", OptimizeOnly, 8)
  RODIN_STAGE_ORBIT_BENCH("Circle", CircleOrbit, "SquaredDistance", SquaredDistance, "TMOPThenOptimize", TMOPThenOptimize, 8)
  RODIN_STAGE_ORBIT_BENCH("Circle", CircleOrbit, "SquaredDistance", SquaredDistance, "OptimizeThenTMOP", OptimizeThenTMOP, 8)
  RODIN_STAGE_ORBIT_BENCH("Circle", CircleOrbit, "SquaredDistance", SquaredDistance, "OptimizeOnly", OptimizeOnly, 20)
  RODIN_STAGE_ORBIT_BENCH("Circle", CircleOrbit, "SquaredDistance", SquaredDistance, "TMOPThenOptimize", TMOPThenOptimize, 20)
  RODIN_STAGE_ORBIT_BENCH("Circle", CircleOrbit, "SquaredDistance", SquaredDistance, "OptimizeThenTMOP", OptimizeThenTMOP, 20)
  RODIN_STAGE_ORBIT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "SquaredDistance", SquaredDistance, "TMOPOnly", TMOPOnly, 8)
  RODIN_STAGE_ORBIT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "SquaredDistance", SquaredDistance, "OptimizeOnly", OptimizeOnly, 8)
  RODIN_STAGE_ORBIT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "SquaredDistance", SquaredDistance, "TMOPThenOptimize", TMOPThenOptimize, 8)
  RODIN_STAGE_ORBIT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "SquaredDistance", SquaredDistance, "OptimizeThenTMOP", OptimizeThenTMOP, 8)
  RODIN_STAGE_ORBIT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "SquaredDistance", SquaredDistance, "OptimizeOnly", OptimizeOnly, 20)
  RODIN_STAGE_ORBIT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "SquaredDistance", SquaredDistance, "TMOPThenOptimize", TMOPThenOptimize, 20)
  RODIN_STAGE_ORBIT_BENCH("WavyCircleFigureEight", WavyCircleFigureEight, "SquaredDistance", SquaredDistance, "OptimizeThenTMOP", OptimizeThenTMOP, 20)

  BENCHMARK(BM_LevelSetPipeline_CarryForwardOrbit)
    ->Name("LevelSetPipeline/CarryForward/Square/SquaredDistance/24x20")
    ->Args({24, 20, static_cast<int>(ShapeCase::SquareOrbit),
      static_cast<int>(MetricCase::SquaredDistance)})
    ->Unit(benchmark::kMillisecond);
  BENCHMARK(BM_LevelSetPipeline_CarryForwardOrbit)
    ->Name("LevelSetPipeline/CarryForward/Square/AreaDistortion/24x20")
    ->Args({24, 20, static_cast<int>(ShapeCase::SquareOrbit),
      static_cast<int>(MetricCase::AreaDistortion)})
    ->Unit(benchmark::kMillisecond);
  BENCHMARK(BM_LevelSetPipeline_CarryForwardOrbit)
    ->Name("LevelSetPipeline/CarryForward/Square/ShapeSize/24x20")
    ->Args({24, 20, static_cast<int>(ShapeCase::SquareOrbit),
      static_cast<int>(MetricCase::ShapeSize)})
    ->Unit(benchmark::kMillisecond);
  BENCHMARK(BM_LevelSetPipeline_CarryForwardOrbit)
    ->Name("LevelSetPipeline/CarryForward/Triangle/SquaredDistance/24x20")
    ->Args({24, 20, static_cast<int>(ShapeCase::TriangleOrbit),
      static_cast<int>(MetricCase::SquaredDistance)})
    ->Unit(benchmark::kMillisecond);
  BENCHMARK(BM_LevelSetPipeline_CarryForwardOrbit)
    ->Name("LevelSetPipeline/CarryForward/Triangle/AreaDistortion/24x20")
    ->Args({24, 20, static_cast<int>(ShapeCase::TriangleOrbit),
      static_cast<int>(MetricCase::AreaDistortion)})
    ->Unit(benchmark::kMillisecond);
  BENCHMARK(BM_LevelSetPipeline_CarryForwardOrbit)
    ->Name("LevelSetPipeline/CarryForward/Triangle/ShapeSize/24x20")
    ->Args({24, 20, static_cast<int>(ShapeCase::TriangleOrbit),
      static_cast<int>(MetricCase::ShapeSize)})
    ->Unit(benchmark::kMillisecond);
  BENCHMARK(BM_LevelSetPipeline_CarryForwardOrbit)
    ->Name("LevelSetPipeline/CarryForward/WavyCircleFigureEight/SquaredDistance/24x20")
    ->Args({24, 20, static_cast<int>(ShapeCase::WavyCircleFigureEight),
      static_cast<int>(MetricCase::SquaredDistance)})
    ->Unit(benchmark::kMillisecond);
  BENCHMARK(BM_LevelSetPipeline_CarryForwardOrbit)
    ->Name("LevelSetPipeline/CarryForward/WavyCircleFigureEight/AreaDistortion/24x20")
    ->Args({24, 20, static_cast<int>(ShapeCase::WavyCircleFigureEight),
      static_cast<int>(MetricCase::AreaDistortion)})
    ->Unit(benchmark::kMillisecond);
  BENCHMARK(BM_LevelSetPipeline_CarryForwardOrbit)
    ->Name("LevelSetPipeline/CarryForward/WavyCircleFigureEight/ShapeSize/24x20")
    ->Args({24, 20, static_cast<int>(ShapeCase::WavyCircleFigureEight),
      static_cast<int>(MetricCase::ShapeSize)})
    ->Unit(benchmark::kMillisecond);

  RODIN_WEIGHT_ORBIT_BENCH("Circle", CircleOrbit, 50)
  RODIN_WEIGHT_ORBIT_BENCH("Circle", CircleOrbit, 100)
  RODIN_WEIGHT_ORBIT_BENCH("Circle", CircleOrbit, 200)
  RODIN_WEIGHT_ORBIT_BENCH("Circle", CircleOrbit, 500)
  RODIN_WEIGHT_ORBIT_BENCH("Circle", CircleOrbit, 1000)

#undef RODIN_CUT_EXACT_BENCH
#undef RODIN_CUT_QUALITY_BENCH
#undef RODIN_TMOP_METRIC_BENCH
#undef RODIN_COARSEN_BENCH
#undef RODIN_ORBIT_BENCH
#undef RODIN_RESOLUTION_ORBIT_BENCH
#undef RODIN_CURVED_FRESH_BENCH
#undef RODIN_CURVED_ORBIT_BENCH
#undef RODIN_WEIGHT_FRESH_BENCH
#undef RODIN_WEIGHT_ORBIT_BENCH
#undef RODIN_STAGE_FRESH_BENCH
#undef RODIN_STAGE_ORBIT_BENCH
}
