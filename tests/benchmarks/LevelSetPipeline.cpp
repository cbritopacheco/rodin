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
#include <Rodin/Adaptation/TargetMatrixOptimization/IsoparametricGeometry.h>
#include <Rodin/Assembly/Default.h>
#include <Rodin/Geometry.h>
#include <Rodin/MMG.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
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

    enum class CurvedTargetCase : int
    {
      ProjectedInterface = 0,
      CurvedQuality005 = 1,
      CurvedQuality010 = 2,
      CurvedQuality020 = 3,
      ProjectedQuality005 = 4,
      ProjectedQuality010 = 5
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
      Index macroBoundaryCuts = 0;
      Index macroBoundaryCutRejectedByQuality = 0;
      Index sampledStencilFallbacks = 0;
      Index noCutQualityFallbacks = 0;
      Index graphOnlyFallbacks = 0;
      Index cutStencilOrder = 1;
      Index splits = 0;
      Index collapses = 0;
      Index featureCollapses = 0;
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
      Index curvedOverlapSamples = 0;
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
      Real tmopTargetMaxMetric = 0;
      Real tmopTargetMinDetT = 0;
      Index tmopTargetInvalidSamples = 0;
      Index tmopIterations = 0;
      Index inverted = 0;
      Index tmopFailures = 0;
      Index reclassifiedCells = 0;
    };

    struct TopologyReport
    {
      Index cellCount = 0;
      Index interfaceEdgeCount = 0;
      Index uncutCellCount = 0;
      Index snappedCrossingCount = 0;
      Index pathologicalCutCount = 0;
      Index macroBoundaryCutCount = 0;
      Index macroBoundaryCutRejectedByQualityCount = 0;
      Index sampledStencilFallbackCount = 0;
      Index noCutQualityFallbackCount = 0;
      Index graphOnlyFallbackCount = 0;
      Index referenceStencilOrder = 0;
      Real minOutputCellQuality = 0;
      Real maxInterfaceDeviation = 0;
    };

    struct TopologyResult
    {
      LocalMesh mesh;
      TopologyReport report;
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
      Real macroBoundaryCuts = 0;
      Real macroBoundaryCutRejectedByQuality = 0;
      Real sampledStencilFallbacks = 0;
      Real noCutQualityFallbacks = 0;
      Real graphOnlyFallbacks = 0;
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
      Real curvedOverlapSamples = 0;
      Real bestCurvedOverlapSamples = std::numeric_limits<Real>::infinity();
      Real worstCurvedOverlapSamples = 0;
      Real finalCurvedOverlapSamples = 0;
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
      Real tmopTargetMaxMetric = 0;
      Real maxTMOPTargetMetric = -std::numeric_limits<Real>::infinity();
      Real finalTMOPTargetMetric = 0;
      Real tmopTargetMinDetT = std::numeric_limits<Real>::infinity();
      Real finalTMOPTargetMinDetT = 0;
      Real tmopTargetInvalidSamples = 0;
      Real worstTMOPTargetInvalidSamples = 0;
      Real finalTMOPTargetInvalidSamples = 0;
      Real tmopIterations = 0;
      Real inverted = 0;
      Real tmopFailures = 0;
      Real featureCollapses = 0;
      Real reclassifiedCells = 0;
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
        macroBoundaryCuts += static_cast<Real>(counters.macroBoundaryCuts);
        macroBoundaryCutRejectedByQuality += static_cast<Real>(
            counters.macroBoundaryCutRejectedByQuality);
        sampledStencilFallbacks += static_cast<Real>(
            counters.sampledStencilFallbacks);
        noCutQualityFallbacks += static_cast<Real>(
            counters.noCutQualityFallbacks);
        graphOnlyFallbacks += static_cast<Real>(counters.graphOnlyFallbacks);
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
        curvedOverlapSamples +=
          static_cast<Real>(counters.curvedOverlapSamples);
        bestCurvedOverlapSamples = std::min(
            bestCurvedOverlapSamples,
            static_cast<Real>(counters.curvedOverlapSamples));
        worstCurvedOverlapSamples = std::max(
            worstCurvedOverlapSamples,
            static_cast<Real>(counters.curvedOverlapSamples));
        finalCurvedOverlapSamples =
          static_cast<Real>(counters.curvedOverlapSamples);
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
        tmopTargetMaxMetric += counters.tmopTargetMaxMetric;
        maxTMOPTargetMetric =
          std::max(maxTMOPTargetMetric, counters.tmopTargetMaxMetric);
        finalTMOPTargetMetric = counters.tmopTargetMaxMetric;
        if (std::isfinite(counters.tmopTargetMinDetT))
          tmopTargetMinDetT =
            std::min(tmopTargetMinDetT, counters.tmopTargetMinDetT);
        finalTMOPTargetMinDetT = counters.tmopTargetMinDetT;
        tmopTargetInvalidSamples +=
          static_cast<Real>(counters.tmopTargetInvalidSamples);
        worstTMOPTargetInvalidSamples = std::max(
            worstTMOPTargetInvalidSamples,
            static_cast<Real>(counters.tmopTargetInvalidSamples));
        finalTMOPTargetInvalidSamples =
          static_cast<Real>(counters.tmopTargetInvalidSamples);
        tmopIterations += static_cast<Real>(counters.tmopIterations);
        inverted += static_cast<Real>(counters.inverted);
        tmopFailures += static_cast<Real>(counters.tmopFailures);
        featureCollapses += static_cast<Real>(counters.featureCollapses);
        reclassifiedCells += static_cast<Real>(counters.reclassifiedCells);
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
        st.counters["avg_cut_macro"] =
          benchmark::Counter(macroBoundaryCuts * inv);
        st.counters["avg_cut_macro_rejected"] =
          benchmark::Counter(macroBoundaryCutRejectedByQuality * inv);
        st.counters["avg_cut_sampled_fallback"] =
          benchmark::Counter(sampledStencilFallbacks * inv);
        st.counters["avg_cut_nocut_fallback"] =
          benchmark::Counter(noCutQualityFallbacks * inv);
        st.counters["avg_cut_graph_only_fallback"] =
          benchmark::Counter(graphOnlyFallbacks * inv);
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
        if (hasP2Diagnostics)
        {
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
        st.counters["avg_curved_overlap"] =
          benchmark::Counter(curvedOverlapSamples * inv);
        st.counters["final_curved_overlap"] =
          benchmark::Counter(finalCurvedOverlapSamples);
        st.counters["best_curved_overlap"] =
          benchmark::Counter(std::isfinite(bestCurvedOverlapSamples) ? bestCurvedOverlapSamples : Real(0));
        st.counters["worst_curved_overlap"] =
          benchmark::Counter(worstCurvedOverlapSamples);

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
        if (hasP2Diagnostics)
        {
          st.counters["avg_tmop_target_max_mu"] =
            benchmark::Counter(tmopTargetMaxMetric * inv);
          st.counters["final_tmop_target_max_mu"] =
            benchmark::Counter(finalTMOPTargetMetric);
          st.counters["worst_tmop_target_max_mu"] =
            benchmark::Counter(std::isfinite(maxTMOPTargetMetric)
                ? maxTMOPTargetMetric : Real(0));
          st.counters["min_tmop_target_detT"] =
            benchmark::Counter(std::isfinite(tmopTargetMinDetT)
                ? tmopTargetMinDetT : Real(0));
          st.counters["final_tmop_target_detT"] =
            benchmark::Counter(finalTMOPTargetMinDetT);
          st.counters["avg_tmop_target_invalid"] =
            benchmark::Counter(tmopTargetInvalidSamples * inv);
          st.counters["final_tmop_target_invalid"] =
            benchmark::Counter(finalTMOPTargetInvalidSamples);
          st.counters["worst_tmop_target_invalid"] =
            benchmark::Counter(worstTMOPTargetInvalidSamples);
        }
        st.counters["avg_tmop_iterations"] =
          benchmark::Counter(tmopIterations * inv);
        st.counters["avg_inverted"] = benchmark::Counter(inverted * inv);
        st.counters["avg_tmop_failed"] =
          benchmark::Counter(tmopFailures * inv);
        st.counters["avg_feature_collapses"] =
          benchmark::Counter(featureCollapses * inv);
        st.counters["avg_reclassified_cells"] =
          benchmark::Counter(reclassifiedCells * inv);
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
      Real targetMaxMetric = 0;
      Real targetMinDetT = 0;
      Index targetInvalidSamples = 0;
      bool hasCurvedMetrics = false;
      Real curvedFitRms = 0;
      Real curvedFitMax = 0;
      Real curvedMinQuality = 0;
      Real curvedMinJacobian = 0;
      Index curvedInvalidSamples = 0;
      Index curvedOverlapSamples = 0;
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
      Index overlapSamples = 0;  ///< Curved interface samples in non-edge-adjacent neighbours; canary for cell-cell overlap (global-injectivity failure missed by per-cell det/quality).
    };

    inline bool pointInLinearTriangle(
        const Math::SpatialPoint& p,
        const Math::SpatialPoint& a,
        const Math::SpatialPoint& b,
        const Math::SpatialPoint& c)
    {
      const Real d = (b[0] - a[0]) * (c[1] - a[1])
                   - (b[1] - a[1]) * (c[0] - a[0]);
      if (std::abs(d) <= Real(1e-30))
        return false;
      const Real inv = Real(1) / d;
      const Real wb = ((c[1] - a[1]) * (p[0] - a[0])
                    - (c[0] - a[0]) * (p[1] - a[1])) * inv;
      const Real wc = ((b[0] - a[0]) * (p[1] - a[1])
                    - (b[1] - a[1]) * (p[0] - a[0])) * inv;
      const Real wa = Real(1) - wb - wc;
      const Real eps = Real(1e-9);
      return wa > eps && wb > eps && wc > eps;
    }

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
        const auto& qf = QF::PolytopeQuadratureFormula::get(
            4, Polytope::Type::Triangle);
        for (size_t q = 0; q < qf.getSize(); ++q)
        {
          Math::SpatialMatrix<Real> J;
          trans.jacobian(J, qf.getPoint(q));
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

      // Cell-cell overlap canary: sample every Interface edge's curved
      // transformation at 5 interior points and check membership in cells
      // that share a vertex with the edge but NOT the edge itself. A hit
      // means the curved interface intrudes into a non-incident neighbour
      // -- global-injectivity failure missed by per-cell det/quality.
      const_cast<LocalMesh&>(mesh).getConnectivity().compute(0, 2);
      const_cast<LocalMesh&>(mesh).getConnectivity().compute(1, 2);
      static const std::array<Real, 5> sOv = {{
        Real(0.1), Real(0.25), Real(0.5), Real(0.75), Real(0.9) }};
      for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
      {
        const auto attr = mesh.getAttribute(1, e);
        if (!attr || *attr != Interface)
          continue;
        const auto& ep = conn.getPolytope(1, e);
        std::array<Index, 2> edgeAdj{
          std::numeric_limits<Index>::max(),
          std::numeric_limits<Index>::max() };
        {
          size_t k = 0;
          for (Index ci : conn.getIncidence({ 1, 2 }, e))
            if (k < 2) edgeAdj[k++] = ci;
        }
        std::vector<Index> nbrs;
        for (Index vEnd : { ep(0), ep(1) })
          for (Index ci : conn.getIncidence({ 0, 2 }, vEnd))
            if (ci != edgeAdj[0] && ci != edgeAdj[1])
              nbrs.push_back(ci);
        std::sort(nbrs.begin(), nbrs.end());
        nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
        if (nbrs.empty())
          continue;
        const auto edgeIt = mesh.getPolytope(1, e);
        for (Real s : sOv)
        {
          Math::SpatialPoint sample;
          edgeIt->getTransformation().transform(
              sample, Math::SpatialPoint{ s });
          for (Index ci : nbrs)
          {
            const auto& cell = conn.getPolytope(2, ci);
            const auto Va = mesh.getVertexCoordinates(cell(0));
            const auto Vb = mesh.getVertexCoordinates(cell(1));
            const auto Vc = mesh.getVertexCoordinates(cell(2));
            if (pointInLinearTriangle(sample, Va, Vb, Vc))
            {
              ++stats.overlapSamples;
              break;
            }
          }
        }
      }
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

    Index markUncutAnalyticInterfaceEdges(
        LocalMesh& mesh,
        ShapeCase shape,
        Real t)
    {
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      Index count = 0;
      const auto& conn = mesh.getConnectivity();
      for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
      {
        const auto attr = mesh.getAttribute(1, e);
        if (attr && *attr == Boundary)
          continue;
        const auto& edge = conn.getPolytope(1, e);
        const Real f0 =
          levelSetValue(shape, mesh.getVertexCoordinates(edge(0)), t);
        const Real f1 =
          levelSetValue(shape, mesh.getVertexCoordinates(edge(1)), t);
        if (f0 == Real(0) || f1 == Real(0) || f0 * f1 < Real(0))
        {
          mesh.setAttribute({1, e}, Interface);
          ++count;
        }
      }
      return count;
    }

    TopologyResult makeUncutLevelSetProxy(
        LocalMesh mesh,
        ShapeCase shape,
        Real t)
    {
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      annotateBoundary(mesh);
      const Index interfaceEdges =
        markUncutAnalyticInterfaceEdges(mesh, shape, t);

      TopologyResult result;
      result.mesh = std::move(mesh);
      const auto stats = meshStats(result.mesh);
      result.report.cellCount = result.mesh.getCellCount();
      result.report.minOutputCellQuality = stats.minQuality;
      result.report.uncutCellCount = result.mesh.getCellCount();
      result.report.referenceStencilOrder = 0;
      const auto& conn = result.mesh.getConnectivity();
      for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
      {
        const auto attr = result.mesh.getAttribute(1, e);
        if (!attr || *attr != Interface)
          continue;
        ++result.report.interfaceEdgeCount;
        const auto& edge = conn.getPolytope(1, e);
        for (Index v : {edge(0), edge(1)})
          result.report.maxInterfaceDeviation =
            std::max(result.report.maxInterfaceDeviation,
                std::abs(levelSetValue(
                    shape, result.mesh.getVertexCoordinates(v), t)));
      }
      result.report.macroBoundaryCutCount = interfaceEdges;
      return result;
    }

    void resetCellAttributes(LocalMesh& mesh, Attribute attr)
    {
      for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
        mesh.setAttribute({ 2, c }, attr);
    }

    TopologyResult makeExternalCutResult(
        LocalMesh mesh,
        ShapeCase shape,
        Real t)
    {
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
      annotateBoundary(mesh);

      TopologyResult result;
      result.mesh = std::move(mesh);
      const auto stats = meshStats(result.mesh);
      result.report.cellCount = result.mesh.getCellCount();
      result.report.minOutputCellQuality = stats.minQuality;
      result.report.uncutCellCount = result.mesh.getCellCount();
      result.report.referenceStencilOrder = 0;
      const auto& conn = result.mesh.getConnectivity();
      for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
      {
        const auto attr = result.mesh.getAttribute(1, e);
        if (attr && *attr == Interface)
        {
          ++result.report.interfaceEdgeCount;
          const auto& edge = conn.getPolytope(1, e);
          for (Index v : { edge(0), edge(1) })
            result.report.maxInterfaceDeviation =
              std::max(result.report.maxInterfaceDeviation,
                  std::abs(levelSetValue(
                      shape, result.mesh.getVertexCoordinates(v), t)));
        }
      }
      return result;
    }

    bool localVertexMoveValid(
        const LocalMesh& mesh,
        Index vertex,
        const Math::SpatialPoint& candidate,
        Real minQuality = Real(1e-5))
    {
      const auto& conn = mesh.getConnectivity();
      for (Index ci : conn.getIncidence({ 0, 2 }, vertex))
      {
        const auto& cell = conn.getPolytope(2, ci);
        std::array<Math::SpatialPoint, 3> x{{
          mesh.getVertexCoordinates(cell(0)),
          mesh.getVertexCoordinates(cell(1)),
          mesh.getVertexCoordinates(cell(2)) }};
        for (size_t k = 0; k < 3; ++k)
          if (cell(k) == vertex)
            x[k] = candidate;
        if (!(signedTriangleArea(x[0], x[1], x[2]) > Real(1e-14)))
          return false;
        if (triangleQuality(x[0], x[1], x[2]) < minQuality)
          return false;
      }
      return true;
    }

    void snapInterfaceVerticesToAnalytic(
        LocalMesh& mesh,
        ShapeCase shape,
        Real t,
        Real maxMove)
    {
      mesh.getConnectivity().compute(0, 2);
      mesh.getConnectivity().compute(1, 0);
      std::vector<char> mask(mesh.getVertexCount(), 0);
      const auto& conn = mesh.getConnectivity();
      for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
      {
        const auto attr = mesh.getAttribute(1, e);
        if (!attr || *attr != Interface)
          continue;
        const auto& edge = conn.getPolytope(1, e);
        mask[edge(0)] = 1;
        mask[edge(1)] = 1;
      }
      for (Index v = 0; v < static_cast<Index>(mesh.getVertexCount()); ++v)
      {
        if (v >= static_cast<Index>(mask.size()) || !mask[v])
          continue;
        const auto x = mesh.getVertexCoordinates(v);
        auto p = projectToLevelSet(shape, t, x);
        const Real d = (p - x).norm();
        if (d > maxMove && d > Real(0))
          p = x + (maxMove / d) * (p - x);
        if (localVertexMoveValid(mesh, v, p))
          mesh.setVertexCoordinates(v, p);
      }
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 0);
    }

    TopologyResult cutShapeWithMMG(
        LocalMesh& background,
        ShapeCase shape,
        Real t,
        size_t resolution)
    {
      const Real h = Real(1) / static_cast<Real>(resolution - 1);
      LocalMesh base(background);
      base.getConnectivity().compute(2, 1);
      base.getConnectivity().compute(1, 0);
      resetCellAttributes(base, Negative);
      annotateBoundary(base);

      MMG::Mesh mmgMesh(std::move(base));
      P1<Real, LocalMesh> phiSpace(mmgMesh);
      MMG::RealGridFunction phi(phiSpace);
      for (Index v = 0; v < mmgMesh.getVertexCount(); ++v)
        phi[v] = levelSetValue(shape, mmgMesh.getVertexCoordinates(v), t);

      MMG::Mesh fitted = MMG::LevelSetDiscretizer()
        .split(Negative, { Negative, Positive })
        .setBoundaryReference(Interface)
        .setHMin(Real(0.35) * h)
        .setHMax(Real(1.35) * h)
        .setHausdorff(Real(0.05) * h)
        .setGradation(Real(1.25))
        .setAngleDetection(false)
        .discretize(phi);
      snapInterfaceVerticesToAnalytic(
          fitted, shape, t, Real(0.35) * h);

      return makeExternalCutResult(std::move(fitted), shape, t);
    }

    void optimizeMMGTopology(
        LocalMesh& mesh,
        size_t resolution)
    {
      const Real h = Real(1) / static_cast<Real>(resolution - 1);
      MMG::Mesh mmgMesh(std::move(mesh));
      MMG::Optimizer()
        .setHMin(Real(0.35) * h)
        .setHMax(Real(1.35) * h)
        .setHausdorff(Real(0.05) * h)
        .setGradation(Real(1.25))
        .setAngleDetection(false)
        .optimize(mmgMesh);
      mmgMesh.getConnectivity().compute(2, 1);
      mmgMesh.getConnectivity().compute(1, 0);
      annotateBoundary(mmgMesh);
      mesh = LocalMesh(std::move(mmgMesh));
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

    Index reclassifyCellsByAnalyticPhi(
        LocalMesh& mesh,
        ShapeCase shape,
        Real t)
    {
      Index changed = 0;
      const auto& conn = mesh.getConnectivity();
      for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
      {
        if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
          continue;
        const auto cellIt = mesh.getPolytope(2, c);
        Math::SpatialPoint x;
        cellIt->getTransformation().transform(
            x, Math::SpatialPoint{Real(1) / Real(3), Real(1) / Real(3)});
        const Attribute next =
          levelSetValue(shape, x, t) <= Real(0) ? Negative : Positive;
        const auto old = mesh.getAttribute(2, c);
        if (!old || *old != next)
        {
          mesh.setAttribute({2, c}, next);
          ++changed;
        }
      }
      return changed;
    }

    template <class Metric>
    bool solveCurvedTMOPWithMetric(
        LocalMesh& mesh,
        ShapeCase shape,
        Real t,
        Metric metric,
        Real qualityWeight,
        Index maxIterations,
        CurvedTargetCase targetCase = CurvedTargetCase::ProjectedInterface,
        TMOPSolveStats* tmopStats = nullptr)
    {
      auto phiValue = [shape, t](const Math::SpatialPoint& x)
      {
        return levelSetValue(shape, x, t);
      };
      auto phiGradient = [shape, t](const Math::SpatialPoint& x)
      {
        return levelSetGradient(shape, x, t);
      };
      try
      {
        Variational::RealH1Element<2> fe(Polytope::Type::Triangle);
        upgradeTransformations(mesh, fe, Interface);  // affine P2 lift

        auto projectToInterface = [shape, t](const Math::SpatialPoint& x)
        {
          return projectToLevelSet(shape, t, x);
        };
        auto solveWithTarget = [&](const auto& target)
        {
          LocalMesh beforeSolve(mesh);
          QualityTerm quality(metric, target, qualityWeight);
          quality.setQuadratureOrder(4);
          DeviationTerm deviation(Real(0.25));
          AnalyticLevelSetFitTerm fit(
              phiValue, phiGradient, Optional<Attribute>(Interface), Real(2));
          AnalyticLevelSetFitTerm bfit(
              boxBoundaryValue, boxBoundaryGradient,
              Optional<Attribute>(Boundary), Real(1));

          VectorH1<2, LocalMesh> space(
              std::integral_constant<size_t, 2>{}, mesh, 2);
          GridFunction u(space);
          u.getData().setZero();
          TrialFunction du(space);
          TestFunction  v (space);

          auto makeResidual = [&] {
            return quality.residual(u, v) + deviation.residual(u, v)
                 + fit.residual(u, v)     + bfit.residual(u, v);
          };
          auto makeTangent = [&] {
            return quality.tangent(u, du, v) + deviation.tangent(du, v)
                 + fit.tangent(u, du, v)     + bfit.tangent(u, du, v);
          };
          auto energy = [&] {
            return quality.energy(u) + deviation.energy(u)
                 + fit.energy(u)     + bfit.energy(u);
          };

          IsoparametricTMOPParameters params;
          params.maxIterations = maxIterations;
          const auto report = solveIsoparametricTMOP(
              mesh, fe, u, du, v, makeResidual, makeTangent, energy,
              Interface, params);
          if (tmopStats)
          {
            tmopStats->iterations += report.iterations;
            const auto targetStats =
              targetQualityMetrics(mesh, metric, target);
            tmopStats->targetMaxMetric = targetStats.maxMetric;
            tmopStats->targetMinDetT = targetStats.minDetT;
            tmopStats->targetInvalidSamples = targetStats.invalidSamples;
            const auto curvedStats = curvedInterfaceStats(mesh, shape, t);
            tmopStats->curvedFitRms = curvedStats.fitRms;
            tmopStats->curvedFitMax = curvedStats.fitMax;
            tmopStats->curvedMinQuality = curvedStats.minQuality;
            tmopStats->curvedMinJacobian = curvedStats.minJacobian;
            tmopStats->curvedInvalidSamples =
              curvedStats.invalidJacobianSamples;
            tmopStats->curvedOverlapSamples = curvedStats.overlapSamples;
            tmopStats->hasCurvedMetrics = true;
          }
          syncLinearBackbone(mesh, fe);
          demoteTransformations(mesh);

          Index invalid = 0;
          Index overlap = 0;
          Real minJac = 0;
          Real minQ = 0;
          if (tmopStats && tmopStats->hasCurvedMetrics)
          {
            invalid = tmopStats->curvedInvalidSamples;
            overlap = tmopStats->curvedOverlapSamples;
            minJac = tmopStats->curvedMinJacobian;
            minQ = tmopStats->curvedMinQuality;
          }
          else
          {
            const auto curved = curvedInterfaceStats(mesh, shape, t);
            invalid = curved.invalidJacobianSamples;
            overlap = curved.overlapSamples;
            minJac = curved.minJacobian;
            minQ = curved.minQuality;
          }
          const auto finalStats = meshStats(mesh);
          const bool valid = invalid == 0
            && overlap == 0
            && minJac > Real(0)
            && std::isfinite(minJac)
            && minQ > Real(0)
            && std::isfinite(minQ)
            && finalStats.inverted == 0
            && std::isfinite(finalStats.coverage)
            && finalStats.coverage > Real(0.8)
            && finalStats.coverage < Real(1.2)
            && (report.converged || report.iterations > 0);
          if (!valid)
          {
            mesh = std::move(beforeSolve);
            mesh.getConnectivity().compute(2, 1);
            mesh.getConnectivity().compute(1, 0);
            if (tmopStats)
              tmopStats->hasCurvedMetrics = false;
            return false;
          }
          return true;
        };

        switch (targetCase)
        {
          case CurvedTargetCase::ProjectedInterface:
          {
            ProjectedInterfaceTargetJacobian target(
                mesh, Interface, projectToInterface);
            return solveWithTarget(target);
          }
          case CurvedTargetCase::CurvedQuality005:
            return solveWithTarget(CurvedQualityTargetJacobian(mesh, Real(0.05)));
          case CurvedTargetCase::CurvedQuality010:
            return solveWithTarget(CurvedQualityTargetJacobian(mesh, Real(0.10)));
          case CurvedTargetCase::CurvedQuality020:
            return solveWithTarget(CurvedQualityTargetJacobian(mesh, Real(0.20)));
          case CurvedTargetCase::ProjectedQuality005:
          {
            ProjectedQualityTargetJacobian target(
                mesh, Interface, projectToInterface, Real(0.05));
            return solveWithTarget(target);
          }
          case CurvedTargetCase::ProjectedQuality010:
          {
            ProjectedQualityTargetJacobian target(
                mesh, Interface, projectToInterface, Real(0.10));
            return solveWithTarget(target);
          }
        }
        return false;
      }
      catch (...)
      {
        return false;
      }
    }

    bool solveCurvedTMOP(
        LocalMesh& mesh,
        ShapeCase shape,
        Real t,
        MetricCase metricCase = MetricCase::ShapeSize,
        CurvedTargetCase targetCase = CurvedTargetCase::ProjectedInterface,
        TMOPSolveStats* stats = nullptr)
    {
      switch (metricCase)
      {
        case MetricCase::SquaredDistance:
          return solveCurvedTMOPWithMetric(
              mesh, shape, t, SquaredDistanceMetric{}, 0.04, 4,
              targetCase, stats);
        case MetricCase::AreaDistortion:
          return solveCurvedTMOPWithMetric(
              mesh, shape, t, AreaDistortionMetric{}, 0.04, 4,
              targetCase, stats);
        case MetricCase::ShapeSize:
          return solveCurvedTMOPWithMetric(
              mesh, shape, t, ShapeSizeBlendMetric(Real(0.5)), 0.03, 4,
              targetCase, stats);
        case MetricCase::ShapeDistortion:
          return solveCurvedTMOPWithMetric(
              mesh, shape, t, ShapeDistortionMetric{}, 0.02, 1,
              targetCase, stats);
      }
      return false;
    }

    StageCounters collectStageCounters(
        const TopologyResult& cut,
        const LocalMesh& finalMesh,
        ShapeCase shape,
        bool tmopSucceeded,
        Real t,
        bool hasP2Diagnostics = false)
    {
      StageCounters counters;
      const auto stats = meshStats(finalMesh);
      const auto fit = interfaceFit(finalMesh, shape, t);
      counters.cutCells = cut.mesh.getCellCount();
      counters.finalCells = finalMesh.getCellCount();
      counters.interfaceEdges = cut.report.interfaceEdgeCount;
      counters.uncutCells = cut.report.uncutCellCount;
      counters.snappedCrossings = cut.report.snappedCrossingCount;
      counters.pathologicalCuts = cut.report.pathologicalCutCount;
      counters.macroBoundaryCuts = cut.report.macroBoundaryCutCount;
      counters.macroBoundaryCutRejectedByQuality =
        cut.report.macroBoundaryCutRejectedByQualityCount;
      counters.sampledStencilFallbacks = cut.report.sampledStencilFallbackCount;
      counters.noCutQualityFallbacks = cut.report.noCutQualityFallbackCount;
      counters.graphOnlyFallbacks = cut.report.graphOnlyFallbackCount;
      counters.cutStencilOrder = cut.report.referenceStencilOrder;
      counters.cutMinQuality = cut.report.minOutputCellQuality;
      counters.finalMinQuality = stats.minQuality;
      counters.fitRms = fit.first;
      counters.fitMax = fit.second;
      counters.hasP2Diagnostics = hasP2Diagnostics;
      if (hasP2Diagnostics)
      {
        const auto curved = curvedInterfaceStats(finalMesh, shape, t);
        counters.curvedOverlapSamples = curved.overlapSamples;
        counters.curvedFitRms = curved.fitRms;
        counters.curvedFitMax = curved.fitMax;
        counters.curvedMinQuality = curved.minQuality;
        counters.curvedMinJacobian = curved.minJacobian;
        counters.curvedInvalidSamples = curved.invalidJacobianSamples;
      }
      counters.coverage = stats.coverage;
      counters.signedArea = stats.signedArea;
      counters.maxInterfaceDeviation = cut.report.maxInterfaceDeviation;
      counters.inverted = stats.inverted;
      counters.tmopFailures = tmopSucceeded ? 0 : 1;

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
      st.counters["cut_stencil_order"] =
        benchmark::Counter(counters.cutStencilOrder);
      st.counters["cut_macro"] =
        benchmark::Counter(counters.macroBoundaryCuts);
      st.counters["cut_macro_rejected"] =
        benchmark::Counter(counters.macroBoundaryCutRejectedByQuality);
      st.counters["cut_sampled_fallback"] =
        benchmark::Counter(counters.sampledStencilFallbacks);
      st.counters["cut_nocut_fallback"] =
        benchmark::Counter(counters.noCutQualityFallbacks);
      st.counters["cut_graph_only_fallback"] =
        benchmark::Counter(counters.graphOnlyFallbacks);
      st.counters["feature_collapses"] =
        benchmark::Counter(counters.featureCollapses);
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
      if (counters.hasP2Diagnostics)
      {
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
        st.counters["tmop_target_max_mu"] =
          benchmark::Counter(counters.tmopTargetMaxMetric);
        st.counters["avg_tmop_target_max_mu"] =
          benchmark::Counter(counters.tmopTargetMaxMetric);
        st.counters["final_tmop_target_max_mu"] =
          benchmark::Counter(counters.tmopTargetMaxMetric);
        st.counters["min_tmop_target_detT"] =
          benchmark::Counter(counters.tmopTargetMinDetT);
        st.counters["final_tmop_target_detT"] =
          benchmark::Counter(counters.tmopTargetMinDetT);
        st.counters["tmop_target_invalid"] =
          benchmark::Counter(counters.tmopTargetInvalidSamples);
        st.counters["avg_tmop_target_invalid"] =
          benchmark::Counter(counters.tmopTargetInvalidSamples);
        st.counters["final_tmop_target_invalid"] =
          benchmark::Counter(counters.tmopTargetInvalidSamples);
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
      st.counters["reclassified_cells"] =
        benchmark::Counter(counters.reclassifiedCells);
      st.counters["splits"] = benchmark::Counter(counters.splits);
      st.counters["collapses"] = benchmark::Counter(counters.collapses);
      st.counters["swaps"] = benchmark::Counter(counters.swaps);
      st.counters["smooths"] = benchmark::Counter(counters.smooths);
      st.counters["feature_smooths"] =
        benchmark::Counter(counters.featureSmooths);
    }

    StageCounters setCounters(
        benchmark::State& st,
        const TopologyResult& cut,
        const LocalMesh& finalMesh,
        ShapeCase shape,
        bool tmopSucceeded,
        Real t,
        Real cutSeconds = 0,
        Real optimizeSeconds = 0,
        Real tmopSeconds = 0,
        const TMOPSolveStats* tmopStats = nullptr,
        bool hasP2Diagnostics = false,
        Index reclassifiedCells = 0)
    {
      auto counters = collectStageCounters(
          cut, finalMesh, shape, tmopSucceeded, t,
          hasP2Diagnostics);
      counters.cutSeconds = cutSeconds;
      counters.optimizeSeconds = optimizeSeconds;
      counters.tmopSeconds = tmopSeconds;
      counters.reclassifiedCells = reclassifiedCells;
      if (tmopStats)
      {
        counters.tmopAssemblySeconds = tmopStats->assemblySeconds;
        counters.tmopSolveSeconds = tmopStats->solveSeconds;
        counters.tmopMeritSeconds = tmopStats->meritSeconds;
        counters.tmopTargetMaxMetric = tmopStats->targetMaxMetric;
        counters.tmopTargetMinDetT = tmopStats->targetMinDetT;
        counters.tmopTargetInvalidSamples = tmopStats->targetInvalidSamples;
        counters.tmopIterations = tmopStats->iterations;
        if (tmopStats->hasCurvedMetrics)
        {
          counters.hasP2Diagnostics = true;
          counters.curvedFitRms = tmopStats->curvedFitRms;
          counters.curvedFitMax = tmopStats->curvedFitMax;
          counters.curvedMinQuality = tmopStats->curvedMinQuality;
          counters.curvedMinJacobian = tmopStats->curvedMinJacobian;
          counters.curvedInvalidSamples = tmopStats->curvedInvalidSamples;
          counters.curvedOverlapSamples = tmopStats->curvedOverlapSamples;
        }
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

    static void BM_LevelSetPipeline_CurvedTargetCarryForwardOrbit(
        benchmark::State& st)
    {
      const auto resolution = static_cast<size_t>(st.range(0));
      const auto steps = static_cast<Index>(st.range(1));
      const auto shape = static_cast<ShapeCase>(st.range(2));
      const auto metric = static_cast<MetricCase>(st.range(3));
      const auto target = static_cast<CurvedTargetCase>(st.range(4));
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
          auto topology =
            makeUncutLevelSetProxy(LocalMesh(background), shape, t);
          const Real cutSeconds = elapsedSeconds(cutStart);
          auto curved = topology.mesh;
          curved.getConnectivity().compute(2, 1);
          curved.getConnectivity().compute(1, 0);

          TMOPSolveStats tmopStats;
          const auto tmopStart = BenchClock::now();
          const bool tmopOk =
            solveCurvedTMOP(curved, shape, t, metric, target, &tmopStats);
          const Real tmopSeconds = elapsedSeconds(tmopStart);
          const Index reclassified =
            reclassifyCellsByAnalyticPhi(curved, shape, t);
          if (!tmopOk)
            tmopFailures++;

          counters = setCounters(st, topology, curved, shape, tmopOk, t,
              cutSeconds, 0, tmopSeconds, &tmopStats, true,
              reclassified);
          averages.add(counters);

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
      st.counters["curved_target_case"] =
        benchmark::Counter(static_cast<int>(target));
      st.counters["production_order"] = benchmark::Counter(1);
      st.counters["tmop_failed"] = benchmark::Counter(counters.tmopFailures);
      st.SetItemsProcessed(st.iterations() * counters.finalCells);
    }

    static void BM_LevelSetPipeline_CurvedMMGCarryForwardOrbit(
        benchmark::State& st)
    {
      const auto resolution = static_cast<size_t>(st.range(0));
      const auto steps = static_cast<Index>(st.range(1));
      const auto shape = static_cast<ShapeCase>(st.range(2));
      const auto metric = static_cast<MetricCase>(st.range(3));
      const auto target = static_cast<CurvedTargetCase>(st.range(4));
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
          auto topology = cutShapeWithMMG(background, shape, t, resolution);
          const Real cutSeconds = elapsedSeconds(cutStart);
          auto optimized = topology.mesh;

          const auto optStart = BenchClock::now();
          optimizeMMGTopology(optimized, resolution);
          const Real optimizeSeconds = elapsedSeconds(optStart);

          optimized.getConnectivity().compute(2, 1);
          optimized.getConnectivity().compute(1, 0);
          auto curved = optimized;

          TMOPSolveStats tmopStats;
          const auto tmopStart = BenchClock::now();
          const bool tmopOk =
            solveCurvedTMOP(curved, shape, t, metric, target, &tmopStats);
          const Real tmopSeconds = elapsedSeconds(tmopStart);
          const Index reclassified =
            reclassifyCellsByAnalyticPhi(curved, shape, t);
          if (!tmopOk)
            tmopFailures++;

          counters = setCounters(st, topology, curved, shape, tmopOk, t,
              cutSeconds, optimizeSeconds, tmopSeconds, &tmopStats, true,
              reclassified);
          averages.add(counters);

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
      st.counters["mmg_topology"] = benchmark::Counter(1);
      st.counters["curved_target_case"] =
        benchmark::Counter(static_cast<int>(target));
      st.counters["production_order"] = benchmark::Counter(1);
      st.counters["tmop_failed"] = benchmark::Counter(counters.tmopFailures);
      st.SetItemsProcessed(st.iterations() * counters.finalCells);
    }
  }

#define RODIN_NO_CUT_P2_TMOP_RELABEL_ORBIT(SHAPE_LABEL, SHAPE_ENUM, RESOLUTION, STEPS) \
  BENCHMARK(BM_LevelSetPipeline_CurvedTargetCarryForwardOrbit) \
    ->Name("LevelSetPipeline/NoCutP2TMOPRelabel/" SHAPE_LABEL "/ShapeSize/Projected/" #RESOLUTION "x" #STEPS) \
    ->Args({RESOLUTION, STEPS, static_cast<int>(ShapeCase::SHAPE_ENUM), \
      static_cast<int>(MetricCase::ShapeSize), \
      static_cast<int>(CurvedTargetCase::ProjectedInterface)}) \
    ->Unit(benchmark::kMillisecond);

#define RODIN_MMG_CUT_P2_TMOP_RELABEL_ORBIT(SHAPE_LABEL, SHAPE_ENUM, RESOLUTION, STEPS) \
  BENCHMARK(BM_LevelSetPipeline_CurvedMMGCarryForwardOrbit) \
    ->Name("LevelSetPipeline/MMGCutP2TMOPRelabel/" SHAPE_LABEL "/ShapeSize/Projected/" #RESOLUTION "x" #STEPS) \
    ->Args({RESOLUTION, STEPS, static_cast<int>(ShapeCase::SHAPE_ENUM), \
      static_cast<int>(MetricCase::ShapeSize), \
      static_cast<int>(CurvedTargetCase::ProjectedInterface)}) \
    ->Unit(benchmark::kMillisecond);

  RODIN_NO_CUT_P2_TMOP_RELABEL_ORBIT("Circle", CircleOrbit, 5, 8)
  RODIN_NO_CUT_P2_TMOP_RELABEL_ORBIT("Circle", CircleOrbit, 10, 8)
  RODIN_NO_CUT_P2_TMOP_RELABEL_ORBIT("Circle", CircleOrbit, 16, 8)
  RODIN_NO_CUT_P2_TMOP_RELABEL_ORBIT("Circle", CircleOrbit, 24, 8)
  RODIN_NO_CUT_P2_TMOP_RELABEL_ORBIT("WavyCircleFigureEight", WavyCircleFigureEight, 5, 8)
  RODIN_NO_CUT_P2_TMOP_RELABEL_ORBIT("WavyCircleFigureEight", WavyCircleFigureEight, 10, 8)
  RODIN_NO_CUT_P2_TMOP_RELABEL_ORBIT("WavyCircleFigureEight", WavyCircleFigureEight, 16, 8)
  RODIN_NO_CUT_P2_TMOP_RELABEL_ORBIT("WavyCircleFigureEight", WavyCircleFigureEight, 24, 8)

  RODIN_MMG_CUT_P2_TMOP_RELABEL_ORBIT("Circle", CircleOrbit, 5, 8)
  RODIN_MMG_CUT_P2_TMOP_RELABEL_ORBIT("Circle", CircleOrbit, 10, 8)
  RODIN_MMG_CUT_P2_TMOP_RELABEL_ORBIT("Circle", CircleOrbit, 16, 8)
  RODIN_MMG_CUT_P2_TMOP_RELABEL_ORBIT("Circle", CircleOrbit, 24, 8)
  RODIN_MMG_CUT_P2_TMOP_RELABEL_ORBIT("WavyCircleFigureEight", WavyCircleFigureEight, 5, 8)
  RODIN_MMG_CUT_P2_TMOP_RELABEL_ORBIT("WavyCircleFigureEight", WavyCircleFigureEight, 10, 8)
  RODIN_MMG_CUT_P2_TMOP_RELABEL_ORBIT("WavyCircleFigureEight", WavyCircleFigureEight, 16, 8)
  RODIN_MMG_CUT_P2_TMOP_RELABEL_ORBIT("WavyCircleFigureEight", WavyCircleFigureEight, 24, 8)

#undef RODIN_NO_CUT_P2_TMOP_RELABEL_ORBIT
#undef RODIN_MMG_CUT_P2_TMOP_RELABEL_ORBIT
}
