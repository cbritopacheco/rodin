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
#include <vector>

#include <benchmark/benchmark.h>

#include <Rodin/Adaptation.h>
#include <Rodin/Adaptation/TargetMatrixOptimization/IsoparametricGeometry.h>
#include <Rodin/Adaptation/TargetMatrixOptimization/Metrics.h>
#include <Rodin/Assembly/Default.h>
#include <Rodin/Geometry.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Solver/SparseLU.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Adaptation::TargetMatrixOptimization;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace Rodin::Tests::Benchmarks
{
  namespace
  {
    constexpr Real Pi = 3.14159265358979323846;
    constexpr Attribute Interface = 30;
    constexpr Attribute Negative = 1;
    constexpr Attribute Positive = 10;
    constexpr size_t CurvedQuadratureOrder = 4;
    constexpr size_t ClassificationQuadratureOrder = 8;

    using CellSpace = P0<Real, LocalMesh>;
    using ScalarData = Math::Vector<Real>;
    using CellField = GridFunction<CellSpace, ScalarData>;

    enum class TermCase : int
    {
      QualityDeviation = 0,
      FitDeviation = 1,
      PhaseDeviation = 2,
      FitPhaseDeviation = 3,
      AllTerms = 4
    };

    enum class WeightCase : int
    {
      Balanced = 0,
      Fit2x = 1,
      Phase2x = 2
    };

    struct Weights
    {
      Real quality = 1;
      Real fit = 1;
      Real phase = Real(0.1);
      Real deviation = 1;
      Real targetBlend = Real(0.1);
      Real phaseEpsilonFactor = Real(0.5);
    };

    Weights weightsFor(WeightCase c)
    {
      Weights w;
      switch (c)
      {
        case WeightCase::Balanced:
          return w;
        case WeightCase::Fit2x:
          w.fit = 2;
          return w;
        case WeightCase::Phase2x:
          w.phase = Real(0.2);
          return w;
      }
      return w;
    }

    template <int Term>
    constexpr bool usesQuality()
    {
      return Term == static_cast<int>(TermCase::QualityDeviation)
          || Term == static_cast<int>(TermCase::AllTerms);
    }

    template <int Term>
    constexpr bool usesFit()
    {
      return Term == static_cast<int>(TermCase::FitDeviation)
          || Term == static_cast<int>(TermCase::FitPhaseDeviation)
          || Term == static_cast<int>(TermCase::AllTerms);
    }

    template <int Term>
    constexpr bool usesPhase()
    {
      return Term == static_cast<int>(TermCase::PhaseDeviation)
          || Term == static_cast<int>(TermCase::FitPhaseDeviation)
          || Term == static_cast<int>(TermCase::AllTerms);
    }

    Math::SpatialPoint center(Real t)
    {
      return Math::SpatialPoint{
        Real(0.5) + Real(0.18) * std::sin(Real(2) * Pi * t),
        Real(0.5) + Real(0.13) * std::sin(Real(4) * Pi * t + Real(0.35))
      };
    }

    Real phiAt(const Math::SpatialPoint& x, Real t)
    {
      const auto c = center(t);
      const Real dx = x[0] - c[0];
      const Real dy = x[1] - c[1];
      return std::hypot(dx, dy) - Real(0.17);
    }

    Math::SpatialPoint gradPhiAt(const Math::SpatialPoint& x, Real t)
    {
      const auto c = center(t);
      const Real dx = x[0] - c[0];
      const Real dy = x[1] - c[1];
      const Real rho = std::hypot(dx, dy);
      if (rho <= Real(1e-14))
        return Math::SpatialPoint{ Real(1), Real(0) };
      return Math::SpatialPoint{ dx / rho, dy / rho };
    }

    Math::SpatialMatrix<Real> hessianPhiAt(const Math::SpatialPoint& x, Real t)
    {
      Math::SpatialMatrix<Real> h(2, 2);
      h.setZero();
      const auto c = center(t);
      const Real dx = x[0] - c[0];
      const Real dy = x[1] - c[1];
      const Real rho2 = dx * dx + dy * dy;
      const Real rho = std::sqrt(rho2);
      if (rho <= Real(1e-14))
        return h;

      const Real rho3 = rho2 * rho;
      h(0, 0) = dy * dy / rho3;
      h(0, 1) = -dx * dy / rho3;
      h(1, 0) = h(0, 1);
      h(1, 1) = dx * dx / rho3;
      return h;
    }

    class LevelSetSignFunction final
      : public RealFunctionBase<LevelSetSignFunction>
    {
      public:
        LevelSetSignFunction(Real t, size_t order)
          : m_t(t),
            m_order(order)
        {}

        Real getValue(const Point& p) const
        {
          return phiAt(p.getPhysicalCoordinates(), m_t) <= Real(0)
            ? Real(-1)
            : Real(1);
        }

        Optional<size_t> getOrder(const Polytope&) const noexcept
        {
          return m_order;
        }

        LevelSetSignFunction* copy() const noexcept override
        {
          return new LevelSetSignFunction(*this);
        }

      private:
        Real m_t;
        size_t m_order;
    };

    struct CellClassificationStats
    {
      Index negative = 0;
      Index positive = 0;
      Index changed = 0;
    };

    std::vector<int> captureCellSigns(const LocalMesh& mesh)
    {
      std::vector<int> out(mesh.getCellCount(), 0);
      for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
      {
        const auto attr = mesh.getAttribute(2, c);
        if (!attr)
          continue;
        if (*attr == Negative)
          out[static_cast<size_t>(c)] = -1;
        else if (*attr == Positive)
          out[static_cast<size_t>(c)] = 1;
      }
      return out;
    }

    Index countChangedSigns(
        const std::vector<int>& before,
        const std::vector<int>& after)
    {
      const size_t n = std::min(before.size(), after.size());
      Index changed = 0;
      for (size_t i = 0; i < n; ++i)
        if (before[i] != after[i])
          ++changed;
      changed += static_cast<Index>(before.size() > after.size()
        ? before.size() - after.size()
        : after.size() - before.size());
      return changed;
    }

    CellClassificationStats classifyCellsByProjection(
        LocalMesh& mesh,
        Real t,
        CellField& classification)
    {
      const auto before = captureCellSigns(mesh);

      CellSpace p0(mesh);
      TrialFunction cellValue(p0);
      TestFunction test(p0);
      LevelSetSignFunction sign(t, ClassificationQuadratureOrder);

      Problem projection(cellValue, test);
      projection = Integral(cellValue, test) - Integral(sign, test);
      Solver::SparseLU(projection).solve();
      classification.getData() = cellValue.getSolution().getData();

      CellClassificationStats stats;
      for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
      {
        const Attribute attr =
          classification.getData()[c] <= Real(0) ? Negative : Positive;
        mesh.setAttribute({ 2, c }, attr);
        if (attr == Negative)
          ++stats.negative;
        else
          ++stats.positive;
      }
      stats.changed = countChangedSigns(before, captureCellSigns(mesh));
      return stats;
    }

    void clearInterfaceEdges(LocalMesh& mesh)
    {
      mesh.getConnectivity().compute(1, 0);
      const auto& conn = mesh.getConnectivity();
      for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
      {
        const auto attr = mesh.getAttribute(1, e);
        if (attr && *attr == Interface)
          mesh.setAttribute({ 1, e }, Optional<Attribute>{});
      }
    }

    std::vector<char> captureInterfaceSkeleton(const LocalMesh& mesh)
    {
      const auto& conn = mesh.getConnectivity();
      std::vector<char> out(conn.getCount(1), 0);
      for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
      {
        const auto attr = mesh.getAttribute(1, e);
        if (attr && *attr == Interface)
          out[static_cast<size_t>(e)] = 1;
      }
      return out;
    }

    Index countChangedSkeleton(
        const std::vector<char>& before,
        const std::vector<char>& after)
    {
      const size_t n = std::min(before.size(), after.size());
      Index changed = 0;
      for (size_t i = 0; i < n; ++i)
        if (before[i] != after[i])
          ++changed;
      changed += static_cast<Index>(before.size() > after.size()
        ? before.size() - after.size()
        : after.size() - before.size());
      return changed;
    }

    struct InterfaceStats
    {
      Index edgeCount = 0;
      Index maxDegree = 0;
      Index branchVertices = 0;
    };

    InterfaceStats computeInterfaceStats(const LocalMesh& mesh)
    {
      InterfaceStats stats;
      const auto& conn = mesh.getConnectivity();
      std::vector<Index> degree(mesh.getVertexCount(), 0);
      for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
      {
        const auto attr = mesh.getAttribute(1, e);
        if (!attr || *attr != Interface)
          continue;
        ++stats.edgeCount;
        const auto& edge = conn.getPolytope(1, e);
        ++degree[edge(0)];
        ++degree[edge(1)];
      }
      for (Index d : degree)
      {
        stats.maxDegree = std::max(stats.maxDegree, d);
        if (d > 2)
          ++stats.branchVertices;
      }
      return stats;
    }

    InterfaceStats markAttributeJumpInterfaceEdges(LocalMesh& mesh)
    {
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 2);
      mesh.getConnectivity().compute(1, 0);
      const auto& conn = mesh.getConnectivity();
      for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
      {
        const auto adj = conn.getIncidence({ 1, 2 }, e);
        if (adj.size() != 2)
          continue;
        const auto a0 = mesh.getAttribute(2, adj[0]);
        const auto a1 = mesh.getAttribute(2, adj[1]);
        if (!a0 || !a1)
          continue;
        const bool jump = (*a0 == Negative && *a1 == Positive)
                       || (*a0 == Positive && *a1 == Negative);
        if (jump)
          mesh.setAttribute({ 1, e }, Interface);
      }
      return computeInterfaceStats(mesh);
    }

    template <class Element>
    void syncInterfaceEdgesFromCells(LocalMesh& mesh, Attribute interfaceAttribute)
    {
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 2);
      mesh.getConnectivity().compute(1, 0);

      auto& conn = mesh.getConnectivity();
      Element feSeg(Polytope::Type::Segment);
      const auto& refSeg = Element::getNodes(Polytope::Type::Segment);
      static const std::array<Math::SpatialPoint, 3> triRef = {{
        Math::SpatialPoint{ Real(0), Real(0) },
        Math::SpatialPoint{ Real(1), Real(0) },
        Math::SpatialPoint{ Real(0), Real(1) } }};

      for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
      {
        const auto attr = mesh.getAttribute(1, e);
        if (!attr || *attr != interfaceAttribute)
          continue;

        const auto adj = conn.getIncidence({ 1, 2 }, e);
        if (adj.empty())
          continue;

        const auto& edge = conn.getPolytope(1, e);
        const auto cellIndex = adj[0];
        const auto& cell = conn.getPolytope(2, cellIndex);

        size_t localA = 3;
        size_t localB = 3;
        for (size_t i = 0; i < 3; ++i)
        {
          if (cell(i) == edge(0))
            localA = i;
          if (cell(i) == edge(1))
            localB = i;
        }
        if (localA >= 3 || localB >= 3)
          continue;

        auto cellIt = mesh.getPolytope(2, cellIndex);
        PointCloud pc(2, refSeg.size());
        for (size_t j = 0; j < refSeg.size(); ++j)
        {
          const Real s = refSeg[j][0];
          const Math::SpatialPoint rc{
            (Real(1) - s) * triRef[localA][0] + s * triRef[localB][0],
            (Real(1) - s) * triRef[localA][1] + s * triRef[localB][1]
          };
          const auto x = Point(*cellIt, rc).getPhysicalCoordinates();
          pc(0, j) = x[0];
          pc(1, j) = x[1];
        }

        mesh.setPolytopeTransformation(
            { size_t(1), e },
            new ParametricTransformation<Element>(pc, feSeg));
      }
    }

    struct LabelReport
    {
      CellClassificationStats cells;
      InterfaceStats interface;
      Index changedInterfaceEdges = 0;
    };

    template <class Element>
    LabelReport labelAndBuildInterface(
        LocalMesh& mesh,
        Real t,
        CellField& classification)
    {
      const auto before = captureInterfaceSkeleton(mesh);
      clearInterfaceEdges(mesh);

      LabelReport report;
      report.cells = classifyCellsByProjection(mesh, t, classification);
      report.interface = markAttributeJumpInterfaceEdges(mesh);
      syncInterfaceEdgesFromCells<Element>(mesh, Interface);
      report.changedInterfaceEdges =
        countChangedSkeleton(before, captureInterfaceSkeleton(mesh));
      return report;
    }

    Real totalCellMeasure(LocalMesh& mesh)
    {
      Real measure = 0;
      const auto& conn = mesh.getConnectivity();
      const auto& qf = QF::PolytopeQuadratureFormula::get(
          CurvedQuadratureOrder, Polytope::Type::Triangle);
      for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
      {
        if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
          continue;
        const auto cellIt = mesh.getPolytope(2, c);
        for (size_t q = 0; q < qf.getSize(); ++q)
        {
          Math::SpatialMatrix<Real> J;
          cellIt->getTransformation().jacobian(J, qf.getPoint(q));
          if (J.rows() == 2 && J.cols() == 2)
            measure += qf.getWeight(q) * std::abs(J.determinant());
        }
      }
      return measure;
    }

    Real edgeMeasure(LocalMesh& mesh, Attribute attribute)
    {
      Real measure = 0;
      const auto& conn = mesh.getConnectivity();
      const auto& qf = QF::PolytopeQuadratureFormula::get(
          CurvedQuadratureOrder, Polytope::Type::Segment);
      for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
      {
        const auto attr = mesh.getAttribute(1, e);
        if (!attr || *attr != attribute)
          continue;
        const auto edgeIt = mesh.getPolytope(1, e);
        for (size_t q = 0; q < qf.getSize(); ++q)
        {
          Math::SpatialMatrix<Real> J;
          edgeIt->getTransformation().jacobian(J, qf.getPoint(q));
          if (J.rows() == 2 && J.cols() == 1)
            measure += qf.getWeight(q) * std::hypot(J(0, 0), J(1, 0));
        }
      }
      return measure;
    }

    template <size_t Order>
    LocalMesh makeMesh(size_t resolution)
    {
      const Real h = Real(1) / static_cast<Real>(resolution - 1);
      LocalMesh mesh =
        LocalMesh::UniformGrid(Polytope::Type::Triangle, { resolution, resolution });
      mesh.scale(h);
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 2);
      mesh.getConnectivity().compute(1, 0);

      RealH1Element<Order> fe(Polytope::Type::Triangle);
      upgradeTransformations(mesh, fe);
      return mesh;
    }

    struct SolveStats
    {
      bool converged = false;
      bool moved = false;
      size_t iterations = 0;
      Real initialResidual = 0;
      Real finalResidual = 0;
      Real finalStepNorm = 0;
      Real fitEnergyInitial = 0;
      Real fitEnergyFinal = 0;
      Real phaseEnergyInitial = 0;
      Real phaseEnergyFinal = 0;
      Index wrongSideInitial = 0;
      Index wrongSideFinal = 0;
    };

    template <size_t Order, int Term>
    SolveStats solveTMOPStep(
        LocalMesh& mesh,
        Real t,
        Real h,
        const Weights& weights)
    {
      SolveStats stats;

      auto phiValue = [t](const Math::SpatialPoint& x) { return phiAt(x, t); };
      auto phiGradient = [t](const Math::SpatialPoint& x) { return gradPhiAt(x, t); };
      auto phiHessian = [t](const Math::SpatialPoint& x) { return hessianPhiAt(x, t); };

      ShapeSizeBlendMetric metric(Real(0.5));
      CurvedQualityTargetJacobian target(mesh, weights.targetBlend);

      VectorH1<Order, LocalMesh> space(
          std::integral_constant<size_t, Order>{}, mesh, 2);
      GridFunction u(space);
      u.getData().setZero();
      TrialFunction du(space);
      TestFunction v(space);

      QualityTerm quality(metric, target, weights.quality);
      quality.setQuadratureOrder(CurvedQuadratureOrder);

      DeviationTerm deviation(weights.deviation);

      AnalyticLevelSetFitTerm fit(
          phiValue,
          phiGradient,
          phiHessian,
          Optional<Attribute>(Interface),
          weights.fit);
      fit
        .setQuadratureOrder(CurvedQuadratureOrder)
        .setNormalization(std::max(edgeMeasure(mesh, Interface), Real(1e-12)));

      VolumetricPhaseConsistencyTerm phase(
          phiValue,
          phiGradient,
          phiHessian,
          Negative,
          Positive,
          weights.phase);
      phase
        .setQuadratureOrder(CurvedQuadratureOrder)
        .setEpsilon(std::max(weights.phaseEpsilonFactor * h, Real(1e-12)))
        .setMargin(Real(1))
        .setNormalization(std::max(totalCellMeasure(mesh), Real(1e-12)));

      if constexpr (usesFit<Term>())
        stats.fitEnergyInitial = fit.energy(u);
      if constexpr (usesPhase<Term>())
      {
        stats.phaseEnergyInitial = phase.energy(u);
        stats.wrongSideInitial = phase.countWrongSideQuadrature(mesh);
      }

      Problem problem(du, v);
      if constexpr (Term == static_cast<int>(TermCase::QualityDeviation))
      {
        problem = quality.tangent(u, du, v)
                + deviation.tangent(du, v)
                + quality.residual(u, v)
                + deviation.residual(u, v);
      }
      else if constexpr (Term == static_cast<int>(TermCase::FitDeviation))
      {
        problem = fit.tangent(u, du, v)
                + deviation.tangent(du, v)
                + fit.residual(u, v)
                + deviation.residual(u, v);
      }
      else if constexpr (Term == static_cast<int>(TermCase::PhaseDeviation))
      {
        problem = phase.tangent(u, du, v)
                + deviation.tangent(du, v)
                + phase.residual(u, v)
                + deviation.residual(u, v);
      }
      else if constexpr (Term == static_cast<int>(TermCase::FitPhaseDeviation))
      {
        problem = fit.tangent(u, du, v)
                + phase.tangent(u, du, v)
                + deviation.tangent(du, v)
                + fit.residual(u, v)
                + phase.residual(u, v)
                + deviation.residual(u, v);
      }
      else
      {
        problem = quality.tangent(u, du, v)
                + fit.tangent(u, du, v)
                + phase.tangent(u, du, v)
                + deviation.tangent(du, v)
                + quality.residual(u, v)
                + fit.residual(u, v)
                + phase.residual(u, v)
                + deviation.residual(u, v);
      }

      Solver::SparseLU linearSolver(problem);
      Solver::NewtonSolver newton(linearSolver);
      newton
        .setMaxIterations(20)
        .setDampingFactor(Real(1))
        .setAbsoluteTolerance(Real(1e-10))
        .setRelativeTolerance(Real(1e-8))
        .setStepTolerance(Real(1e-10));
      newton.solve(u);

      const auto& nr = newton.getReport();
      stats.converged = nr.converged;
      stats.iterations = nr.iterations;
      stats.initialResidual = nr.initial_residual;
      stats.finalResidual = nr.final_residual;
      stats.finalStepNorm = nr.final_step_norm;

      if constexpr (usesFit<Term>())
        stats.fitEnergyFinal = fit.energy(u);
      if constexpr (usesPhase<Term>())
        stats.phaseEnergyFinal = phase.energy(u);

      if (nr.converged && u.getData().allFinite())
      {
        RealH1Element<Order> fe(Polytope::Type::Triangle);
        moveMesh(mesh, u, fe, Interface);
        stats.moved = true;
      }

      if constexpr (usesPhase<Term>())
        stats.wrongSideFinal = phase.countWrongSideQuadrature(mesh);
      return stats;
    }

    struct PipelineStats
    {
      Index cells = 0;
      Index negative = 0;
      Index positive = 0;
      Index changedCells = 0;
      Index interfaceEdges = 0;
      Index changedInterfaceEdges = 0;
      Index interfaceMaxDegree = 0;
      Index branchVertices = 0;
      Index convergedSteps = 0;
      Index movedSteps = 0;
      Index newtonIterations = 0;
      Real finalResidual = 0;
      Real finalStepNorm = 0;
      Real fitEnergyInitial = 0;
      Real fitEnergyFinal = 0;
      Real phaseEnergyInitial = 0;
      Real phaseEnergyFinal = 0;
      Index wrongSideInitial = 0;
      Index wrongSideFinal = 0;
      Real curvedFitRms = 0;
      Real curvedFitMax = 0;
      Real curvedQmin = 0;
      Real curvedMinDet = 0;
      Index curvedInvalid = 0;
      Index curvedOverlap = 0;
      Real targetMaxMetric = 0;
    };

    template <size_t Order, int Term>
    PipelineStats runPipeline(
        size_t resolution,
        Index steps,
        WeightCase weightCase)
    {
      const Real h = Real(1) / static_cast<Real>(resolution - 1);
      const auto weights = weightsFor(weightCase);
      LocalMesh mesh = makeMesh<Order>(resolution);
      RealH1Element<Order> fe(Polytope::Type::Triangle);

      PipelineStats totals;
      for (Index step = 0; step < steps; ++step)
      {
        const Real t = static_cast<Real>(step)
          / static_cast<Real>(std::max<Index>(1, steps - 1));

        CellSpace cellSpace(mesh);
        CellField classification(cellSpace);
        labelAndBuildInterface<RealH1Element<Order>>(mesh, t, classification);

        const auto solve = solveTMOPStep<Order, Term>(mesh, t, h, weights);
        const auto label =
          labelAndBuildInterface<RealH1Element<Order>>(mesh, t, classification);

        auto phiValue = [t](const Math::SpatialPoint& x) { return phiAt(x, t); };
        const auto curved = curvedMetrics(mesh, phiValue, Interface);
        const auto targetStats = targetQualityMetrics(
            mesh,
            ShapeSizeBlendMetric(Real(0.5)),
            CurvedQualityTargetJacobian(mesh, weights.targetBlend),
            Polytope::Type::Triangle,
            static_cast<Index>(CurvedQuadratureOrder));

        totals.cells = static_cast<Index>(mesh.getCellCount());
        totals.negative += label.cells.negative;
        totals.positive += label.cells.positive;
        totals.changedCells += label.cells.changed;
        totals.interfaceEdges += label.interface.edgeCount;
        totals.changedInterfaceEdges += label.changedInterfaceEdges;
        totals.interfaceMaxDegree =
          std::max(totals.interfaceMaxDegree, label.interface.maxDegree);
        totals.branchVertices += label.interface.branchVertices;
        totals.convergedSteps += solve.converged ? 1 : 0;
        totals.movedSteps += solve.moved ? 1 : 0;
        totals.newtonIterations += static_cast<Index>(solve.iterations);
        totals.finalResidual += solve.finalResidual;
        totals.finalStepNorm += solve.finalStepNorm;
        totals.fitEnergyInitial += solve.fitEnergyInitial;
        totals.fitEnergyFinal += solve.fitEnergyFinal;
        totals.phaseEnergyInitial += solve.phaseEnergyInitial;
        totals.phaseEnergyFinal += solve.phaseEnergyFinal;
        totals.wrongSideInitial += solve.wrongSideInitial;
        totals.wrongSideFinal += solve.wrongSideFinal;
        totals.curvedFitRms += curved.fitRms;
        totals.curvedFitMax = std::max(totals.curvedFitMax, curved.fitMax);
        totals.curvedQmin += curved.qmin;
        totals.curvedMinDet += curved.minDet;
        totals.curvedInvalid += curved.invalidJacobianSamples;
        totals.curvedOverlap += curved.overlapSamples;
        totals.targetMaxMetric += targetStats.maxMetric;
      }
      benchmark::DoNotOptimize(mesh.getCellCount());
      return totals;
    }

    void publish(
        benchmark::State& st,
        const PipelineStats& stats,
        size_t resolution,
        Index steps,
        WeightCase weightCase,
        size_t order,
        int term)
    {
      const Real inv = Real(1) / static_cast<Real>(std::max<Index>(1, steps));
      st.counters["resolution"] = benchmark::Counter(resolution);
      st.counters["steps"] = benchmark::Counter(steps);
      st.counters["geometry_order"] = benchmark::Counter(order);
      st.counters["displacement_order"] = benchmark::Counter(order);
      st.counters["term_case"] = benchmark::Counter(term);
      st.counters["weight_case"] =
        benchmark::Counter(static_cast<int>(weightCase));
      st.counters["avg_negative_cells"] =
        benchmark::Counter(static_cast<Real>(stats.negative) * inv);
      st.counters["avg_positive_cells"] =
        benchmark::Counter(static_cast<Real>(stats.positive) * inv);
      st.counters["avg_changed_cells"] =
        benchmark::Counter(static_cast<Real>(stats.changedCells) * inv);
      st.counters["avg_interface_edges"] =
        benchmark::Counter(static_cast<Real>(stats.interfaceEdges) * inv);
      st.counters["avg_changed_interface_edges"] =
        benchmark::Counter(static_cast<Real>(stats.changedInterfaceEdges) * inv);
      st.counters["interface_max_degree"] =
        benchmark::Counter(stats.interfaceMaxDegree);
      st.counters["avg_branch_vertices"] =
        benchmark::Counter(static_cast<Real>(stats.branchVertices) * inv);
      st.counters["converged_steps"] =
        benchmark::Counter(stats.convergedSteps);
      st.counters["moved_steps"] =
        benchmark::Counter(stats.movedSteps);
      st.counters["avg_newton_iterations"] =
        benchmark::Counter(static_cast<Real>(stats.newtonIterations) * inv);
      st.counters["avg_final_residual"] =
        benchmark::Counter(stats.finalResidual * inv);
      st.counters["avg_final_step_norm"] =
        benchmark::Counter(stats.finalStepNorm * inv);
      st.counters["avg_fit_energy_initial"] =
        benchmark::Counter(stats.fitEnergyInitial * inv);
      st.counters["avg_fit_energy_final"] =
        benchmark::Counter(stats.fitEnergyFinal * inv);
      st.counters["avg_phase_energy_initial"] =
        benchmark::Counter(stats.phaseEnergyInitial * inv);
      st.counters["avg_phase_energy_final"] =
        benchmark::Counter(stats.phaseEnergyFinal * inv);
      st.counters["avg_wrong_side_initial"] =
        benchmark::Counter(static_cast<Real>(stats.wrongSideInitial) * inv);
      st.counters["avg_wrong_side_final"] =
        benchmark::Counter(static_cast<Real>(stats.wrongSideFinal) * inv);
      st.counters["avg_curved_fit_rms"] =
        benchmark::Counter(stats.curvedFitRms * inv);
      st.counters["max_curved_fit"] =
        benchmark::Counter(stats.curvedFitMax);
      st.counters["avg_curved_qmin"] =
        benchmark::Counter(stats.curvedQmin * inv);
      st.counters["avg_curved_min_det"] =
        benchmark::Counter(stats.curvedMinDet * inv);
      st.counters["curved_invalid"] =
        benchmark::Counter(stats.curvedInvalid);
      st.counters["curved_overlap"] =
        benchmark::Counter(stats.curvedOverlap);
      st.counters["avg_target_max_metric"] =
        benchmark::Counter(stats.targetMaxMetric * inv);
      st.SetItemsProcessed(st.iterations() * stats.cells * steps);
    }

    template <size_t Order, int Term>
    static void BM_LevelSetPipeline_TMOP(benchmark::State& st)
    {
      const auto resolution = static_cast<size_t>(st.range(0));
      const auto steps = static_cast<Index>(st.range(1));
      const auto weightCase = static_cast<WeightCase>(st.range(2));

      PipelineStats stats;
      for (auto _ : st)
      {
        stats = runPipeline<Order, Term>(resolution, steps, weightCase);
        benchmark::DoNotOptimize(stats.cells);
      }
      publish(st, stats, resolution, steps, weightCase, Order, Term);
    }
  }

#define RODIN_TMOP_BENCHMARK(ORDER, ORDER_LABEL, TERM, TERM_LABEL, WEIGHT, WEIGHT_LABEL, RESOLUTION, STEPS) \
  BENCHMARK_TEMPLATE(BM_LevelSetPipeline_TMOP, ORDER, static_cast<int>(TermCase::TERM)) \
    ->Name("LevelSetPipeline/TMOP/" ORDER_LABEL "/" TERM_LABEL "/" WEIGHT_LABEL "/" #RESOLUTION "x" #STEPS) \
    ->Args({ RESOLUTION, STEPS, static_cast<int>(WeightCase::WEIGHT) }) \
    ->Unit(benchmark::kMillisecond);

#define RODIN_TMOP_WEIGHT_MATRIX(ORDER, ORDER_LABEL, TERM, TERM_LABEL, RESOLUTION, STEPS) \
  RODIN_TMOP_BENCHMARK(ORDER, ORDER_LABEL, TERM, TERM_LABEL, Balanced, "Balanced", RESOLUTION, STEPS) \
  RODIN_TMOP_BENCHMARK(ORDER, ORDER_LABEL, TERM, TERM_LABEL, Fit2x, "Fit2x", RESOLUTION, STEPS) \
  RODIN_TMOP_BENCHMARK(ORDER, ORDER_LABEL, TERM, TERM_LABEL, Phase2x, "Phase2x", RESOLUTION, STEPS)

#define RODIN_TMOP_TERM_MATRIX(ORDER, ORDER_LABEL, RESOLUTION, STEPS) \
  RODIN_TMOP_WEIGHT_MATRIX(ORDER, ORDER_LABEL, QualityDeviation, "QualityDeviation", RESOLUTION, STEPS) \
  RODIN_TMOP_WEIGHT_MATRIX(ORDER, ORDER_LABEL, FitDeviation, "FitDeviation", RESOLUTION, STEPS) \
  RODIN_TMOP_WEIGHT_MATRIX(ORDER, ORDER_LABEL, PhaseDeviation, "PhaseDeviation", RESOLUTION, STEPS) \
  RODIN_TMOP_WEIGHT_MATRIX(ORDER, ORDER_LABEL, FitPhaseDeviation, "FitPhaseDeviation", RESOLUTION, STEPS) \
  RODIN_TMOP_WEIGHT_MATRIX(ORDER, ORDER_LABEL, AllTerms, "AllTerms", RESOLUTION, STEPS)

#define RODIN_TMOP_ORDER_MATRIX(ORDER, ORDER_LABEL) \
  RODIN_TMOP_TERM_MATRIX(ORDER, ORDER_LABEL, 8, 4) \
  RODIN_TMOP_TERM_MATRIX(ORDER, ORDER_LABEL, 16, 4)

  RODIN_TMOP_ORDER_MATRIX(1, "P1_P1")
  RODIN_TMOP_ORDER_MATRIX(2, "P2_P2")

#undef RODIN_TMOP_ORDER_MATRIX
#undef RODIN_TMOP_TERM_MATRIX
#undef RODIN_TMOP_WEIGHT_MATRIX
#undef RODIN_TMOP_BENCHMARK
}
