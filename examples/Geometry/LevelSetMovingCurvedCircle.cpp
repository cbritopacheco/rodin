/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @example Minimal P2 TMOP fitting of a moving curved level-set circle.
 *
 * Pipeline:
 *   P2 curved mesh from the start
 *   -> classify cells with a high-order L2 projection into P0
 *   -> mark material interface edges
 *   -> one undamped Newton solve for the P2 displacement
 *   -> move the mesh, relabel, and write diagnostics.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include <Rodin/Adaptation.h>
#include <Rodin/Adaptation/TargetMatrixOptimization/IsoparametricGeometry.h>
#include <Rodin/Adaptation/TargetMatrixOptimization/Metrics.h>
#include <Rodin/Assembly/Default.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include <Rodin/Solver/SparseLU.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Adaptation::TargetMatrixOptimization;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

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
  using P1Space = P1<Real, LocalMesh>;
  using P1Field = GridFunction<P1Space, ScalarData>;
  using P2Element = RealH1Element<2>;

  Math::SpatialPoint center(Real t)
  {
    return Math::SpatialPoint{
      Real(0.5) + Real(0.18) * std::sin(Real(2) * Pi * t),
      Real(0.5) + Real(0.13) * std::sin(Real(4) * Pi * t + Real(0.35))
    };
  }

  Real radius()
  {
    return Real(0.17);
  }

  Real phiAt(const Math::SpatialPoint& x, Real t)
  {
    const auto c = center(t);
    const Real dx = x[0] - c[0];
    const Real dy = x[1] - c[1];
    return std::hypot(dx, dy) - radius();
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

  struct CellClassificationStats
  {
    Index negative = 0;
    Index positive = 0;
    Index changed = 0;
  };

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

  void fillClassificationView(
      const LocalMesh& mesh,
      const CellField& classification,
      P1Field& view)
  {
    const auto& conn = mesh.getConnectivity();
    for (Index v = 0; v < static_cast<Index>(mesh.getVertexCount()); ++v)
    {
      Real sum = 0;
      Index count = 0;
      for (Index c : conn.getIncidence({ 0, 2 }, v))
      {
        if (c >= 0 && c < classification.getData().size())
        {
          sum += classification.getData()[c];
          ++count;
        }
      }
      view[v] = count > 0 ? sum / static_cast<Real>(count) : Real(0);
    }
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
  void syncInterfaceEdgesFromCells(
      LocalMesh& mesh,
      const Element&,
      Attribute interfaceAttribute)
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
      const Element& fe,
      Real t,
      CellField& classification)
  {
    const auto before = captureInterfaceSkeleton(mesh);
    clearInterfaceEdges(mesh);

    LabelReport report;
    report.cells = classifyCellsByProjection(mesh, t, classification);
    report.interface = markAttributeJumpInterfaceEdges(mesh);
    syncInterfaceEdgesFromCells(mesh, fe, Interface);
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

  struct SolveReport
  {
    bool moved = false;
    bool minimizerConverged = false;
    size_t minimizerIterations = 0;
    Index acceptedSteps = 0;
    Index rejectedSteps = 0;
    Real initialEnergy = 0;
    Real finalEnergy = 0;
    Real initialGradientNorm = 0;
    Real finalGradientNorm = 0;
    Real finalStepNorm = 0;
    Real finalAlpha = 0;
    Real finalDirectionalDerivative = 0;
    Real fitEnergyInitial = 0;
    Real fitEnergyFinal = 0;
    Real phaseEnergyInitial = 0;
    Real phaseEnergyFinal = 0;
    Index wrongSideInitial = 0;
    Index wrongSideFinal = 0;
    Math::Vector<Real> displacement;
  };

  SolveReport solveTMOPStep(
      LocalMesh& mesh,
      Real t,
      Real h,
      Real qualityWeight,
      Real fitWeight,
      Real phaseWeight,
      Real deviationWeight,
      Real targetBlend,
      Real phaseEpsilonFactor,
      Index maxIterations,
      const P2Element& fe)
  {
    SolveReport report;

    auto phiValue = [t](const Math::SpatialPoint& x) { return phiAt(x, t); };
    auto phiGradient = [t](const Math::SpatialPoint& x) { return gradPhiAt(x, t); };
    auto phiHessian = [t](const Math::SpatialPoint& x) { return hessianPhiAt(x, t); };

    ShapeSizeBlendMetric metric(Real(0.5));
    CurvedQualityTargetJacobian target(mesh, targetBlend);

    VectorH1<2, LocalMesh> space(std::integral_constant<size_t, 2>{}, mesh, 2);
    GridFunction u(space);
    u.getData().setZero();
    TrialFunction p(space);
    TestFunction v(space);

    QualityTerm quality(metric, target, qualityWeight);
    quality.setQuadratureOrder(CurvedQuadratureOrder);

    DeviationTerm deviation(deviationWeight);
    AnalyticLevelSetFitTerm fit(
        phiValue,
        phiGradient,
        phiHessian,
        Optional<Attribute>(Interface),
        fitWeight);
    fit
      .setQuadratureOrder(CurvedQuadratureOrder)
      .setNormalization(std::max(edgeMeasure(mesh, Interface), Real(1e-12)));

    VolumetricPhaseConsistencyTerm phase(
        phiValue,
        phiGradient,
        phiHessian,
        Negative,
        Positive,
        phaseWeight);
    phase
      .setQuadratureOrder(CurvedQuadratureOrder)
      .setEpsilon(std::max(phaseEpsilonFactor * h, Real(1e-12)))
      .setMargin(Real(1))
      .setNormalization(std::max(totalCellMeasure(mesh), Real(1e-12)));

    report.fitEnergyInitial = fit.energy(u);
    report.phaseEnergyInitial = phase.energy(u);
    report.wrongSideInitial = phase.countWrongSideQuadrature(mesh);

    auto makeResidual = [&]()
    {
      return quality.residual(u, v)
           + deviation.residual(u, v)
           + fit.residual(u, v)
           + phase.residual(u, v);
    };
    auto energy = [&]()
    {
      return quality.energy(u)
           + deviation.energy(u)
           + fit.energy(u)
           + phase.energy(u);
    };

    IsoparametricTMOPMinimizerParameters params;
    params.maxIterations = std::max<Index>(1, maxIterations);
    params.gradientTolerance = Real(1e-10);
    params.stepTolerance = Real(1e-10);
    params.energyTolerance = Real(1e-14);
    params.preconditionerLength = std::max(Real(0.5) * h, Real(1e-8));
    params.minDetFloor = Real(0);
    auto admissible = [&]()
    {
      return isTargetAdmissible(
          mesh,
          u,
          fe,
          target,
          params.minDetFloor,
          Polytope::Type::Triangle,
          static_cast<Index>(CurvedQuadratureOrder));
    };

    const auto min = minimizeIsoparametricTMOP(
        mesh, fe, u, p, v, makeResidual, energy, admissible, Interface, params);

    report.moved = min.acceptedAnyStep;
    report.minimizerConverged = min.converged;
    report.minimizerIterations = min.iterations;
    report.acceptedSteps = min.acceptedSteps;
    report.rejectedSteps = min.rejectedSteps;
    report.initialEnergy = min.initialEnergy;
    report.finalEnergy = min.finalEnergy;
    report.initialGradientNorm = min.initialGradientNorm;
    report.finalGradientNorm = min.finalGradientNorm;
    report.finalStepNorm = min.finalStepNorm;
    report.finalAlpha = min.finalAlpha;
    report.finalDirectionalDerivative = min.finalDirectionalDerivative;
    report.fitEnergyFinal = fit.energy(u);
    report.phaseEnergyFinal = phase.energy(u);
    report.displacement = u.getData();
    report.wrongSideFinal = phase.countWrongSideQuadrature(mesh);
    return report;
  }
}

int main(int argc, char** argv)
{
  size_t resolution = 20;
  Index steps = 20;
  Real fitWeight = Real(1);
  Real phaseWeight = Real(0.1);
  Real qualityWeight = Real(1);
  Real deviationWeight = Real(1);
  Real targetBlend = Real(0.1);
  Real phaseEpsilonFactor = Real(0.5);
  Index maxMinimizerIterations = 50;

  if (argc > 1)
    resolution = static_cast<size_t>(std::max(3, std::atoi(argv[1])));
  if (argc > 2)
    steps = static_cast<Index>(std::max(2, std::atoi(argv[2])));
  if (argc > 3)
    fitWeight = std::max(Real(0), static_cast<Real>(std::atof(argv[3])));
  if (argc > 4)
    phaseWeight = std::max(Real(0), static_cast<Real>(std::atof(argv[4])));
  if (argc > 5)
    qualityWeight = std::max(Real(0), static_cast<Real>(std::atof(argv[5])));
  if (argc > 6)
    deviationWeight = std::max(Real(0), static_cast<Real>(std::atof(argv[6])));
  if (argc > 7)
    targetBlend = std::max(Real(0), std::min(Real(1),
          static_cast<Real>(std::atof(argv[7]))));
  if (argc > 8)
    phaseEpsilonFactor = std::max(Real(1e-8), static_cast<Real>(std::atof(argv[8])));
  if (argc > 9)
    maxMinimizerIterations = static_cast<Index>(std::max(1, std::atoi(argv[9])));

  const Real h = Real(1) / static_cast<Real>(resolution - 1);

  LocalMesh mesh =
    LocalMesh::UniformGrid(Polytope::Type::Triangle, { resolution, resolution });
  mesh.scale(h);
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 2);
  mesh.getConnectivity().compute(1, 0);

  P2Element fe(Polytope::Type::Triangle);
  upgradeTransformations(mesh, fe);

  IO::XDMF xdmf("LevelSetMovingCurvedCircle");
  auto grid = xdmf.grid("p2-tmop");

  std::cout << std::setprecision(12);

  for (Index step = 0; step < steps; ++step)
  {
    const Real t = static_cast<Real>(step)
      / static_cast<Real>(std::max<Index>(1, steps - 1));

    CellSpace cellSpace(mesh);
    CellField classification(cellSpace);
    const auto pre = labelAndBuildInterface(mesh, fe, t, classification);

    const auto report = solveTMOPStep(
        mesh,
        t,
        h,
        qualityWeight,
        fitWeight,
        phaseWeight,
        deviationWeight,
        targetBlend,
        phaseEpsilonFactor,
        maxMinimizerIterations,
        fe);

    const auto post = labelAndBuildInterface(mesh, fe, t, classification);
    const auto phiValue = [t](const Math::SpatialPoint& x) { return phiAt(x, t); };
    const auto curved = curvedMetrics(mesh, phiValue, Interface);
    const auto targetStats = targetQualityMetrics(
        mesh,
        ShapeSizeBlendMetric(Real(0.5)),
        CurvedQualityTargetJacobian(mesh, targetBlend),
        Polytope::Type::Triangle,
        static_cast<Index>(CurvedQuadratureOrder));

    P1<Real, LocalMesh> p1(mesh);
    GridFunction phi(p1);
    RealFunction phiFunction([t](const Point& p)
    {
      return phiAt(p.getPhysicalCoordinates(), t);
    });
    phi.project(phiFunction);

    GridFunction classificationView(p1);
    fillClassificationView(mesh, classification, classificationView);

    VectorH1<2, LocalMesh> displacementSpace(
        std::integral_constant<size_t, 2>{}, mesh, 2);
    GridFunction displacement(displacementSpace);
    if (report.displacement.size() == displacement.getData().size())
      displacement.getData() = report.displacement;
    else
      displacement.getData().setZero();

    grid.setMesh(mesh, IO::XDMF::MeshPolicy::Transient);
    grid.clear();
    grid.add("phi", phi, IO::XDMF::Center::Node);
    grid.add("classification", classificationView, IO::XDMF::Center::Node);
    grid.add("displacement", displacement, IO::XDMF::Center::Node);
    xdmf.write(t).flush();

    std::cout
      << "step=" << step
      << " t=" << t
      << " negative=" << post.cells.negative
      << " positive=" << post.cells.positive
      << " changed_cells=" << post.cells.changed
      << " interface_edges=" << post.interface.edgeCount
      << " changed_interface_edges=" << post.changedInterfaceEdges
      << " interface_max_degree=" << post.interface.maxDegree
      << " interface_branch_vertices=" << post.interface.branchVertices
      << " minimizer_converged=" << (report.minimizerConverged ? 1 : 0)
      << " minimizer_iterations=" << report.minimizerIterations
      << " accepted_steps=" << report.acceptedSteps
      << " rejected_steps=" << report.rejectedSteps
      << " energy_initial=" << report.initialEnergy
      << " energy_final=" << report.finalEnergy
      << " gradient_initial=" << report.initialGradientNorm
      << " gradient_final=" << report.finalGradientNorm
      << " step_norm_final=" << report.finalStepNorm
      << " alpha_final=" << report.finalAlpha
      << " directional_derivative_final=" << report.finalDirectionalDerivative
      << " moved=" << (report.moved ? 1 : 0)
      << " wrong_side_initial=" << report.wrongSideInitial
      << " wrong_side_final=" << report.wrongSideFinal
      << " fit_energy_initial=" << report.fitEnergyInitial
      << " fit_energy_final=" << report.fitEnergyFinal
      << " phase_energy_initial=" << report.phaseEnergyInitial
      << " phase_energy_final=" << report.phaseEnergyFinal
      << " curved_fit_rms=" << curved.fitRms
      << " curved_fit_max=" << curved.fitMax
      << " curved_qmin=" << curved.qmin
      << " curved_min_det=" << curved.minDet
      << " curved_invalid=" << curved.invalidJacobianSamples
      << " curved_overlap=" << curved.overlapSamples
      << " target_min_detJ=" << targetStats.minDetJ
      << " target_min_detT=" << targetStats.minDetT
      << " target_max_metric=" << targetStats.maxMetric
      << " pre_negative=" << pre.cells.negative
      << " pre_positive=" << pre.cells.positive
      << '\n';
  }

  xdmf.close();
  return 0;
}
