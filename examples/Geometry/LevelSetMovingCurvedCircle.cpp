/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */

/**
 * @example No-cut P2 TMOP fitting of a moving curved level-set circle
 *
 * Pure Rodin no-cut pipeline:
 *   classify cells by analytic phi
 *   -> mark attribute-jump edges as Interface
 *   -> P2 TMOP with quality + phase consistency + edge fit + deviation
 *   -> reclassify/rebuild Interface
 *   -> carry the moved mesh forward.
 *
 * There is no Rodin cut and no MMG step in this example. Interface edges are
 * material-boundary edges, i.e. edges between Negative and Positive cells.
 * They are never raw sign-changing/straddling edges.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

#include <Eigen/SparseLU>

#include <Rodin/Adaptation.h>
#include <Rodin/Adaptation/TargetMatrixOptimization/IsoparametricGeometry.h>
#include <Rodin/Adaptation/TargetMatrixOptimization/Metrics.h>
#include <Rodin/Assembly/Default.h>
#include <Rodin/Geometry.h>
#include <Rodin/IO.h>
#include <Rodin/QF/PolytopeQuadratureFormula.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Adaptation::TargetMatrixOptimization;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

namespace
{
  constexpr Real Pi = 3.14159265358979323846;
  constexpr Attribute Interface = 30;
  constexpr Attribute Boundary = 40;
  constexpr Attribute Negative = 1;
  constexpr Attribute Positive = 2;
  constexpr Index CurvedQuadratureOrder = 4;
  constexpr Real CurvedDetFloor = 1e-12;

  Math::SpatialPoint center(Real t)
  {
    return Math::SpatialPoint{
      Real(0.5) + Real(0.18) * std::sin(Real(2) * Pi * t),
      Real(0.5) + Real(0.13) * std::sin(Real(4) * Pi * t + Real(0.35))
    };
  }

  Real wavyRadius(Real theta, Real t)
  {
    return Real(0.17)
      + Real(0.035) * std::sin(
          Real(5) * theta + Real(0.65) * std::sin(Real(2) * Pi * t));
  }

  Real wavyRadiusDerivative(Real theta, Real t)
  {
    return Real(0.175) * std::cos(
        Real(5) * theta + Real(0.65) * std::sin(Real(2) * Pi * t));
  }

  Real phiAt(const Math::SpatialPoint& x, Real t)
  {
    const auto c = center(t);
    const Real dx = x[0] - c[0];
    const Real dy = x[1] - c[1];
    const Real theta = std::atan2(dy, dx);
    return std::hypot(dx, dy) - wavyRadius(theta, t);
  }

  Math::SpatialPoint gradPhiAt(const Math::SpatialPoint& x, Real t)
  {
    const auto c = center(t);
    const Real dx = x[0] - c[0];
    const Real dy = x[1] - c[1];
    const Real rho2 = dx * dx + dy * dy;
    const Real rho = std::sqrt(rho2);
    if (rho <= Real(1e-14))
      return Math::SpatialPoint{1, 0};
    const Real theta = std::atan2(dy, dx);
    const Real dr = wavyRadiusDerivative(theta, t);
    return Math::SpatialPoint{
      dx / rho + dr * dy / rho2,
      dy / rho - dr * dx / rho2
    };
  }

  Math::SpatialPoint projectToInterface(const Math::SpatialPoint& x, Real t)
  {
    Math::SpatialPoint p = x;
    for (int i = 0; i < 3; ++i)
    {
      const Real f = phiAt(p, t);
      const auto g = gradPhiAt(p, t);
      const Real gg = g[0] * g[0] + g[1] * g[1];
      if (gg < Real(1e-30))
        break;
      p = Math::SpatialPoint{p[0] - f * g[0] / gg,
                             p[1] - f * g[1] / gg};
    }
    return p;
  }

  Real boxBoundaryValue(const Math::SpatialPoint& x)
  {
    const std::array<Real, 4> d = {{
      std::abs(x[0]), std::abs(x[0] - Real(1)),
      std::abs(x[1]), std::abs(x[1] - Real(1)) }};
    const auto side =
      static_cast<size_t>(std::min_element(d.begin(), d.end()) - d.begin());
    if (side == 0) return x[0];
    if (side == 1) return x[0] - Real(1);
    if (side == 2) return x[1];
    return x[1] - Real(1);
  }

  Math::SpatialPoint boxBoundaryGradient(const Math::SpatialPoint& x)
  {
    const std::array<Real, 4> d = {{
      std::abs(x[0]), std::abs(x[0] - Real(1)),
      std::abs(x[1]), std::abs(x[1] - Real(1)) }};
    const auto side =
      static_cast<size_t>(std::min_element(d.begin(), d.end()) - d.begin());
    if (side < 2)
      return Math::SpatialPoint{1, 0};
    return Math::SpatialPoint{0, 1};
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
        mesh.setAttribute({1, e}, Boundary);
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
        mesh.setAttribute({1, e}, Optional<Attribute>{});
    }
  }

  std::vector<int> captureCellSigns(const LocalMesh& mesh)
  {
    std::vector<int> out(mesh.getCellCount(), 0);
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      const auto attr = mesh.getAttribute(2, c);
      if (!attr) continue;
      if (*attr == Negative) out[static_cast<size_t>(c)] = -1;
      else if (*attr == Positive) out[static_cast<size_t>(c)] = 1;
    }
    return out;
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

  Index countChangedCellSigns(
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

  struct CellClassificationStats
  {
    Index negative = 0;
    Index positive = 0;
    Index changed = 0;
  };

  CellClassificationStats classifyCellsByAnalyticPhi(LocalMesh& mesh, Real t)
  {
    const auto before = captureCellSigns(mesh);
    CellClassificationStats stats;
    const auto& conn = mesh.getConnectivity();
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
        continue;
      const auto cellIt = mesh.getPolytope(2, c);
      Math::SpatialPoint x;
      cellIt->getTransformation().transform(
          x, Math::SpatialPoint{Real(1) / Real(3), Real(1) / Real(3)});
      const Attribute next = phiAt(x, t) <= Real(0) ? Negative : Positive;
      mesh.setAttribute({2, c}, next);
      if (next == Negative) ++stats.negative;
      else ++stats.positive;
    }
    stats.changed = countChangedCellSigns(before, captureCellSigns(mesh));
    return stats;
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
      if (d > 2) ++stats.branchVertices;
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
      const auto edgeAttr = mesh.getAttribute(1, e);
      if (edgeAttr && *edgeAttr == Boundary)
        continue;
      const auto adj = conn.getIncidence({1, 2}, e);
      if (adj.size() != 2)
        continue;
      const auto a0 = mesh.getAttribute(2, adj[0]);
      const auto a1 = mesh.getAttribute(2, adj[1]);
      if (!a0 || !a1)
        continue;
      const bool jump = (*a0 == Negative && *a1 == Positive)
                     || (*a0 == Positive && *a1 == Negative);
      if (jump)
        mesh.setAttribute({1, e}, Interface);
    }
    return computeInterfaceStats(mesh);
  }

  struct RelabelReport
  {
    CellClassificationStats cells;
    InterfaceStats interface;
    Index changedInterfaceEdges = 0;
  };

  RelabelReport relabelAndRebuildInterface(LocalMesh& mesh, Real t)
  {
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 2);
    mesh.getConnectivity().compute(1, 0);
    const auto oldSkeleton = captureInterfaceSkeleton(mesh);
    clearInterfaceEdges(mesh);
    RelabelReport report;
    report.cells = classifyCellsByAnalyticPhi(mesh, t);
    report.interface = markAttributeJumpInterfaceEdges(mesh);
    report.changedInterfaceEdges =
      countChangedSkeleton(oldSkeleton, captureInterfaceSkeleton(mesh));
    return report;
  }

  Real totalCellMeasure(LocalMesh& mesh)
  {
    Real measure = 0;
    const auto& conn = mesh.getConnectivity();
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
        continue;
      const auto cellIt = mesh.getPolytope(2, c);
      const auto& qf = QF::PolytopeQuadratureFormula::get(
          CurvedQuadratureOrder, Polytope::Type::Triangle);
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
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto attr = mesh.getAttribute(1, e);
      if (!attr || *attr != attribute)
        continue;
      const auto edgeIt = mesh.getPolytope(1, e);
      for (Real s : {Real(0.25), Real(0.75)})
      {
        Math::SpatialMatrix<Real> J;
        edgeIt->getTransformation().jacobian(J, Math::SpatialPoint{s});
        if (J.rows() == 2 && J.cols() == 1)
          measure += Real(0.5) * std::hypot(J(0, 0), J(1, 0));
      }
    }
    return measure;
  }

  struct LinearQualityStats
  {
    Real qmin = std::numeric_limits<Real>::infinity();
    Index inverted = 0;
  };

  Real triangleQuality(
      const Math::SpatialPoint& x0,
      const Math::SpatialPoint& x1,
      const Math::SpatialPoint& x2)
  {
    const Real orient =
      (x1[0] - x0[0]) * (x2[1] - x0[1])
    - (x1[1] - x0[1]) * (x2[0] - x0[0]);
    if (!(orient > Real(0)) || !std::isfinite(orient))
      return 0;
    const Real area = Real(0.5) * orient;
    const Real den = (x1 - x0).squaredNorm()
                   + (x2 - x1).squaredNorm()
                   + (x0 - x2).squaredNorm();
    if (!(den > Real(0)) || !std::isfinite(den))
      return 0;
    return Real(4) * std::sqrt(Real(3)) * area / den;
  }

  LinearQualityStats linearQuality(LocalMesh& mesh)
  {
    LinearQualityStats stats;
    const auto& conn = mesh.getConnectivity();
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
        continue;
      const auto& cell = conn.getPolytope(2, c);
      const auto x0 = mesh.getVertexCoordinates(cell(0));
      const auto x1 = mesh.getVertexCoordinates(cell(1));
      const auto x2 = mesh.getVertexCoordinates(cell(2));
      const Real orient =
        (x1[0] - x0[0]) * (x2[1] - x0[1])
      - (x1[1] - x0[1]) * (x2[0] - x0[0]);
      if (!(orient > Real(0)))
        ++stats.inverted;
      stats.qmin = std::min(stats.qmin, triangleQuality(x0, x1, x2));
    }
    if (!std::isfinite(stats.qmin)) stats.qmin = 0;
    return stats;
  }

  struct SolveReport
  {
    bool accepted = false;
    Index outerIterations = 0;
    Index newtonIterations = 0;
    Index changedCells = 0;
    Index changedInterfaceEdges = 0;
    Index wrongSideInitial = 0;
    Index wrongSideFinal = 0;
    Real lastAlpha = 0;
    Real fitEnergyInitial = 0;
    Real fitEnergyFinal = 0;
    Real phaseEnergyInitial = 0;
    Real phaseEnergyFinal = 0;
    CurvedMetrics curved;
  };

  SolveReport solveNoCutTMOP(
      LocalMesh& mesh,
      Real t,
      Real h,
      Real qualityWeight,
      Real fitWeight,
      Real phaseWeight,
      Real deviationWeight,
      Real targetBlend,
      Real phaseEpsilonFactor,
      Real phaseMargin,
      Index maxOuterIterations,
      Index maxNewtonIterations)
  {
    SolveReport out;
    auto phiValue = [t](const Math::SpatialPoint& x) { return phiAt(x, t); };
    auto phiGradient = [t](const Math::SpatialPoint& x) { return gradPhiAt(x, t); };
    auto projector = [t](const Math::SpatialPoint& x)
    {
      return projectToInterface(x, t);
    };

    Variational::RealH1Element<2> fe(Polytope::Type::Triangle);
    for (Index outer = 0; outer < maxOuterIterations; ++outer)
    {
      ++out.outerIterations;
      const auto pre = relabelAndRebuildInterface(mesh, t);
      if (pre.interface.edgeCount == 0)
        break;

      upgradeTransformations(mesh, fe, Interface);
      auto target = ProjectedQualityTargetJacobian(
          mesh, Interface, projector, targetBlend);
      ShapeSizeBlendMetric metric(Real(0.5));

      VectorH1<2, LocalMesh> space(std::integral_constant<size_t, 2>{}, mesh, 2);
      GridFunction u(space);
      u.getData().setZero();
      TrialFunction du(space);
      TestFunction v(space);

      QualityTerm quality(metric, target, qualityWeight);
      quality.setQuadratureOrder(CurvedQuadratureOrder);
      DeviationTerm deviation(deviationWeight);
      AnalyticLevelSetFitTerm fit(
          phiValue, phiGradient, Optional<Attribute>(Interface), fitWeight);
      fit.setNormalization(std::max(edgeMeasure(mesh, Interface), Real(1e-12)));
      AnalyticLevelSetFitTerm boundaryFit(
          boxBoundaryValue, boxBoundaryGradient,
          Optional<Attribute>(Boundary), Real(1));
      boundaryFit.setNormalization(std::max(edgeMeasure(mesh, Boundary), Real(1e-12)));
      VolumetricPhaseConsistencyTerm phase(
          phiValue, phiGradient, Negative, Positive, phaseWeight);
      phase
        .setQuadratureOrder(CurvedQuadratureOrder)
        .setEpsilon(std::max(phaseEpsilonFactor * h, Real(1e-12)))
        .setMargin(phaseMargin)
        .setNormalization(std::max(totalCellMeasure(mesh), Real(1e-12)));

      if (outer == 0)
      {
        out.fitEnergyInitial = fit.energy(mesh);
        out.phaseEnergyInitial = phase.energy(mesh);
        out.wrongSideInitial = phase.countWrongSideQuadrature(mesh);
      }

      auto makeResidual = [&]()
      {
        return quality.residual(u, v)
             + deviation.residual(u, v)
             + fit.residual(u, v)
             + boundaryFit.residual(u, v)
             + phase.residual(u, v);
      };
      auto makeTangent = [&]()
      {
        return quality.tangent(u, du, v)
             + deviation.tangent(du, v)
             + fit.tangent(u, du, v)
             + boundaryFit.tangent(u, du, v)
             + phase.tangent(u, du, v);
      };
      auto energy = [&]()
      {
        return quality.energy(u)
             + deviation.energy(u)
             + fit.energy(u)
             + boundaryFit.energy(u)
             + phase.energy(u);
      };

      bool hadAcceptedStep = false;
      try
      {
        for (Index it = 0; it < maxNewtonIterations; ++it)
        {
          ++out.newtonIterations;
          LinearForm R(v);
          R = makeResidual();
          R.assemble();
          const auto r = R.getVector();
          if (!r.allFinite())
            break;
          if (r.norm() <= Real(1e-10))
            break;

          BilinearForm J(du, v);
          J = makeTangent();
          J.assemble();

          Eigen::SparseLU<std::decay_t<decltype(J.getOperator())>> lu;
          lu.compute(J.getOperator());
          if (lu.info() != Eigen::Success)
            break;
          const Math::Vector<Real> dx = lu.solve(-r);
          if (lu.info() != Eigen::Success || !dx.allFinite())
            break;

          const Math::Vector<Real> u0 = u.getData();
          const Real e0 = energy();
          Real alpha = 1;
          bool acceptedStep = false;
          for (int ls = 0; ls < 30; ++ls)
          {
            u.getData() = u0 + alpha * dx;
            const Real e = energy();
            if (std::isfinite(e) && e <= e0 * (Real(1) + Real(1e-12))
                && isCurvedMoveValid(mesh, u, fe, CurvedDetFloor))
            {
              acceptedStep = true;
              break;
            }
            alpha *= Real(0.5);
          }
          out.lastAlpha = alpha;
          if (!acceptedStep)
          {
            u.getData() = u0;
            break;
          }
          hadAcceptedStep = true;
          if (alpha * dx.norm() <= Real(1e-10))
            break;
        }
      }
      catch (...)
      {
        hadAcceptedStep = false;
      }

      if (!hadAcceptedStep || !u.getData().allFinite())
        break;

      LocalMesh beforeMove(mesh);
      moveMesh(mesh, u, fe, Interface);
      out.curved = curvedMetrics(mesh, phiValue, Interface);
      const auto lin = linearQuality(mesh);
      const bool valid = out.curved.invalidJacobianSamples == 0
        && out.curved.overlapSamples == 0
        && out.curved.minDet > CurvedDetFloor
        && out.curved.qmin > Real(0)
        && lin.inverted == 0;
      if (!valid)
      {
        mesh = std::move(beforeMove);
        break;
      }

      out.accepted = true;
      out.fitEnergyFinal = fit.energy(mesh);
      out.phaseEnergyFinal = phase.energy(mesh);
      out.wrongSideFinal = phase.countWrongSideQuadrature(mesh);

      syncLinearBackbone(mesh, fe);
      demoteTransformations(mesh);
      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 2);
      mesh.getConnectivity().compute(1, 0);

      const auto post = relabelAndRebuildInterface(mesh, t);
      out.changedCells += post.cells.changed;
      out.changedInterfaceEdges += post.changedInterfaceEdges;
      if (post.cells.changed == 0 && post.changedInterfaceEdges == 0)
        break;
    }
    return out;
  }

  LocalMesh interfacePolyline(const LocalMesh& mesh)
  {
    std::vector<Math::SpatialPoint> points;
    std::vector<std::array<Index, 2>> edges;
    const auto& conn = mesh.getConnectivity();
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto attr = mesh.getAttribute(1, e);
      if (!attr || *attr != Interface)
        continue;
      const auto& edge = conn.getPolytope(1, e);
      const Index base = static_cast<Index>(points.size());
      points.push_back(mesh.getVertexCoordinates(edge(0)));
      points.push_back(mesh.getVertexCoordinates(edge(1)));
      edges.push_back({{base, base + 1}});
    }
    auto builder = LocalMesh::Builder();
    builder.initialize(2).nodes(points.size()).reserve(1, edges.size());
    for (const auto& p : points)
      builder.vertex(p);
    for (const auto& e : edges)
      builder.polytope(Polytope::Type::Segment, {e[0], e[1]});
    auto out = builder.finalize();
    for (Index e = 0; e < out.getPolytopeCount(1); ++e)
      out.setAttribute({1, e}, Interface);
    return out;
  }
}

int main(int argc, char** argv)
{
  size_t resolution = 40;
  Index steps = 20;
  Real fitWeight = 2;
  Real qualityWeight = 1;
  Real deviationWeight = Real(0.25);
  Real targetBlend = Real(0.05);
  Real phaseWeight = Real(0.05);
  Real phaseEpsilonFactor = Real(0.5);
  Real phaseMargin = Real(1);
  Index tmopMaxIterations = 8;
  Index maxOuterRelabelIterations = 4;

  if (argc > 1)
    resolution = static_cast<size_t>(std::max(3, std::atoi(argv[1])));
  if (argc > 2)
    steps = static_cast<Index>(std::max(2, std::atoi(argv[2])));
  if (argc > 3)
    fitWeight = std::max(Real(0), static_cast<Real>(std::atof(argv[3])));
  if (argc > 4)
    qualityWeight = std::max(Real(0), static_cast<Real>(std::atof(argv[4])));
  if (argc > 5)
    deviationWeight = std::max(Real(0), static_cast<Real>(std::atof(argv[5])));
  if (argc > 6)
    targetBlend = std::max(Real(0), std::min(Real(1),
          static_cast<Real>(std::atof(argv[6]))));
  if (argc > 7)
    phaseWeight = std::max(Real(0), static_cast<Real>(std::atof(argv[7])));
  if (argc > 8)
    phaseEpsilonFactor = std::max(Real(1e-8), static_cast<Real>(std::atof(argv[8])));
  if (argc > 9)
    phaseMargin = static_cast<Real>(std::atof(argv[9]));
  if (argc > 10)
    tmopMaxIterations = static_cast<Index>(std::max(1, std::atoi(argv[10])));
  if (argc > 11)
    maxOuterRelabelIterations = static_cast<Index>(std::max(1, std::atoi(argv[11])));

  const Real h = Real(1) / static_cast<Real>(resolution - 1);

  LocalMesh background =
    LocalMesh::UniformGrid(Polytope::Type::Triangle, { resolution, resolution });
  background.scale(h);
  background.getConnectivity().compute(2, 1);
  background.getConnectivity().compute(1, 2);
  background.getConnectivity().compute(1, 0);

  IO::XDMF xdmf("LevelSetMovingCurvedCircle");

  Real sumFit = 0;
  Real sumQmin = 0;
  Real sumPhaseInitial = 0;
  Real sumPhaseFinal = 0;
  Index sumAccepted = 0;
  Index sumChangedCells = 0;
  Index sumChangedEdges = 0;
  Index sumWrongSideInitial = 0;
  Index sumWrongSideFinal = 0;

  for (Index step = 0; step < steps; ++step)
  {
    const Real t = static_cast<Real>(step)
      / static_cast<Real>(std::max<Index>(1, steps - 1));

    LocalMesh mesh(background);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 2);
    mesh.getConnectivity().compute(1, 0);
    annotateBoundary(mesh);
    relabelAndRebuildInterface(mesh, t);

    const auto report = solveNoCutTMOP(
        mesh, t, h,
        qualityWeight, fitWeight, phaseWeight, deviationWeight,
        targetBlend, phaseEpsilonFactor, phaseMargin,
        maxOuterRelabelIterations, tmopMaxIterations);

    relabelAndRebuildInterface(mesh, t);
    const auto interfaceStats = computeInterfaceStats(mesh);
    const auto lin = linearQuality(mesh);
    auto phiValue = [t](const Math::SpatialPoint& x) { return phiAt(x, t); };
    const auto curved = report.accepted
      ? report.curved
      : curvedMetrics(mesh, phiValue, Interface);

    P1<Real, LocalMesh> outputSpace(mesh);
    GridFunction outputPhi(outputSpace);
    for (Index i = 0; i < mesh.getVertexCount(); ++i)
      outputPhi[i] = phiAt(mesh.getVertexCoordinates(i), t);

    auto grid = xdmf.grid("no-cut-p2-tmop-relabel");
    grid.setMesh(mesh, IO::XDMF::MeshPolicy::Transient);
    grid.clear();
    grid.add("phi", outputPhi, IO::XDMF::Center::Node);

    auto iface = interfacePolyline(mesh);
    auto curveGrid = xdmf.grid("attribute-jump-interface");
    curveGrid.setMesh(iface, IO::XDMF::MeshPolicy::Transient);
    curveGrid.clear();
    xdmf.write(t).flush();

    sumFit += curved.fitRms;
    sumQmin += lin.qmin;
    sumPhaseInitial += report.phaseEnergyInitial;
    sumPhaseFinal += report.phaseEnergyFinal;
    sumAccepted += report.accepted ? 1 : 0;
    sumChangedCells += report.changedCells;
    sumChangedEdges += report.changedInterfaceEdges;
    sumWrongSideInitial += report.wrongSideInitial;
    sumWrongSideFinal += report.wrongSideFinal;

    std::cout << "step " << step
              << " t=" << t
              << " cells=" << mesh.getCellCount()
              << " interface_edges=" << interfaceStats.edgeCount
              << " interface_max_degree=" << interfaceStats.maxDegree
              << " interface_branch_vertices=" << interfaceStats.branchVertices
              << " qmin=" << lin.qmin
              << " inverted=" << lin.inverted
              << " curved_fit_rms=" << curved.fitRms
              << " curved_fit_max=" << curved.fitMax
              << " curved_qmin=" << curved.qmin
              << " curved_min_det=" << curved.minDet
              << " curved_invalid=" << curved.invalidJacobianSamples
              << " curved_overlap=" << curved.overlapSamples
              << " outer_iterations=" << report.outerIterations
              << " newton_iterations=" << report.newtonIterations
              << " tmop_accepted=" << (report.accepted ? 1 : 0)
              << " last_alpha=" << report.lastAlpha
              << " changed_cells=" << report.changedCells
              << " changed_interface_edges=" << report.changedInterfaceEdges
              << " wrong_side_initial=" << report.wrongSideInitial
              << " wrong_side_final=" << report.wrongSideFinal
              << " fit_energy_initial=" << report.fitEnergyInitial
              << " fit_energy_final=" << report.fitEnergyFinal
              << " phase_energy_initial=" << report.phaseEnergyInitial
              << " phase_energy_final=" << report.phaseEnergyFinal
              << " fit_weight=" << fitWeight
              << " phase_weight=" << phaseWeight
              << " quality_weight=" << qualityWeight
              << " deviation_weight=" << deviationWeight
              << " target_blend=" << targetBlend
              << std::endl;

    clearInterfaceEdges(mesh);
    mesh.flush();
    background = std::move(mesh);
  }

  xdmf.close();

  const Real n = static_cast<Real>(std::max<Index>(steps, 1));
  std::cout << "AVERAGE over " << steps << " steps:"
            << " accepted=" << static_cast<Real>(sumAccepted) / n
            << " qmin=" << sumQmin / n
            << " curved_fit_rms=" << sumFit / n
            << " changed_cells=" << static_cast<Real>(sumChangedCells) / n
            << " changed_interface_edges=" << static_cast<Real>(sumChangedEdges) / n
            << " wrong_side_initial=" << static_cast<Real>(sumWrongSideInitial) / n
            << " wrong_side_final=" << static_cast<Real>(sumWrongSideFinal) / n
            << " phase_energy_initial=" << sumPhaseInitial / n
            << " phase_energy_final=" << sumPhaseFinal / n
            << std::endl;

  return 0;
}
