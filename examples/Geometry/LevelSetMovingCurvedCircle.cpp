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
#include <type_traits>
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
  constexpr Attribute Positive = 10;
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
    (void) theta;
    (void) t;
    return Real(0.17);
  }

  Real wavyRadiusDerivative(Real theta, Real t)
  {
    (void) theta;
    (void) t;
    return Real(0);
  }

  Real wavyRadiusSecondDerivative(Real theta, Real t)
  {
    (void) theta;
    (void) t;
    return Real(0);
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

  Math::SpatialMatrix<Real> hessianPhiAt(const Math::SpatialPoint& x, Real t)
  {
    Math::SpatialMatrix<Real> h(2, 2);
    const auto c = center(t);
    const Real dx = x[0] - c[0];
    const Real dy = x[1] - c[1];
    const Real rho2 = dx * dx + dy * dy;
    const Real rho = std::sqrt(rho2);
    if (rho <= Real(1e-14))
    {
      h(0, 0) = h(0, 1) = h(1, 0) = h(1, 1) = Real(0);
      return h;
    }

    const Real theta = std::atan2(dy, dx);
    const Real dr = wavyRadiusDerivative(theta, t);
    const Real ddr = wavyRadiusSecondDerivative(theta, t);

    const Real rho3 = rho2 * rho;
    const Real rho4 = rho2 * rho2;

    const Real hRho00 = dy * dy / rho3;
    const Real hRho01 = -dx * dy / rho3;
    const Real hRho11 = dx * dx / rho3;

    const Real hTheta00 = Real(2) * dx * dy / rho4;
    const Real hTheta01 = (dy * dy - dx * dx) / rho4;
    const Real hTheta11 = -Real(2) * dx * dy / rho4;

    const Real dThetaDx = -dy / rho2;
    const Real dThetaDy = dx / rho2;
    const Real gradTheta0 = dThetaDx;
    const Real gradTheta1 = dThetaDy;

    const Real hR00 = ddr * gradTheta0 * gradTheta0 + dr * hTheta00;
    const Real hR01 = ddr * gradTheta0 * gradTheta1 + dr * hTheta01;
    const Real hR11 = ddr * gradTheta1 * gradTheta1 + dr * hTheta11;

    h(0, 0) = hRho00 - hR00;
    h(0, 1) = hRho01 - hR01;
    h(1, 0) = h(0, 1);
    h(1, 1) = hRho11 - hR11;
    return h;
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
    return x[0] * (Real(1) - x[0]) * x[1] * (Real(1) - x[1]);
  }

  Math::SpatialMatrix<Real> boxBoundaryHessian(const Math::SpatialPoint& x)
  {
    Math::SpatialMatrix<Real> h(2, 2);
    h(0, 0) = -Real(2) * x[1] * (Real(1) - x[1]);
    h(1, 1) = -Real(2) * x[0] * (Real(1) - x[0]);
    h(0, 1) = (Real(1) - Real(2) * x[0]) * (Real(1) - Real(2) * x[1]);
    h(1, 0) = h(0, 1);
    return h;
  }

  Math::SpatialPoint boxBoundaryGradient(const Math::SpatialPoint& x)
  {
    return Math::SpatialPoint{
      (Real(1) - Real(2) * x[0]) * x[1] * (Real(1) - x[1]),
      (Real(1) - Real(2) * x[1]) * x[0] * (Real(1) - x[0])
    };
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

  Math::SpatialPoint triangleCorner(std::uint8_t k)
  {
    if (k == 1) return Math::SpatialPoint{Real(1), Real(0)};
    if (k == 2) return Math::SpatialPoint{Real(0), Real(1)};
    return Math::SpatialPoint{Real(0), Real(0)};
  }

    template <class Cell, class Edge>
    std::array<std::uint8_t, 2> localEdgeVertices(
      const Cell& cell,
      const Edge& edge)
  {
    for (std::uint8_t a = 0; a < 3; ++a)
      for (std::uint8_t b = a + 1; b < 3; ++b)
        if ((cell(a) == edge(0) && cell(b) == edge(1))
            || (cell(a) == edge(1) && cell(b) == edge(0)))
          return {{a, b}};
    return {{0, 1}};
  }

  Math::SpatialPoint triangleEdgeMidpoint(
      std::uint8_t a,
      std::uint8_t b)
  {
    const auto ca = triangleCorner(a);
    const auto cb = triangleCorner(b);
    return Math::SpatialPoint{
      Real(0.5) * (ca[0] + cb[0]),
      Real(0.5) * (ca[1] + cb[1])};
  }

  template <class Element>
  size_t scalarNodeNear(
      Geometry::Polytope::Type geometry,
      const Math::SpatialPoint& rc)
  {
    const auto& ref = Element::getNodes(geometry);
    size_t best = 0;
    Real bestD = std::numeric_limits<Real>::infinity();
    for (size_t j = 0; j < ref.size(); ++j)
    {
      Real d = 0;
      for (std::uint8_t k = 0; k < rc.size(); ++k)
      {
        const Real diff = ref[j][k] - rc[k];
        d += diff * diff;
      }
      if (d < bestD)
      {
        bestD = d;
        best = j;
      }
    }
    return best;
  }

  template <class Element>
  bool validP2CellPointCloud(
      const Geometry::PointCloud& pc,
      const Element& fe)
  {
    Geometry::ParametricTransformation<Element> trans(pc, fe);
    const auto& qf = QF::PolytopeQuadratureFormula::get(
        CurvedQuadratureOrder, Polytope::Type::Triangle);
    for (size_t q = 0; q < qf.getSize(); ++q)
    {
      Math::SpatialMatrix<Real> J;
      trans.jacobian(J, qf.getPoint(q));
      if (J.rows() != 2 || J.cols() != 2)
        return false;
      const Real det = J.determinant();
      if (!std::isfinite(det) || !(det > CurvedDetFloor))
        return false;
    }
    return true;
  }

  template <class Element>
  Geometry::PointCloud currentCellPointCloud(
      const LocalMesh& mesh,
      Index ci)
  {
    auto cellIt = mesh.getPolytope(2, ci);
    const auto& ref = Element::getNodes(Polytope::Type::Triangle);
    Geometry::PointCloud pc(2, ref.size());
    for (size_t j = 0; j < ref.size(); ++j)
    {
      const auto x = Geometry::Point(*cellIt, ref[j]).getPhysicalCoordinates();
      pc(0, j) = x[0];
      pc(1, j) = x[1];
    }
    return pc;
  }

  template <class Element>
  Geometry::PointCloud currentEdgePointCloud(
      const LocalMesh& mesh,
      Index e)
  {
    auto edgeIt = mesh.getPolytope(1, e);
    const auto& ref = Element::getNodes(Polytope::Type::Segment);
    Geometry::PointCloud pc(2, ref.size());
    for (size_t j = 0; j < ref.size(); ++j)
    {
      const auto x = Geometry::Point(*edgeIt, ref[j]).getPhysicalCoordinates();
      pc(0, j) = x[0];
      pc(1, j) = x[1];
    }
    return pc;
  }

  template <class Element>
  void initializeProjectedInterfaceTransformations(
      LocalMesh& mesh,
      const Element& fe,
      Real t)
  {
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
    const auto& conn = mesh.getConnectivity();
    const Index nc = static_cast<Index>(mesh.getCellCount());
    for (Index ci = 0; ci < nc; ++ci)
    {
      if (conn.getGeometry(2, ci) != Polytope::Type::Triangle)
        continue;
      auto affine = currentCellPointCloud<Element>(mesh, ci);
      auto desired = affine;
      const auto& cell = conn.getPolytope(2, ci);
      for (Index e : conn.getIncidence({2, 1}, ci))
      {
        const auto attr = mesh.getAttribute(1, e);
        if (!attr || *attr != Interface)
          continue;
        const auto local =
          localEdgeVertices(cell, conn.getPolytope(1, e));
        const size_t node = scalarNodeNear<Element>(
            Polytope::Type::Triangle,
            triangleEdgeMidpoint(local[0], local[1]));
        const Math::SpatialPoint x{affine(0, node), affine(1, node)};
        const auto xp = projectToInterface(x, t);
        desired(0, node) = xp[0];
        desired(1, node) = xp[1];
      }

      Real alpha = 1;
      Geometry::PointCloud blended = affine;
      for (int attempt = 0; attempt < 20; ++attempt)
      {
        for (size_t j = 0; j < desired.getCount(); ++j)
        {
          blended(0, j) = affine(0, j) + alpha * (desired(0, j) - affine(0, j));
          blended(1, j) = affine(1, j) + alpha * (desired(1, j) - affine(1, j));
        }
        if (validP2CellPointCloud(blended, fe))
          break;
        alpha *= Real(0.5);
      }
      mesh.setPolytopeTransformation(
          {size_t(2), ci},
          new ParametricTransformation<Element>(blended, fe));
    }

    Element feSeg(Polytope::Type::Segment);
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto attr = mesh.getAttribute(1, e);
      if (!attr || *attr != Interface)
        continue;
      auto pc = currentEdgePointCloud<Element>(mesh, e);
      const size_t node = scalarNodeNear<Element>(
          Polytope::Type::Segment, Math::SpatialPoint{Real(0.5)});
      const Math::SpatialPoint x{pc(0, node), pc(1, node)};
      const auto xp = projectToInterface(x, t);
      pc(0, node) = xp[0];
      pc(1, node) = xp[1];
      mesh.setPolytopeTransformation(
          {size_t(1), e},
          new ParametricTransformation<Element>(pc, feSeg));
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

  CellClassificationStats classifyCellsByAnalyticPhi(
      LocalMesh& mesh,
      Real t,
      Real epsilon,
      Real margin)
  {
    const auto before = captureCellSigns(mesh);

    CellClassificationStats stats;
    const auto& conn = mesh.getConnectivity();

    epsilon = std::max(epsilon, Real(1e-12));

    auto psi = [margin](Real z)
    {
      const Real gap = std::max(Real(0), margin - z);
      return Real(0.5) * gap * gap;
    };

    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      if (conn.getGeometry(2, c) != Polytope::Type::Triangle)
        continue;

      const auto cellIt = mesh.getPolytope(2, c);
      const auto& trans = cellIt->getTransformation();

      const auto& qf = QF::PolytopeQuadratureFormula::get(
          CurvedQuadratureOrder,
          Polytope::Type::Triangle);

      Real negativeEnergy = 0;
      Real positiveEnergy = 0;
      Real negativeMeasure = 0;
      Real positiveMeasure = 0;

      for (size_t q = 0; q < qf.getSize(); ++q)
      {
        Math::SpatialPoint x;
        trans.transform(x, qf.getPoint(q));

        Math::SpatialMatrix<Real> J;
        trans.jacobian(J, qf.getPoint(q));

        if (J.rows() != 2 || J.cols() != 2)
          continue;

        const Real det = J.determinant();
        if (!std::isfinite(det) || !(det > Real(0)))
          continue;

        const Real wdet = qf.getWeight(q) * det;
        const Real phi = phiAt(x, t);

        const Real zNegative = -phi / epsilon;
        const Real zPositive =  phi / epsilon;

        negativeEnergy += wdet * psi(zNegative);
        positiveEnergy += wdet * psi(zPositive);

        if (phi <= Real(0))
          negativeMeasure += wdet;
        else
          positiveMeasure += wdet;
      }

      Attribute next;

      if (std::isfinite(negativeEnergy)
          && std::isfinite(positiveEnergy)
          && std::abs(negativeEnergy - positiveEnergy) > Real(1e-14))
      {
        next = negativeEnergy <= positiveEnergy ? Negative : Positive;
      }
      else
      {
        // Tie fallback: majority sign by quadrature measure.
        next = negativeMeasure >= positiveMeasure ? Negative : Positive;
      }

      mesh.setAttribute({2, c}, next);

      if (next == Negative)
        ++stats.negative;
      else
        ++stats.positive;
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

  RelabelReport relabelAndRebuildInterface(
      LocalMesh& mesh,
      Real t,
      Real epsilon,
      Real margin)
  {
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 2);
    mesh.getConnectivity().compute(1, 0);
    const auto oldSkeleton = captureInterfaceSkeleton(mesh);
    clearInterfaceEdges(mesh);
    RelabelReport report;
    report.cells = classifyCellsByAnalyticPhi(mesh, t, epsilon, margin);
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
    Real initialResidual = 0;
    Real finalResidual = 0;
    Real lastQuadraticRate = -1;
    Index regularizationIncrements = 0;
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
      Index maxNewtonIterations)
  {
    SolveReport out;
    auto phiValue = [t](const Math::SpatialPoint& x) { return phiAt(x, t); };
    auto phiGradient = [t](const Math::SpatialPoint& x) { return gradPhiAt(x, t); };
    auto phiHessian = [t](const Math::SpatialPoint& x) { return hessianPhiAt(x, t); };
    auto projector = [t](const Math::SpatialPoint& x)
    {
      return projectToInterface(x, t);
    };

    Variational::RealH1Element<2> fe(Polytope::Type::Triangle);
    ++out.outerIterations;
    const auto pre = relabelAndRebuildInterface(
        mesh,
        t,
        std::max(phaseEpsilonFactor * h, Real(1e-12)),
        phaseMargin);
    if (pre.interface.edgeCount == 0)
      return out;

    upgradeTransformations(mesh, fe, Interface);
    initializeProjectedInterfaceTransformations(mesh, fe, t);
    ProjectedInterfaceTargetJacobian target(mesh, Interface, projector);
    ShapeSizeBlendMetric metric(targetBlend);

    VectorH1<2, LocalMesh> space(std::integral_constant<size_t, 2>{}, mesh, 2);
    GridFunction u(space);
    u.getData().setZero();
    TrialFunction du(space);
    TestFunction v(space);

    QualityTerm quality(metric, target, qualityWeight);
    quality.setQuadratureOrder(CurvedQuadratureOrder);
    DeviationTerm deviation(deviationWeight);
    AnalyticLevelSetFitTerm fit(
      phiValue, phiGradient, phiHessian, Optional<Attribute>(Interface), fitWeight);
    fit.setQuadratureOrder(CurvedQuadratureOrder);
    fit.setNormalization(std::max(edgeMeasure(mesh, Interface), Real(1e-12)));
    AnalyticLevelSetFitTerm bfit(
      boxBoundaryValue, boxBoundaryGradient, boxBoundaryHessian,
      Optional<Attribute>(Boundary), Real(1));
    bfit.setQuadratureOrder(CurvedQuadratureOrder);
    bfit.setNormalization(std::max(edgeMeasure(mesh, Boundary), Real(1e-12)));
    VolumetricPhaseConsistencyTerm phase(
        phiValue, phiGradient, Negative, Positive, phaseWeight);
    phase
      .setQuadratureOrder(CurvedQuadratureOrder)
      .setEpsilon(std::max(phaseEpsilonFactor * h, Real(1e-12)))
      .setMargin(phaseMargin)
      .setNormalization(std::max(totalCellMeasure(mesh), Real(1e-12)));

    out.fitEnergyInitial = fit.energy(mesh);
    out.phaseEnergyInitial = phase.energy(mesh);
    out.wrongSideInitial = phase.countWrongSideQuadrature(mesh);

    auto makeResidual = [&]()
    {
      return quality.residual(u, v)
           + deviation.residual(u, v)
         + fit.residual(u, v)
         + bfit.residual(u, v)
         + phase.residual(u, v);
    };
    auto makeTangent = [&]()
    {
      return quality.tangent(u, du, v)
           + deviation.tangent(du, v)
         + fit.tangent(u, du, v)
         + bfit.tangent(u, du, v)
         + phase.tangent(u, du, v);
    };
    auto energy = [&]()
    {
      return quality.energy(u)
           + deviation.energy(u)
         + fit.energy(u)
         + bfit.energy(u)
         + phase.energy(u);
    };

    IsoparametricTMOPParameters params;
    params.maxIterations = maxNewtonIterations;
    params.residualTolerance = Real(1e-10);
    params.stepTolerance = Real(1e-10);
    params.minDetFloor = Real(0);
    params.printIterations = false;
    const auto tmop = solveIsoparametricTMOP(
      mesh,
      fe,
      u,
      du,
      v,
      makeResidual,
      makeTangent,
      energy,
      Interface,
      params);

    out.newtonIterations = tmop.iterations;
    out.initialResidual = tmop.initialResidual;
    out.finalResidual = tmop.finalResidual;
    out.lastQuadraticRate = tmop.lastQuadraticRate;
    out.lastAlpha = tmop.lastAcceptedAlpha;
    const bool residualImproved =
      std::isfinite(tmop.initialResidual)
      && std::isfinite(tmop.finalResidual)
      && (tmop.finalResidual <= tmop.initialResidual);
    out.accepted = tmop.converged || (tmop.acceptedStep && residualImproved);

    if (!out.accepted)
    {
      Real previousResidual = -1;
      bool hadAcceptedStep = false;
      for (Index it = 0; it < maxNewtonIterations; ++it)
      {
        ++out.newtonIterations;
        LinearForm R(v);
        R = makeResidual();
        R.assemble();
        const auto r = R.getVector();
        if (!r.allFinite())
          break;
        const Real rNorm = r.norm();
        if (it == 0)
          out.initialResidual = rNorm;
        out.finalResidual = rNorm;
        if (previousResidual > Real(0) && std::isfinite(previousResidual))
          out.lastQuadraticRate = rNorm / (previousResidual * previousResidual);
        if (rNorm <= Real(1e-10))
          break;

        BilinearForm J(du, v);
        J = makeTangent();
        J.assemble();

        const auto A = J.getOperator();
        const auto rhs = -r;
        const Math::Vector<Real> u0 = u.getData();
        const Real e0 = energy();
        Eigen::SparseLU<std::decay_t<decltype(A)>> lu;
        lu.compute(A);
        if (lu.info() != Eigen::Success)
          break;

        const auto dx = lu.solve(rhs);
        if (lu.info() != Eigen::Success || !dx.allFinite())
          break;

        bool acceptedStep = false;
        Real alphaNewton = 1;
        Real stepNorm = 0;
        for (int attempt = 0; attempt < 20; ++attempt)
        {
          u.getData() = u0 + alphaNewton * dx;
          if (!u.getData().allFinite())
          {
            alphaNewton *= Real(0.5);
            continue;
          }
          const Real eTrial = energy();
          LinearForm RTrial(v);
          RTrial = makeResidual();
          RTrial.assemble();
          const auto rTrial = RTrial.getVector();
          if (!rTrial.allFinite())
          {
            alphaNewton *= Real(0.5);
            continue;
          }
          const Real rTrialNorm = rTrial.norm();
          const bool residualDecrease =
            rTrialNorm <= (Real(1) - Real(1e-4) * alphaNewton) * rNorm;
          if (std::isfinite(eTrial) && eTrial <= e0 && residualDecrease)
          {
            stepNorm = (alphaNewton * dx).norm();
            out.finalResidual = rTrialNorm;
            acceptedStep = true;
            break;
          }
          alphaNewton *= Real(0.5);
        }

        if (!acceptedStep)
        {
          u.getData() = u0;
          break;
        }

        hadAcceptedStep = true;
        out.lastAlpha = alphaNewton;
        previousResidual = rNorm;
        if (stepNorm <= Real(1e-10))
          break;
      }

      if (!hadAcceptedStep || !u.getData().allFinite())
        return out;

      LocalMesh beforeMove(mesh);
      Real alpha = 1;
      bool acceptedMove = false;
      for (int attempt = 0; attempt < 20; ++attempt)
      {
        LocalMesh candidate(beforeMove);
        VectorH1<2, LocalMesh> candidateSpace(std::integral_constant<size_t, 2>{}, candidate, 2);
        GridFunction scaledU(candidateSpace);
        scaledU.getData() = u.getData();
        scaledU.getData() *= alpha;
        moveMesh(candidate, scaledU, fe, Interface);

        out.curved = curvedMetrics(candidate, phiValue, Interface);
        const bool curvedValid =
          out.curved.invalidJacobianSamples == 0
          && out.curved.overlapSamples == 0
          && std::isfinite(out.curved.minDet)
          && std::isfinite(out.curved.qmin);
        if (curvedValid)
        {
          mesh = std::move(candidate);
          acceptedMove = true;
          break;
        }
        alpha *= Real(0.5);
      }

      if (!acceptedMove)
      {
        mesh = std::move(beforeMove);
        return out;
      }

      out.accepted = true;
      out.lastAlpha = alpha;
    }

      mesh.getConnectivity().compute(2, 1);
      mesh.getConnectivity().compute(1, 2);
      mesh.getConnectivity().compute(1, 0);

      const auto post = relabelAndRebuildInterface(
          mesh,
          t,
          std::max(phaseEpsilonFactor * h, Real(1e-12)),
          phaseMargin);
        initializeProjectedInterfaceTransformations(mesh, fe, t);
        out.fitEnergyFinal = fit.energy(mesh);
        out.phaseEnergyFinal = phase.energy(mesh);
        out.wrongSideFinal = phase.countWrongSideQuadrature(mesh);
      out.changedCells += post.cells.changed;
      out.changedInterfaceEdges += post.changedInterfaceEdges;
    return out;
  }
}

int main(int argc, char** argv)
{
  size_t resolution = 20;
  Index steps = 20;
  Real fitWeight = Real(30);
  Real qualityWeight = Real(0.03);
  Real deviationWeight = Real(0.25);
  Real targetBlend = Real(0.5);
  Real phaseWeight = Real(0);
  Real phaseEpsilonFactor = Real(0.5);
  Real phaseMargin = Real(1);
  Index tmopMaxIterations = 20;

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
    relabelAndRebuildInterface(
        mesh,
        t,
        std::max(phaseEpsilonFactor * h, Real(1e-12)),
        phaseMargin);

    const auto report = solveNoCutTMOP(
        mesh, t, h,
        qualityWeight, fitWeight, phaseWeight, deviationWeight,
      targetBlend, phaseEpsilonFactor, phaseMargin,
      tmopMaxIterations);

    relabelAndRebuildInterface(
        mesh,
        t,
        std::max(phaseEpsilonFactor * h, Real(1e-12)),
        phaseMargin);
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
              << " residual_initial=" << report.initialResidual
              << " residual_final=" << report.finalResidual
              << " qrate_last=" << report.lastQuadraticRate
              << std::endl;

clearInterfaceEdges(mesh);
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
