/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

#include <Rodin/Adaptation.h>
#include <Rodin/Adaptation/TargetMatrixOptimization/IsoparametricGeometry.h>
#include <Rodin/Adaptation/TargetMatrixOptimization/Metrics.h>
#include <Rodin/Geometry.h>
#include <Rodin/Geometry/LevelSetDiscretizerTriangles.h>
#include <Rodin/Geometry/TriangleMeshOptimizer.h>
#include <Rodin/IO.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Adaptation::TargetMatrixOptimization;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

// =====================================================================
//  Analytic moving level set + box boundary
// =====================================================================
namespace
{
  constexpr Real Pi = 3.14159265358979323846;
  constexpr Attribute Interface = 30;
  constexpr Attribute Boundary  = 40;
  constexpr Attribute Negative  = 1;
  constexpr Attribute Positive  = 2;

  Math::SpatialPoint center(Real t)
  {
    return Math::SpatialPoint{
      Real(0.5) + Real(0.18) * std::sin(Real(2) * Pi * t),
      Real(0.5) + Real(0.13) * std::sin(Real(4) * Pi * t + Real(0.35)) };
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
    if (rho <= Real(1e-14)) return Math::SpatialPoint{1, 0};
    const Real theta = std::atan2(dy, dx);
    const Real dr = wavyRadiusDerivative(theta, t);
    return Math::SpatialPoint{
      dx / rho + dr * dy / rho2,
      dy / rho - dr * dx / rho2 };
  }

  Real boxValue(const Math::SpatialPoint& x)
  {
    const std::array<Real, 4> d = {{
      std::abs(x[0]), std::abs(x[0] - Real(1)),
      std::abs(x[1]), std::abs(x[1] - Real(1)) }};
    const auto side = static_cast<size_t>(
        std::min_element(d.begin(), d.end()) - d.begin());
    if (side == 0) return x[0];
    if (side == 1) return x[0] - Real(1);
    if (side == 2) return x[1];
    return x[1] - Real(1);
  }

  Math::SpatialPoint boxGradient(const Math::SpatialPoint& x)
  {
    const std::array<Real, 4> d = {{
      std::abs(x[0]), std::abs(x[0] - Real(1)),
      std::abs(x[1]), std::abs(x[1] - Real(1)) }};
    const auto side = static_cast<size_t>(
        std::min_element(d.begin(), d.end()) - d.begin());
    return side < 2
      ? Math::SpatialPoint{1, 0}
      : Math::SpatialPoint{0, 1};
  }
}

// =====================================================================
//  Pipeline-level utilities. None of these are TMOP-specific.
// =====================================================================
namespace
{
  void annotateBoxBoundary(LocalMesh& mesh)
  {
    const auto& conn = mesh.getConnectivity();
    const Real eps = Real(1e-12);
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto& edge = conn.getPolytope(1, e);
      const auto x0 = mesh.getVertexCoordinates(edge(0));
      const auto x1 = mesh.getVertexCoordinates(edge(1));
      auto onBox = [&](const Math::SpatialPoint& x)
      {
        return std::abs(x[0]) <= eps || std::abs(x[0] - Real(1)) <= eps
            || std::abs(x[1]) <= eps || std::abs(x[1] - Real(1)) <= eps;
      };
      if (onBox(x0) && onBox(x1))
        mesh.setAttribute({ 1, e }, Boundary);
    }
  }

  std::vector<char> interfaceVertexMask(const LocalMesh& mesh)
  {
    std::vector<char> mask(mesh.getVertexCount(), 0);
    const auto& conn = mesh.getConnectivity();
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto attr = mesh.getAttribute(1, e);
      if (!attr || *attr != Interface) continue;
      const auto& edge = conn.getPolytope(1, e);
      mask[edge(0)] = 1;
      mask[edge(1)] = 1;
    }
    return mask;
  }

  // Topological flood-fill side classification. Seeded by the deepest phi
  // value, propagated cell-to-cell, flipping side only across interface
  // edges. Independent of per-cell phi sampling -> correct on curved
  // interfaces.
  template <class PhiFn>
  Index classifySides(LocalMesh& mesh, PhiFn phi)
  {
    auto& conn = mesh.getConnectivity();
    conn.compute(2, 1);
    conn.compute(1, 2);
    const Index nc = mesh.getCellCount();
    std::vector<char> side(static_cast<size_t>(nc), 0);

    auto isIface = [&](Index e)
    {
      const auto a = mesh.getAttribute(1, e);
      return a && *a == Interface;
    };
    auto centroidPhi = [&](Index c)
    {
      const auto& cell = conn.getPolytope(2, c);
      const auto x0 = mesh.getVertexCoordinates(cell(0));
      const auto x1 = mesh.getVertexCoordinates(cell(1));
      const auto x2 = mesh.getVertexCoordinates(cell(2));
      return phi(Math::SpatialPoint{
          (x0[0] + x1[0] + x2[0]) / Real(3),
          (x0[1] + x1[1] + x2[1]) / Real(3) });
    };

    std::vector<Index> stack;
    auto label = [&](Index s, char sv)
    {
      side[s] = sv;
      stack.push_back(s);
      while (!stack.empty())
      {
        const Index c = stack.back();
        stack.pop_back();
        for (Index e : conn.getIncidence({ 2, 1 }, c))
        {
          const bool flip = isIface(e);
          for (Index c2 : conn.getIncidence({ 1, 2 }, e))
          {
            if (c2 == c || side[c2] != 0) continue;
            side[c2] = flip ? static_cast<char>(3 - side[c]) : side[c];
            stack.push_back(c2);
          }
        }
      }
    };

    Index seed = 0;
    Real best = -std::numeric_limits<Real>::infinity();
    for (Index c = 0; c < nc; ++c)
    {
      const Real p = centroidPhi(c);
      if (p > best) { best = p; seed = c; }
    }
    label(seed, best < 0 ? char(1) : char(2));

    Index fallback = 0;
    for (Index c = 0; c < nc; ++c)
      if (side[c] == 0)
      { ++fallback; label(c, centroidPhi(c) < 0 ? char(1) : char(2)); }

    for (Index c = 0; c < nc; ++c)
      mesh.setAttribute({ 2, c }, side[c] == 1 ? Negative : Positive);
    return fallback;
  }

  void demoteInterface(LocalMesh& mesh)
  {
    auto& conn = mesh.getConnectivity();
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto attr = mesh.getAttribute(1, e);
      if (attr && *attr == Interface)
        mesh.setAttribute({ 1, e }, Attribute(0));
    }
  }

  // Linear-shadow diagnostic: corner-based qmin + inversion count.
  std::pair<Real, Index> linearQuality(const LocalMesh& mesh)
  {
    const auto& conn = mesh.getConnectivity();
    Real qmin = std::numeric_limits<Real>::infinity();
    Index inverted = 0;
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      const auto& cell = conn.getPolytope(2, c);
      const auto x0 = mesh.getVertexCoordinates(cell(0));
      const auto x1 = mesh.getVertexCoordinates(cell(1));
      const auto x2 = mesh.getVertexCoordinates(cell(2));
      const Real orient =
          (x1[0] - x0[0]) * (x2[1] - x0[1])
        - (x1[1] - x0[1]) * (x2[0] - x0[0]);
      if (orient <= Real(0)) ++inverted;
      const Real area = Real(0.5) * std::abs(orient);
      const Real den = (x1 - x0).squaredNorm()
                     + (x2 - x1).squaredNorm()
                     + (x0 - x2).squaredNorm();
      if (den > 0)
        qmin = std::min(qmin, Real(4) * std::sqrt(Real(3)) * area / den);
    }
    if (!std::isfinite(qmin)) qmin = 0;
    return { qmin, inverted };
  }
}

// =====================================================================
//  Pipeline: cut -> coarsen -> IsoparametricGeometry build -> Newton ->
//  syncLinearBackbone -> classify -> curved XDMF -> diagnostics.
// =====================================================================
int main(int, char**)
{
  constexpr size_t resolution = 20;
  constexpr Index  steps      = 8;
  constexpr Real   h          = Real(1) / static_cast<Real>(resolution - 1);

  LocalMesh background =
    LocalMesh::UniformGrid(Polytope::Type::Triangle, { resolution, resolution });
  background.scale(h);

  IO::XDMF xdmf("LevelSetMovingCircle");

  Real avgLinQmin = 0, avgLinFit = 0;
  Real avgCurvedQmin = 0, avgCurvedFitRms = 0, avgCurvedFitMax = 0;
  Real avgCurvedMinDet = 0;
  Real avgTargetMaxMu = 0, avgTargetMinDetT = 0;
  Index sumInverted = 0, sumOverlap = 0, sumCurvedInvalid = 0;
  Index sumTargetInvalid = 0, diagSteps = 0;

  for (Index step = 0; step < steps; ++step)
  {
    const Real t = static_cast<Real>(step)
      / static_cast<Real>(std::max<Index>(1, steps - 1));
    auto phi = [t](const Math::SpatialPoint& x) { return phiAt(x, t); };
    auto gradPhi = [t](const Math::SpatialPoint& x) { return gradPhiAt(x, t); };

    background.getConnectivity().compute(2, 1);
    background.getConnectivity().compute(1, 0);

    // === 1. Cut. =====================================================
    P1<Real, LocalMesh> phiSpace(background);
    GridFunction phiH(phiSpace);
    for (Index v = 0; v < background.getVertexCount(); ++v)
      phiH[v] = phi(background.getVertexCoordinates(v));

    auto cut = LevelSetDiscretizerTriangles(phiH)
      .setSignTolerance(1e-12)
      .setInterfaceAttribute(Interface)
      .setNegativeCellAttribute(Negative)
      .setPositiveCellAttribute(Positive)
      .setCrossingSnapTolerance(0.03)
      .setMinCutQuality(0.2)
      .discretize();
    LocalMesh mesh = std::move(cut.mesh);
    annotateBoxBoundary(mesh);

    // The analytic phi=0 projector is shared between stages 2b and 3.
    auto projectToInterface = [&](const Math::SpatialPoint& x)
    {
      Math::SpatialPoint p = x;
      for (int i = 0; i < 2; ++i)
      {
        const Real f = phi(p);
        const auto g = gradPhi(p);
        const Real gg = g[0] * g[0] + g[1] * g[1];
        if (gg < Real(1e-30)) break;
        p = Math::SpatialPoint{
          p[0] - f * g[0] / gg, p[1] - f * g[1] / gg };
      }
      return p;
    };

    // === 2a. Topology only: split / collapse / swap. =================
    //
    // setProtectedVertices stops the topology operators from removing or
    // relabelling any interface vertex. No vertex relocation is requested
    // here -- that's stage 2b.
    {
      auto params = TriangleMeshOptimizerParameters::levelSetCarryForward(h);
      params.featureSmoothingPasses = 0;
      params.smoothingPasses        = 0;
      TriangleMeshOptimizer()
        .setParameters(params)
        .setProtectedVertices(interfaceVertexMask(mesh))
        .optimize(mesh);
    }
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);

    // === 2b. Interface tangential redistribution. ====================
    //
    // Laplacian-with-projection along the interface chain; disperses
    // sliver corners so the geometric solve starts from a clean topology.
    InterfaceTangentialRelaxation()
      .setProjector(projectToInterface)
      .setFeatureAttribute(Interface)
      .setPasses(3)
      .setRelaxation(0.3)
      .relax(mesh);

    mesh.getConnectivity().compute(1, 2);
    const auto pre = linearQuality(mesh);

    // === 3. Upgrade to P2 -> solve TMOP -> move P2 mesh by u. ========
    //
    // Three free functions, FES-INDEPENDENT in their data flow:
    //   upgradeTransformations<RealH1Element<2>>    -- linear -> affine P2
    //   solveIsoparametricTMOP                      -- Newton + line search
    //                                                  + interface constraint
    //                                                  + curved-det gate
    //                                                  + moveMesh(mesh, u) at the end
    //   syncLinearBackbone                          -- P2 corners -> linear vertex coords
    //
    // The TMOP terms are composed with Rodin's existing variational form
    // language (operator+). The driver does not template on N terms.
    Variational::RealH1Element<2> fe(Polytope::Type::Triangle);
    upgradeTransformations(mesh, fe, Interface);

    VectorH1<2, LocalMesh> space(
        std::integral_constant<size_t, 2>{}, mesh, 2);
    GridFunction u(space);
    TrialFunction du(space);
    TestFunction  v (space);

    // Curved TMOP recipe:
    //  * metric  = shape+size blend barrier (gamma=0.5), NOT pure shape
    //              (pure shape is scale-rank-deficient -> slivers survive);
    //  * target  = a P2 fitted target whose interface midside nodes are
    //              projected onto analytic phi=0, with a small ideal-element
    //              bias. This keeps the target curved and fit-compatible
    //              while still pushing away from cut sliver geometry;
    //  * fit     = smooth level-set penalty, c_sigma-normalized so the
    //              weight is resolution-portable, ramped by CONTINUATION;
    //  * no deviation term (it fights both fit and quality).
    //
    // c_sigma = total interface length so the penalty is a mean-square fit.
    Real interfaceLen = 0;
    {
      const auto& conn = mesh.getConnectivity();
      for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
      {
        const auto attr = mesh.getAttribute(1, e);
        if (!attr || *attr != Interface) continue;
        const auto& edge = conn.getPolytope(1, e);
        interfaceLen +=
          (mesh.getVertexCoordinates(edge(1))
         - mesh.getVertexCoordinates(edge(0))).norm();
      }
    }
    const Real cSigma = std::max(interfaceLen, Real(1e-12));

    ShapeSizeBlendMetric metric(Real(0.5));
    ProjectedQualityTargetJacobian target(
        mesh, Interface, projectToInterface, Real(0.05));
    QualityTerm quality(metric, target, Real(1));
    AnalyticLevelSetFitTerm fit(
        phi, gradPhi, Optional<Attribute>(Interface), Real(1));
    fit.setNormalization(cSigma);
    AnalyticLevelSetFitTerm bfit(
        boxValue, boxGradient, Optional<Attribute>(Boundary), Real(1));
    bfit.setNormalization(cSigma);

    auto makeResidual = [&] {
      return quality.residual(u, v) + fit.residual(u, v)
           + bfit.residual(u, v);
    };
    auto makeTangent = [&] {
      return quality.tangent(u, du, v) + fit.tangent(u, du, v)
           + bfit.tangent(u, du, v);
    };
    auto energy = [&] {
      return quality.energy(u) + fit.energy(u) + bfit.energy(u);
    };

    // Weight continuation: start loose (Newton recovers quality on a still-
    // loose interface), then ramp the fit weight to tighten the fit on the
    // already-well-shaped mesh -- gets BOTH instead of trading one for the
    // other on a fixed weight. Each pass re-solves u=0 on the moved mesh.
    IsoparametricTMOPSolverReport report;
    for (const Real w : { Real(1), Real(2), Real(4) })
    {
      fit.setWeight(w);
      bfit.setWeight(w);
      u.getData().setZero();
      report = solveIsoparametricTMOP(
          mesh, fe, u, du, v,
          makeResidual, makeTangent, energy,
          Interface);
    }

    const auto curved = curvedMetrics(mesh, phi, Interface);
    const auto targetStats = targetQualityMetrics(mesh, metric, target);
    syncLinearBackbone(mesh, fe);
    demoteTransformations(mesh);

    // === 4. Classify, export, carry forward. =========================
    const Index ifaceGaps = classifySides(mesh, phi);

    P1<Real, LocalMesh> outSpace(mesh);
    GridFunction outPhi(outSpace);
    for (Index v = 0; v < mesh.getVertexCount(); ++v)
      outPhi[v] = phi(mesh.getVertexCoordinates(v));
    auto grid = xdmf.grid("fitted");
    grid.setMesh(mesh, IO::XDMF::MeshPolicy::Transient);
    grid.clear();
    grid.add("phi", outPhi, IO::XDMF::Center::Node);
    xdmf.write(t).flush();

    // === 5. Diagnostics: linear shadow vs curved truth. ==============
    const auto [lqmin, linv] = linearQuality(mesh);
    const auto featureMask = interfaceVertexMask(mesh);
    Real linFit = 0;
    for (Index v = 0; v < mesh.getVertexCount(); ++v)
      if (v < static_cast<Index>(featureMask.size()) && featureMask[v])
        linFit = std::max(linFit,
            std::abs(phi(mesh.getVertexCoordinates(v))));

    avgLinQmin       += lqmin;
    avgLinFit        += linFit;
    avgCurvedQmin    += curved.qmin;
    avgCurvedFitRms  += curved.fitRms;
    avgCurvedFitMax  += curved.fitMax;
    avgCurvedMinDet  += curved.minDet;
    avgTargetMaxMu   += targetStats.maxMetric;
    avgTargetMinDetT += targetStats.minDetT;
    sumInverted      += linv;
    sumCurvedInvalid += curved.invalidJacobianSamples;
    sumOverlap       += curved.overlapSamples;
    sumTargetInvalid += targetStats.invalidSamples;
    ++diagSteps;

    std::cout << "step " << step << " t=" << t
              << " cells=" << mesh.getCellCount()
              << " LIN[qmin=" << lqmin
              << " inv="      << linv
              << " fit="      << linFit << "]"
              << " CURV[qmin=" << curved.qmin
              << " fit_rms="   << curved.fitRms
              << " fit_max="   << curved.fitMax
              << " minDet="    << curved.minDet
              << " invalid="   << curved.invalidJacobianSamples
              << " overlap="   << curved.overlapSamples << "]"
              << " TMOP_TARGET[max_mu=" << targetStats.maxMetric
              << " minDetT=" << targetStats.minDetT
              << " invalid=" << targetStats.invalidSamples << "]"
              << " it="     << report.iterations
              << " r0="     << report.initialResidual
              << " r1="     << report.finalResidual
              << " gaps="   << ifaceGaps
              << " preq="   << pre.first
              << std::endl;

    demoteInterface(mesh);
    mesh.flush();
    background = std::move(mesh);
  }

  xdmf.close();

  if (diagSteps > 0)
  {
    const Real inv = Real(1) / static_cast<Real>(diagSteps);
    std::cout << "AVERAGE over " << diagSteps << " steps:"
              << " LIN[qmin="  << avgLinQmin * inv
              << " inv="       << static_cast<Real>(sumInverted) * inv
              << " fit="       << avgLinFit * inv << "]"
              << " CURV[qmin=" << avgCurvedQmin * inv
              << " fit_rms="   << avgCurvedFitRms * inv
              << " fit_max="   << avgCurvedFitMax * inv
              << " minDet="    << avgCurvedMinDet * inv
              << " invalid="   << static_cast<Real>(sumCurvedInvalid) * inv
              << " overlap="   << static_cast<Real>(sumOverlap) * inv
              << "]"
              << " TMOP_TARGET[max_mu=" << avgTargetMaxMu * inv
              << " minDetT=" << avgTargetMinDetT * inv
              << " invalid=" << static_cast<Real>(sumTargetInvalid) * inv
              << "]" << std::endl;
  }
  return 0;
}
