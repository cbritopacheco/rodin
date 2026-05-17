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

#include <Rodin/Adaptation.h>
#include <Rodin/Assembly/Default.h>
#include <Rodin/Geometry.h>
#include <Rodin/Geometry/LevelSetDiscretizerTriangles.h>
#include <Rodin/Geometry/TriangleMeshOptimizer.h>
#include <Rodin/IO.h>
#include <Rodin/Solver/NewtonSolver.h>
#include <Rodin/Solver/SparseLU.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Adaptation::TargetMatrixOptimization;
using namespace Rodin::Geometry;
using namespace Rodin::Solver;
using namespace Rodin::Variational;

// Moving level-set pipeline, one minimal loop per step:
//   1. project the analytic level set onto the current background,
//   2. robustly cut it (snap + min-cut-quality, no slivers),
//   3. recover fit + quality with one TargetMatrixOptimization Newton solve,
//   4. coarsen/optimize the mesh (native MMG-inspired hmin/hmax operators),
//   5. carry the optimized cut mesh forward as the next background.
int main(int, char**)
{
  constexpr Real Pi = 3.14159265358979323846;
  constexpr Attribute Interface = 30, Boundary = 40, Negative = 1, Positive = 2;
  constexpr size_t resolution = 24;
  constexpr Index steps = 20;
  constexpr Real radius = 0.18;

  auto center = [&](Real t)
  {
    return Math::SpatialPoint{ 0.5 + 0.14 * std::cos(2 * Pi * t),
                               0.5 + 0.11 * std::sin(2 * Pi * t) };
  };
  auto phiAt = [&](const Math::SpatialPoint& x, Real t)
  {
    const auto c = center(t);
    return std::hypot(x[0] - c[0], x[1] - c[1]) - radius;
  };
  auto gradPhiAt = [&](const Math::SpatialPoint& x, Real t)
  {
    const auto c = center(t);
    const Real dx = x[0] - c[0];
    const Real dy = x[1] - c[1];
    const Real r = std::hypot(dx, dy);
    if (r <= Real(1e-14))
      return Math::SpatialPoint{0, 0};
    return Math::SpatialPoint{dx / r, dy / r};
  };
  auto boxBoundaryValue = [](const Math::SpatialPoint& x)
  {
    const std::array<Real, 4> d = {{
      std::abs(x[0]), std::abs(x[0] - Real(1)),
      std::abs(x[1]), std::abs(x[1] - Real(1))
    }};
    const auto side =
      static_cast<size_t>(std::min_element(d.begin(), d.end()) - d.begin());
    if (side == 0) return x[0];
    if (side == 1) return x[0] - Real(1);
    if (side == 2) return x[1];
    return x[1] - Real(1);
  };
  auto boxBoundaryGradient = [](const Math::SpatialPoint& x)
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
  };

  // Worst normalized triangle quality and #inverted cells of a mesh.
  auto meshStats = [](const LocalMesh& m)
  {
    const auto& conn = m.getConnectivity();
    Real qmin = 1e30; Index inverted = 0;
    for (Index c = 0; c < static_cast<Index>(m.getCellCount()); ++c)
    {
      const auto& t = conn.getPolytope(2, c);
      const auto a = m.getVertexCoordinates(t(0));
      const auto b = m.getVertexCoordinates(t(1));
      const auto d = m.getVertexCoordinates(t(2));
      const Real ar = 0.5 * ((b[0]-a[0])*(d[1]-a[1]) - (b[1]-a[1])*(d[0]-a[0]));
      if (ar <= 0) ++inverted;
      const Real den = (b-a).squaredNorm() + (d-b).squaredNorm()
                     + (a-d).squaredNorm();
      if (den > 0)
        qmin = std::min(qmin, 4 * std::sqrt(3.0) * std::abs(ar) / den);
    }
    return std::pair<Real, Index>{ qmin, inverted };
  };

  // Geometric benchmark on mesh geometry (NOT the resampled analytic field):
  //  - interface fit error: |dist(x, circle)| over InterfaceAttribute-edge
  //    vertices, using the exact circle signed distance;
  //  - interior correspondence: area of the Negative region vs the true
  //    disk area pi r^2.
  // domainError: area where the MESHED domain (cell side attribute) disagrees
  // with the exact disk (sign of the analytic circle SDF at the cell
  // centroid). This is the integrated "SDF of the meshed domain vs analytic
  // SDF" check and directly exposes attribute checkerboard / wrong
  // inside-outside, independent of interface-vertex error.
  // coverage = Sigma|cell area|; signedArea = Sigma signed cell area. The
  // domain is the unit square, so a fold-free, single-cover mesh has both
  // ~= 1. coverage > 1 means the mesh laps over itself (overlap/fold);
  // signedArea < coverage means some cells are oppositely oriented.
  struct Bench { Real fitMax = 0, fitRms = 0, interiorArea = 0,
                      domainError = 0, coverage = 0, signedArea = 0; };
  auto benchmark = [&](const LocalMesh& m, Real time) -> Bench
  {
    const auto c = center(time);
    const auto& conn = m.getConnectivity();
    std::vector<char> onIface(m.getVertexCount(), 0);
    for (Index e = 0; e < static_cast<Index>(conn.getCount(1)); ++e)
    {
      const auto a = m.getAttribute(1, e);
      if (!a || *a != Interface) continue;
      const auto& ed = conn.getPolytope(1, e);
      onIface[ed(0)] = 1; onIface[ed(1)] = 1;
    }
    Bench b; Real sq = 0; Index n = 0;
    for (Index v = 0; v < static_cast<Index>(m.getVertexCount()); ++v)
    {
      if (!onIface[v]) continue;
      const auto x = m.getVertexCoordinates(v);
      const Real d = std::abs(std::hypot(x[0]-c[0], x[1]-c[1]) - radius);
      b.fitMax = std::max(b.fitMax, d); sq += d*d; ++n;
    }
    if (n) b.fitRms = std::sqrt(sq / n);
    for (Index k = 0; k < static_cast<Index>(m.getCellCount()); ++k)
    {
      const auto a = m.getCellAttribute(k);
      const auto& t = conn.getPolytope(2, k);
      const auto p0 = m.getVertexCoordinates(t(0));
      const auto p1 = m.getVertexCoordinates(t(1));
      const auto p2 = m.getVertexCoordinates(t(2));
      const Real signed2 =
          (p1[0]-p0[0])*(p2[1]-p0[1]) - (p1[1]-p0[1])*(p2[0]-p0[0]);
      const Real area = 0.5 * std::abs(signed2);
      b.coverage += area;
      b.signedArea += 0.5 * signed2;
      const bool meshInside = a && *a == Negative;
      if (meshInside)
        b.interiorArea += area;
      const Math::SpatialPoint cen{ (p0[0]+p1[0]+p2[0])/3,
                                    (p0[1]+p1[1]+p2[1])/3 };
      const bool trueInside =
        std::hypot(cen[0]-c[0], cen[1]-c[1]) - radius < 0;
      if (meshInside != trueInside)
        b.domainError += area;     // misclassified region (checkerboard)
    }
    return b;
  };
  const Real diskArea = Pi * radius * radius;

  LocalMesh background =
    LocalMesh::UniformGrid(Polytope::Type::Triangle, { resolution, resolution });
  background.scale(1.0 / (resolution - 1));
  Real qminStart = 0;
  Bench cutStart, tmopStart;

  IO::XDMF xdmf("LevelSetMovingCircle");

  for (Index step = 0; step < steps; ++step)
  {
    const Real t = static_cast<Real>(step) / (steps - 1);

    background.getConnectivity().compute(2, 1);
    background.getConnectivity().compute(1, 0);

    P1<Real, LocalMesh> phiSpace(background);
    GridFunction phi(phiSpace);
    for (Index i = 0; i < background.getVertexCount(); ++i)
      phi[i] = phiAt(background.getVertexCoordinates(i), t);

    // Robust cut: snap near-vertex crossings and keep a cell whole rather
    // than emit a sub-0.2-quality sliver. The interface it cannot fit is
    // recovered by the optimization below.
    auto cut = LevelSetDiscretizerTriangles(phi)
      .setSignTolerance(1e-12)
      .setInterfaceAttribute(Interface)
      .setNegativeCellAttribute(Negative)
      .setPositiveCellAttribute(Positive)
      .setCrossingSnapTolerance(0.05)
      .setMinCutQuality(0.2)
      .discretize();

    {
      const auto& cc = cut.mesh.getConnectivity();
      for (Index e = 0; e < static_cast<Index>(cc.getCount(1)); ++e)
      {
        const auto& edge = cc.getPolytope(1, e);
        const auto x0 = cut.mesh.getVertexCoordinates(edge(0));
        const auto x1 = cut.mesh.getVertexCoordinates(edge(1));
        const Real eps = 1e-12;
        const bool x0Side =
          std::abs(x0[0]) <= eps || std::abs(x0[0] - Real(1)) <= eps
          || std::abs(x0[1]) <= eps || std::abs(x0[1] - Real(1)) <= eps;
        const bool x1Side =
          std::abs(x1[0]) <= eps || std::abs(x1[0] - Real(1)) <= eps
          || std::abs(x1[1]) <= eps || std::abs(x1[1] - Real(1)) <= eps;
        if (x0Side && x1Side)
          cut.mesh.setAttribute({ 1, e }, Boundary);
      }
    }

    // Recover fit + quality: strict target-matrix optimization, one Newton
    // solve of quality + deviation + interface-fit terms.
    P1 space(cut.mesh, cut.mesh.getSpaceDimension());
    GridFunction u(space);
    u.getData().setZero();
    TrialFunction du(space);
    TestFunction v(space);

    // Quality fit + analytic interface fit. The quality term is strict
    // target-matrix TMOP with an oriented equilateral same-area target; the
    // analytic fit term keeps the fitted interface on the current circle
    // instead of only preserving the source P1 cut segment.
    auto phiValue = [&, t](const Math::SpatialPoint& x)
    {
      return phiAt(x, t);
    };
    auto phiGradient = [&, t](const Math::SpatialPoint& x)
    {
      return gradPhiAt(x, t);
    };
    QualityTerm quality(
        SquaredDistanceMetric{},
        OrientedEquilateralSameAreaTargetJacobian(cut.mesh),
        0.05);
    DeviationTerm deviation(1.0);
    AnalyticLevelSetFitTerm fit(
        phiValue, phiGradient, Optional<Attribute>(Interface), 1.0);
    // Softly preserve the carried square domain. Without this term, repeated
    // free-boundary quality recovery contracts the background even when every
    // triangle remains positively oriented.
    AnalyticLevelSetFitTerm boundaryFit(
        boxBoundaryValue, boxBoundaryGradient, Optional<Attribute>(Boundary), 1.0);

    Variational::Problem tmop(du, v);
    tmop = quality.tangent(u, du, v)
         + deviation.tangent(du, v)
         + fit.tangent(u, du, v)
         + boundaryFit.tangent(u, du, v)
         + quality.residual(u, v)
         + deviation.residual(u, v)
         + fit.residual(u, v)
         + boundaryFit.residual(u, v);

    SparseLU solver{tmop};
    NewtonSolver newton(solver);
    newton.setMaxIterations(5).setDampingFactor(1.0)
          .setAbsoluteTolerance(1e-10).setRelativeTolerance(1e-8);
    newton.solve(u);

    LocalMesh optimized = cut.mesh;
    {
      const auto& d = u.getData();
      const Index nv = optimized.getVertexCount();
      for (Index vtx = 0; vtx < nv; ++vtx)
      {
        auto x = optimized.getVertexCoordinates(vtx);
        x[0] += d(vtx);
        x[1] += d(vtx + nv);
        optimized.setVertexCoordinates(vtx, x);
      }
    }

    // Per-stage benchmark: after cut, and after TMOP (interface attributes
    // still present here; the optimizer stage does not yet preserve them).
    const Bench bCut = benchmark(cut.mesh, t);
    const Bench bTmop = benchmark(optimized, t);
    if (step == 0) { cutStart = bCut; tmopStart = bTmop; }

    // Protect the fitted interface: freeze every vertex on an
    // InterfaceAttribute edge so the coarsener preserves the fit exactly.
    std::vector<char> frozen(optimized.getVertexCount(), 0);
    {
      const auto& cc = cut.mesh.getConnectivity();
      for (Index e = 0; e < static_cast<Index>(cc.getCount(1)); ++e)
      {
        const auto attr = cut.mesh.getAttribute(1, e);
        if (!attr || *attr != Interface) continue;
        const auto& edge = cc.getPolytope(1, e);
        frozen[edge(0)] = 1;
        frozen[edge(1)] = 1;
      }
    }

    // Optimize/coarsen stage. Keep the freshly fitted interface frozen, allow
    // conservative interior collapses to remove stale cut topology, disable
    // splitting for now, and keep quality-improving swaps enabled. Earlier
    // overlap came from insufficient swap/collapse validity guards; with those
    // guards in place, conservative collapse is the mechanism that prevents
    // carried-forward meshes from accumulating old interface triangles.
    const auto optReport = TriangleMeshOptimizer()
      .setHMin(0.35 / (resolution - 1))
      .setHMax(std::numeric_limits<Real>::infinity())
      .setMinQuality(0.1)
      .setMaxIterations(4)
      .setProtectedVertices(std::move(frozen))
      .optimize(optimized);

    optimized.getConnectivity().compute(2, 1);
    optimized.getConnectivity().compute(1, 0);

    // The optimizer shuffles triangles between regions while keeping per-slot
    // attributes, which scrambles the Negative/Positive side (checkerboard;
    // benchmark localized domainErr to exactly this stage). Since phi is
    // analytic, re-derive the side attribute exactly from the sign of phi at
    // each cell centroid. The interface edges (preserved by the optimizer)
    // are unaffected.
    {
      const auto& oc = optimized.getConnectivity();
      for (Index k = 0; k < static_cast<Index>(optimized.getCellCount()); ++k)
      {
        const auto& tt = oc.getPolytope(2, k);
        const auto q0 = optimized.getVertexCoordinates(tt(0));
        const auto q1 = optimized.getVertexCoordinates(tt(1));
        const auto q2 = optimized.getVertexCoordinates(tt(2));
        const Math::SpatialPoint cen{ (q0[0]+q1[0]+q2[0])/3,
                                      (q0[1]+q1[1]+q2[1])/3 };
        optimized.setAttribute({ 2, k },
          phiAt(cen, t) < 0 ? Negative : Positive);
      }
    }

    // After-optimize stage: interface attributes are preserved through the
    // optimizer now, so the fit is measurable here too.
    const Bench bOpt = benchmark(optimized, t);

    P1<Real, LocalMesh> outSpace(optimized);
    GridFunction outPhi(outSpace);
    for (Index i = 0; i < optimized.getVertexCount(); ++i)
      outPhi[i] = phiAt(optimized.getVertexCoordinates(i), t);

    auto grid = xdmf.grid("fitted");
    grid.setMesh(optimized, IO::XDMF::MeshPolicy::Transient);
    grid.clear();
    grid.add("phi", outPhi, IO::XDMF::Center::Node);

    // Visualize the TMOP displacement on the pre-coarsen (cut) mesh, where u
    // is defined.
    auto dgrid = xdmf.grid("displacement");
    dgrid.setMesh(cut.mesh, IO::XDMF::MeshPolicy::Transient);
    dgrid.clear();
    dgrid.add("u", u, IO::XDMF::Center::Node);

    xdmf.write(t).flush();

    const auto [qmin, inverted] = meshStats(optimized);
    if (step == 0) qminStart = qmin;

    std::cout << "step " << step << " t=" << t
              << " cut_cells=" << cut.mesh.getCellCount()
              << " uncut=" << cut.report.uncutCellCount
              << " final_cells=" << optimized.getCellCount()
              << " opt(s/c/w)=" << optReport.splits << "/"
              << optReport.collapses << "/" << optReport.swaps
              << " qmin=" << qmin << " inverted=" << inverted
              << " | fitRms cut/tmop/opt=" << bCut.fitRms << "/"
              << bTmop.fitRms << "/" << bOpt.fitRms
              << " areaErr cut/tmop/opt="
              << std::abs(bCut.interiorArea - diskArea) << "/"
              << std::abs(bTmop.interiorArea - diskArea) << "/"
              << std::abs(bOpt.interiorArea - diskArea)
              << " domainErr cut/tmop/opt=" << bCut.domainError << "/"
              << bTmop.domainError << "/" << bOpt.domainError
              << " coverage opt=" << bOpt.coverage
              << " signedArea opt=" << bOpt.signedArea
              << " (domain=1)" << std::endl;

    // Round-trip benchmark: center(t) is periodic (t=0 == t=1), so the
    // carried-forward mesh should come home with comparable fit and quality.
    if (step + 1 == steps)
      std::cout << "ROUND-TRIP"
                << " qmin start/end=" << qminStart << "/" << qmin
                << " fitRms cut start/end=" << cutStart.fitRms << "/"
                << bCut.fitRms
                << " fitRms tmop start/end=" << tmopStart.fitRms << "/"
                << bTmop.fitRms
                << " interiorArea start/end (true=" << diskArea << ")="
                << tmopStart.interiorArea << "/" << bTmop.interiorArea
                << " domainErr start/end=" << tmopStart.domainError << "/"
                << bTmop.domainError
                << " coverage start/end=" << tmopStart.coverage << "/"
                << bOpt.coverage
                << " inverted=" << inverted << std::endl;

    // Transient interface: the interface belongs only to THIS time step.
    // Demote it to ordinary mesh before carry-forward so it is not frozen
    // again, does not accumulate as a swept band, and the next step's
    // coarsener can dissolve it. The next cut re-establishes a fresh
    // interface from analytic phi for the new circle position.
    {
      auto& oc = optimized.getConnectivity();
      for (Index e = 0; e < static_cast<Index>(oc.getCount(1)); ++e)
      {
        const auto a = optimized.getAttribute(1, e);
        if (a && *a == Interface)
          optimized.setAttribute({ 1, e }, Attribute(0));
      }
    }

    // Carry forward: the optimizer rebuilds a clean triangle-only mesh, so it
    // is re-cuttable; with the interface demoted, the coarsener bounds growth
    // across the orbit instead of accumulating frozen swept bands.
    background = std::move(optimized);
  }

  xdmf.close();
  return 0;
}
