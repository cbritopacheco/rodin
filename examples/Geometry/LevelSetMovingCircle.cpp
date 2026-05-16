/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
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
  constexpr Attribute Interface = 30, Negative = 1, Positive = 2;
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

  LocalMesh background =
    LocalMesh::UniformGrid(Polytope::Type::Triangle, { resolution, resolution });
  background.scale(1.0 / (resolution - 1));
  Real qminStart = 0;

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

    // Recover fit + quality: strict target-matrix optimization, one Newton
    // solve of quality + deviation + interface-fit terms.
    P1 space(cut.mesh, cut.mesh.getSpaceDimension());
    GridFunction u(space);
    u.getData().setZero();
    TrialFunction du(space);
    TestFunction v(space);

    QualityTerm quality(
        SquaredDistanceMetric{}, EquilateralSameAreaTargetJacobian(cut.mesh), 1.0);
    DeviationTerm deviation(1.0);
    LevelSetFitTerm fit(cut.interfaceGraph, cut.report, 1.0);

    Variational::Problem tmop(du, v);
    tmop = quality.tangent(u, du, v)
         + deviation.tangent(du, v)
         + fit.tangent(u, du, v)
         + quality.residual(u, v)
         + deviation.residual(u, v)
         + fit.residual(u, v);

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

    // Coarsen/optimize: native edge split/collapse/swap with size bounds,
    // interface frozen.
    TriangleMeshOptimizer()
      .setHMin(0.4 / (resolution - 1))
      .setHMax(2.5 / (resolution - 1))
      .setMinQuality(0.1)
      .setMaxIterations(4)
      .setProtectedVertices(std::move(frozen))
      .optimize(optimized);

    optimized.getConnectivity().compute(2, 1);
    optimized.getConnectivity().compute(1, 0);

    P1<Real, LocalMesh> outSpace(optimized);
    GridFunction outPhi(outSpace);
    for (Index i = 0; i < optimized.getVertexCount(); ++i)
      outPhi[i] = phiAt(optimized.getVertexCoordinates(i), t);

    auto grid = xdmf.grid("fitted");
    grid.setMesh(optimized, IO::XDMF::MeshPolicy::Transient);
    grid.clear();
    grid.add("phi", outPhi, IO::XDMF::Center::Node);
    xdmf.write(t).flush();

    const auto [qmin, inverted] = meshStats(optimized);
    if (step == 0) qminStart = qmin;

    std::cout << "step " << step << " t=" << t
              << " cut_cells=" << cut.mesh.getCellCount()
              << " uncut=" << cut.report.uncutCellCount
              << " final_cells=" << optimized.getCellCount()
              << " qmin=" << qmin
              << " inverted=" << inverted << std::endl;

    // Round-trip benchmark: the circle returns to its t=0 position at t=1
    // (center is periodic), so the carried-forward mesh should come home with
    // quality comparable to the start.
    if (step + 1 == steps)
      std::cout << "ROUND-TRIP qmin start=" << qminStart
                << " end=" << qmin
                << " ratio=" << (qminStart > 0 ? qmin / qminStart : 0.0)
                << " inverted=" << inverted << std::endl;

    // Carry forward: the optimizer rebuilds a clean triangle-only mesh, so it
    // is re-cuttable; the coarsener bounds its growth across the orbit.
    background = std::move(optimized);
  }

  xdmf.close();
  return 0;
}
