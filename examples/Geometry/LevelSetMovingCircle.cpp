/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Adaptation/TargetMatrixOptimization/Metrics.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

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

namespace
{
  constexpr Real Pi = 3.14159265358979323846;
  constexpr Attribute Interface = 30;
  constexpr Attribute Boundary = 40;
  constexpr Attribute Negative = 1;
  constexpr Attribute Positive = 2;

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

  Math::SpatialPoint projectToInterface(
      const Math::SpatialPoint& x,
      Real t)
  {
    Math::SpatialPoint p = x;
    for (int i = 0; i < 2; ++i)
    {
      const Real f = phiAt(p, t);
      const auto g = gradPhiAt(p, t);
      const Real gg = g[0] * g[0] + g[1] * g[1];
      if (gg < Real(1e-30))
        break;
      p = Math::SpatialPoint{ p[0] - f * g[0] / gg,
                              p[1] - f * g[1] / gg };
    }
    return p;
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

  std::vector<char> interfaceVertexMask(const LocalMesh& mesh)
  {
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
    return mask;
  }

  template <class GridFunctionType>
  void applyDisplacement(LocalMesh& mesh, const GridFunctionType& u)
  {
    const auto& data = u.getData();
    const Index nv = mesh.getVertexCount();
    for (Index v = 0; v < nv; ++v)
    {
      auto x = mesh.getVertexCoordinates(v);
      x[0] += data(v);
      x[1] += data(v + nv);
      mesh.setVertexCoordinates(v, x);
    }
  }

  std::vector<Math::SpatialPoint> vertexCoordinates(const LocalMesh& mesh)
  {
    std::vector<Math::SpatialPoint> coordinates;
    coordinates.reserve(static_cast<size_t>(mesh.getVertexCount()));
    for (Index v = 0; v < mesh.getVertexCount(); ++v)
      coordinates.push_back(mesh.getVertexCoordinates(v));
    return coordinates;
  }

  void restoreVertexCoordinates(
      LocalMesh& mesh,
      const std::vector<Math::SpatialPoint>& coordinates)
  {
    const Index count = std::min(
        mesh.getVertexCount(),
        static_cast<Index>(coordinates.size()));
    for (Index v = 0; v < count; ++v)
      mesh.setVertexCoordinates(v, coordinates[static_cast<size_t>(v)]);
    mesh.getConnectivity().compute(2, 1);
    mesh.getConnectivity().compute(1, 0);
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

  std::pair<Real, Index> meshQuality(const LocalMesh& mesh)
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
      if (orient <= 0)
        ++inverted;
      const Real area = Real(0.5) * std::abs(orient);
      const Real den = (x1 - x0).squaredNorm()
                     + (x2 - x1).squaredNorm()
                     + (x0 - x2).squaredNorm();
      if (den > 0)
        qmin = std::min(qmin, Real(4) * std::sqrt(Real(3)) * area / den);
    }
    if (!std::isfinite(qmin))
      qmin = 0;
    return {qmin, inverted};
  }
}

int main(int, char**)
{
  constexpr size_t resolution = 32;
  constexpr Index steps = 40;
  constexpr Real h = Real(1) / static_cast<Real>(resolution - 1);

  LocalMesh background =
    LocalMesh::UniformGrid(Polytope::Type::Triangle, { resolution, resolution });
  background.scale(h);

  IO::XDMF xdmf("LevelSetMovingCircle");

  // Diagnostics are accumulated and reported as an average over the whole
  // process, not just the final step.
  Real sumQmin = 0, sumFit = 0;
  Index sumInverted = 0, rejectedTMOP = 0, diagSteps = 0;

  for (Index step = 0; step < steps; ++step)
  {
    const Real t = static_cast<Real>(step)
      / static_cast<Real>(std::max<Index>(1, steps - 1));

    background.getConnectivity().compute(2, 1);
    background.getConnectivity().compute(1, 0);

    P1<Real, LocalMesh> phiSpace(background);
    GridFunction phi(phiSpace);
    for (Index i = 0; i < background.getVertexCount(); ++i)
      phi[i] = phiAt(background.getVertexCoordinates(i), t);

    auto cut = LevelSetDiscretizerTriangles(phi)
      .setSignTolerance(1e-12)
      .setInterfaceAttribute(Interface)
      .setNegativeCellAttribute(Negative)
      .setPositiveCellAttribute(Positive)
      .setCrossingSnapTolerance(0.05)
      .setMinCutQuality(0.2)
      .discretize();
    annotateBoundary(cut.mesh);

    // ===== PRODUCTION / PDE-ready order =====
    //   cut -> optimize/coarsen -> rebuild feature targets -> TMOP (last)
    // Topology-changing operations (coarsen/swap/collapse/feature-smooth)
    // run BEFORE TMOP. TMOP is the final mesh-quality recovery step before
    // the PDE solve / export; nothing changes topology after it. The
    // alternative cut -> TMOP -> optimize order is kept only as a
    // benchmark-matrix diagnostic (StageCase), never as the main pipeline.

    // 1) optimize / coarsen on the raw cut mesh (interface frozen, feature
    //    vertices smoothed tangentially on phi=0).
    LocalMesh optimized = cut.mesh;
    {
      const auto cutInterfaceMask = interfaceVertexMask(optimized);
      const auto parameters =
        TriangleMeshOptimizerParameters::levelSetCarryForward(h);
      TriangleMeshOptimizer()
        .setParameters(parameters)
        .setFeatureProjection([t](const Math::SpatialPoint& x)
          { return projectToInterface(x, t); })
        .setProtectedVertices(cutInterfaceMask)
        .optimize(optimized);
    }
    optimized.getConnectivity().compute(2, 1);
    optimized.getConnectivity().compute(1, 0);

    // 2) rebuild feature/target metadata on the CURRENT (coarsened) mesh.
    //    Cut provenance is stale after topology changes, so the protected
    //    feature mask and the per-cell target Jacobian are rebuilt here from
    //    the coarsened mesh, before TMOP. The fit is analytic in phi(t), so
    //    it is inherently current and needs no rebuild.
    const auto featureMask = interfaceVertexMask(optimized);
    IdealElementTargetJacobian target(optimized);

    // 3) TMOP LAST: final mesh-quality recovery before PDE/export.
    P1 space(optimized, optimized.getSpaceDimension());
    optimized.getConnectivity().compute(1, 2);
    GridFunction u(space);
    u.getData().setZero();
    TrialFunction du(space);
    TestFunction v(space);

    auto phiValue = [t](const Math::SpatialPoint& x) { return phiAt(x, t); };
    auto phiGradient =
      [t](const Math::SpatialPoint& x) { return gradPhiAt(x, t); };

    QualityTerm quality(ShapeSizeMetric{}, target, 0.03);
    DeviationTerm deviation(1.0);
    AnalyticLevelSetFitTerm fit(
        phiValue, phiGradient, Optional<Attribute>(Interface), 1.0);
    AnalyticLevelSetFitTerm boundaryFit(
        boxBoundaryValue, boxBoundaryGradient,
        Optional<Attribute>(Boundary), 1.0);

    Variational::Problem tmop(du, v);
    tmop = quality.tangent(u, du, v)
         + deviation.tangent(du, v)
         + fit.tangent(u, du, v)
         + boundaryFit.tangent(u, du, v)
         + quality.residual(u, v)
         + deviation.residual(u, v)
         + fit.residual(u, v)
         + boundaryFit.residual(u, v)
         // + DirichletBC(du, Zero(space.getMesh().getSpaceDimension()));
         ;

    SparseLU solver{tmop};
    NewtonSolver newton(solver);
    newton.setMaxIterations(5)
          .setDampingFactor(1)
          .setAbsoluteTolerance(1e-10)
          .setRelativeTolerance(1e-8);
    newton.solve(u);
    const auto beforeTMOP = vertexCoordinates(optimized);
    applyDisplacement(optimized, u);

    optimized.getConnectivity().compute(2, 1);
    optimized.getConnectivity().compute(1, 0);
    const auto [tmopQmin, tmopInverted] = meshQuality(optimized);
    const bool acceptedTMOP =
      tmopInverted == 0 && std::isfinite(tmopQmin) && tmopQmin > Real(0);
    if (!acceptedTMOP)
    {
      restoreVertexCoordinates(optimized, beforeTMOP);
      ++rejectedTMOP;
    }
    optimized.getConnectivity().compute(2, 1);
    optimized.getConnectivity().compute(1, 0);
    // Cell side attributes come from the cutter's exact P1 side topology and
    // are preserved through the optimizer. Do not overwrite them with a
    // centroid phi heuristic after TMOP: fitted/interface-adjacent cells can be
    // mis-tagged that way, especially under strong quality weights.

    // 4) export (a PDE solve would run on `optimized` here).
    P1<Real, LocalMesh> outputSpace(optimized);
    GridFunction outputPhi(outputSpace);
    for (Index i = 0; i < optimized.getVertexCount(); ++i)
      outputPhi[i] = phiAt(optimized.getVertexCoordinates(i), t);

    auto grid = xdmf.grid("fitted");
    grid.setMesh(optimized, IO::XDMF::MeshPolicy::Transient);
    grid.clear();
    grid.add("phi", outputPhi, IO::XDMF::Center::Node);
    xdmf.write(t).flush();

    // 5) diagnostics (accumulated for the whole-process average).
    const auto [qmin, inverted] = meshQuality(optimized);
    Real fitErr = 0;
    for (Index i = 0; i < optimized.getVertexCount(); ++i)
      if (i < featureMask.size() && featureMask[i])
        fitErr = std::max(fitErr,
            std::abs(phiAt(optimized.getVertexCoordinates(i), t)));
    sumQmin += qmin;
    sumInverted += inverted;
    sumFit += fitErr;
    ++diagSteps;
    std::cout << "step " << step
              << " t=" << t
              << " cells=" << optimized.getCellCount()
              << " qmin=" << qmin
              << " inverted=" << inverted
              << " fit=" << fitErr
              << " tmop_accepted=" << (acceptedTMOP ? 1 : 0)
              << std::endl;

    demoteInterface(optimized);
    optimized.flush();
    background = std::move(optimized);
  }

  xdmf.close();

  if (diagSteps > 0)
    std::cout << "AVERAGE over " << diagSteps << " steps:"
              << " qmin=" << sumQmin / static_cast<Real>(diagSteps)
              << " inverted="
              << static_cast<Real>(sumInverted) / static_cast<Real>(diagSteps)
              << " fit=" << sumFit / static_cast<Real>(diagSteps)
              << " rejected_tmop=" << rejectedTMOP
              << std::endl;
  return 0;
}
