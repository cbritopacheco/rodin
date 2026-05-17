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

  void reclassifyCells(LocalMesh& mesh, Real t)
  {
    const auto& conn = mesh.getConnectivity();
    for (Index c = 0; c < static_cast<Index>(mesh.getCellCount()); ++c)
    {
      const auto& cell = conn.getPolytope(2, c);
      const auto x0 = mesh.getVertexCoordinates(cell(0));
      const auto x1 = mesh.getVertexCoordinates(cell(1));
      const auto x2 = mesh.getVertexCoordinates(cell(2));
      const Math::SpatialPoint centroid{
        (x0[0] + x1[0] + x2[0]) / Real(3),
        (x0[1] + x1[1] + x2[1]) / Real(3)
      };
      mesh.setAttribute({ 2, c },
          phiAt(centroid, t) < 0 ? Negative : Positive);
    }
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
  constexpr size_t resolution = 48;
  constexpr Index steps = 50;
  constexpr Real h = Real(1) / static_cast<Real>(resolution - 1);

  LocalMesh background =
    LocalMesh::UniformGrid(Polytope::Type::Triangle, { resolution, resolution });
  background.scale(h);

  IO::XDMF xdmf("LevelSetMovingWavyCircle");

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

    P1 space(cut.mesh, cut.mesh.getSpaceDimension());
    GridFunction u(space);
    u.getData().setZero();
    TrialFunction du(space);
    TestFunction v(space);

    auto phiValue = [t](const Math::SpatialPoint& x)
    {
      return phiAt(x, t);
    };
    auto phiGradient = [t](const Math::SpatialPoint& x)
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
    newton.setMaxIterations(5)
          .setDampingFactor(1)
          .setAbsoluteTolerance(1e-10)
          .setRelativeTolerance(1e-8);
    newton.solve(u);

    LocalMesh optimized = cut.mesh;
    applyDisplacement(optimized, u);

    const auto protectedVertices = interfaceVertexMask(optimized);
    const auto parameters = TriangleMeshOptimizerParameters::levelSetCarryForward(h);
    TriangleMeshOptimizer()
      .setParameters(parameters)
      .setFeatureProjection([t](const Math::SpatialPoint& x)
      {
        return projectToInterface(x, t);
      })
      .setProtectedVertices(protectedVertices)
      .optimize(optimized);

    optimized.getConnectivity().compute(2, 1);
    optimized.getConnectivity().compute(1, 0);
    reclassifyCells(optimized, t);

    P1<Real, LocalMesh> outputSpace(optimized);
    GridFunction outputPhi(outputSpace);
    for (Index i = 0; i < optimized.getVertexCount(); ++i)
      outputPhi[i] = phiAt(optimized.getVertexCoordinates(i), t);

    auto grid = xdmf.grid("fitted");
    grid.setMesh(optimized, IO::XDMF::MeshPolicy::Transient);
    grid.clear();
    grid.add("phi", outputPhi, IO::XDMF::Center::Node);
    xdmf.write(t).flush();

    const auto [qmin, inverted] = meshQuality(optimized);
    std::cout << "step " << step
              << " t=" << t
              << " cells=" << optimized.getCellCount()
              << " qmin=" << qmin
              << " inverted=" << inverted
              << std::endl;

    demoteInterface(optimized);
    background = std::move(optimized);
  }

  xdmf.close();
  return 0;
}
