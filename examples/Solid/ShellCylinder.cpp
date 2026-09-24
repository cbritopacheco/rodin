/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file ShellCylinder.cpp
 * @brief 3D-shell model of a cylinder under periodic pressure, after
 * Chapelle, Mardare & Münch, J. Elasticity 76 (2004) 199--246, Section 5.
 *
 * The 3D-shell model is linear elasticity restricted to displacements that
 * are quadratic across the thickness. One layer of curved P2 wedges realises
 * exactly that space: the wedge basis is (P2 triangle) x (P2 segment) and
 * the segment runs through the thickness.
 *
 * The modified model of Section 4.4 truncates the volumetric term to its
 * linear part in the transverse coordinate. It is realised here by the
 * paper's own shortcut: the lambda term is integrated with two Gauss points
 * through the thickness (wedge quadrature order 3), the mu term is not.
 *
 * Geometry: cylinder of radius R and length 2R, load p0 cos(2 theta) acting
 * radially, one eighth computed (x in [0, R], theta in [0, pi/2]).
 * Symmetry planes carry a normal-displacement penalty. The end x = R is free
 * (non-inhibited pure bending) or clamped (inhibited).
 *
 * Output: the compliance of the full cylinder scaled as in the paper,
 * 12 p0^2 R^3 / (E eps^3) for the free case and 12 p0^2 R^3 / (E eps) for
 * the clamped case.
 *
 * Usage: RodinSolidShellCylinder N nu free|clamped standard|modified
 *                                [eps] [penalty]
 */
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <Rodin/Solver/SparseLU.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int argc, char** argv)
{
  if (argc < 5)
  {
    std::cerr << "Usage: " << argv[0]
      << " N nu free|clamped standard|modified [eps] [penalty]\n";
    return 1;
  }

  const size_t n = std::stoul(argv[1]);
  const Real nu = std::stod(argv[2]);
  const bool clamped = std::string(argv[3]) == "clamped";
  const bool modified = std::string(argv[4]) == "modified";
  const Real eps = argc > 5 ? std::stod(argv[5]) : 1e-2;
  const Real penalty = argc > 6 ? std::stod(argv[6]) : 1e8;

  const Real radius = 1, young = 1, p0 = 1;
  const Real thickness = eps * radius;
  const Real mu = young / (2 * (1 + nu));
  const Real lambda = young * nu / ((1 + nu) * (1 - 2 * nu));

  // Boundary attributes
  const Attribute midPlane = 1, end = 2, planeZ = 3, planeY = 4;

  // Index-space grid: i along the axis, j around, k across the thickness.
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, {n + 1, n + 1, 2});
  const size_t d = mesh.getDimension();
  mesh.getConnectivity().compute(d - 1, d);
  mesh.getConnectivity().compute(d, d - 1);
  mesh.getConnectivity().compute(d - 1, 0);

  // Midsurface chart times the thickness coordinate:
  // Phi(x, theta, s) = (x, (R + s) cos theta, (R + s) sin theta).
  const auto chart =
    [&](const Math::SpatialPoint& q)
    {
      const Real x = q(0) / n * radius;
      const Real theta = q(1) / n * M_PI / 2;
      const Real r = radius + (q(2) - 0.5) * thickness;
      Math::SpatialPoint p(3);
      p << x, r * std::cos(theta), r * std::sin(theta);
      return p;
    };

  // Label the boundary from the index-space centroid of each face.
  for (auto it = mesh.getBoundary(); it; ++it)
  {
    Math::SpatialPoint c = Math::SpatialPoint::Zero(3);
    for (const auto& v : it->getVertices())
      c += mesh.getVertexCoordinates(v);
    c /= it->getVertices().size();
    Optional<Attribute> attr;
    if (c(0) < 1e-12)
      attr = midPlane;
    else if (c(0) > n - 1e-12)
      attr = end;
    else if (c(1) < 1e-12)
      attr = planeZ;
    else if (c(1) > n - 1e-12)
      attr = planeY;
    mesh.setAttribute({d - 1, it->getIndex()}, attr);
  }

  // Isoparametric P2 geometry: sample the chart at the geometry nodes.
  const RealH1Element<2> geometryFE(Polytope::Type::Wedge);
  std::vector<PointCloud> nodes;
  nodes.reserve(mesh.getCellCount());
  for (auto it = mesh.getCell(); it; ++it)
  {
    PointCloud pc(3, geometryFE.getCount());
    for (size_t local = 0; local < geometryFE.getCount(); ++local)
    {
      Math::SpatialPoint q;
      it->getTransformation().transform(q, geometryFE.getNode(local));
      const auto p = chart(q);
      for (size_t c = 0; c < 3; ++c)
        pc(c, local) = p(c);
    }
    nodes.push_back(std::move(pc));
  }
  for (Index v = 0; v < mesh.getVertexCount(); ++v)
    mesh.setVertexCoordinates(v, chart(mesh.getVertexCoordinates(v)));
  for (Index i = 0; i < mesh.getCellCount(); ++i)
  {
    mesh.setPolytopeTransformation({d, i},
      new ParametricTransformation<RealH1Element<2>>(std::move(nodes[i]), geometryFE));
  }

  H1 vh(std::integral_constant<size_t, 2>{}, mesh, d);
  TrialFunction u(vh);
  TestFunction  v(vh);

  // Radial load p0 cos(2 theta), spread uniformly over the thickness.
  VectorFunction f{
    0,
    [&](const Point& p) { return p0 / thickness * std::cos(2 * std::atan2(p.z(), p.y())) * p.y() / std::hypot(p.y(), p.z()); },
    [&](const Point& p) { return p0 / thickness * std::cos(2 * std::atan2(p.z(), p.y())) * p.z() / std::hypot(p.y(), p.z()); } };

  const VectorFunction ex{1, 0, 0}, ey{0, 1, 0}, ez{0, 0, 1};

  auto volumetric = Integral(lambda * Div(u), Div(v));
  if (modified)
    volumetric.setOrder(3);

  Problem shell(u, v);
  shell = Integral(mu * Jacobian(u), Jacobian(v))
        + Integral(mu * Jacobian(u), Transpose(Jacobian(v)))
        + volumetric
        + BoundaryIntegral(penalty * Dot(ex, u), Dot(ex, v)).over(midPlane)
        + BoundaryIntegral(penalty * Dot(ez, u), Dot(ez, v)).over(planeZ)
        + BoundaryIntegral(penalty * Dot(ey, u), Dot(ey, v)).over(planeY)
        - Integral(f, v);
  if (clamped)
    shell += DirichletBC(u, VectorFunction{0, 0, 0}).on(end);

  Solver::SparseLU(shell).solve();

  LinearForm load(v);
  load = Integral(f, v);
  load.assemble();
  const Real compliance = 8 * load.getVector().dot(u.getSolution().getData());
  const Real scale = clamped
    ? 12 * p0 * p0 * std::pow(radius, 3) / (young * eps)
    : 12 * p0 * p0 * std::pow(radius, 3) / (young * std::pow(eps, 3));

  std::cout.precision(8);
  std::cout << "N=" << n << " nu=" << nu
    << " " << (clamped ? "clamped" : "free")
    << " " << (modified ? "modified" : "standard")
    << " eps=" << eps << " penalty=" << penalty
    << " dofs=" << vh.getSize()
    << " energy=" << compliance / scale << std::endl;

  return 0;
}
