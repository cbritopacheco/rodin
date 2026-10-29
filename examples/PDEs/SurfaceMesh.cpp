/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [surfacemesh-source]
#include <cmath>
#include <cstdlib>

#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <Rodin/IO/XDMF.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main()
{
  //! [surface-mesh]
  constexpr size_t nodes = 11;
  Mesh volume = Mesh().UniformGrid(Polytope::Type::Tetrahedron, {nodes, nodes, nodes});
  volume.scale(2.0 / (nodes - 1));
  volume.displace(VectorFunction{-1, -1, -1});
  volume.getConnectivity().compute(2, 3);
  Mesh mesh = volume.skin();
  for (Index vertex = 0; vertex < mesh.getVertexCount(); ++vertex)
  {
    Math::SpatialPoint point = mesh.getVertexCoordinates(vertex);
    point.normalize();
    mesh.setVertexCoordinates(vertex, point);
  }
  //! [surface-mesh]

  //! [surface-form]
  P1 vh(mesh);
  TrialFunction u(vh);
  TestFunction v(vh);
  RealFunction source = [](const Point& point) { return 7 * point.x() * point.y(); };
  Problem reactionDiffusion(u, v);
  reactionDiffusion = Integral(Grad(u), Grad(v)) + Integral(u, v) - Integral(source, v);
  //! [surface-form]

  //! [surface-solve]
  Solver::CG solver(reactionDiffusion);
  solver.setTolerance(1e-10).solve();
  if (!solver.success())
    return EXIT_FAILURE;
  IO::XDMF xdmf("SurfaceMesh");
  xdmf.grid().setMesh(mesh).add("u", u.getSolution());
  xdmf.write().flush();
  //! [surface-solve]
  return 0;
}
//! [surfacemesh-source]
