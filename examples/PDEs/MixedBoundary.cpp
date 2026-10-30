/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [mixed-boundary-source]
#include <cmath>
#include <cstdlib>

#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Solver.h>
#include <Rodin/IO/XDMF.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main()
{
  //! [mixed-boundary-mesh]
  constexpr size_t nodes = 17;
  constexpr Attribute fixed = 2;
  constexpr Attribute loaded = 3;
  constexpr Attribute insulated = 4;
  constexpr Real tolerance = 1e-12; // Recognize coordinates on this unit grid.
  Mesh mesh = Mesh().UniformGrid(Polytope::Type::Triangle, {nodes, nodes});
  mesh.scale(1.0 / (nodes - 1));
  mesh.getConnectivity().compute(1, 2);
  for (auto edge = mesh.getFace(); !edge.end(); ++edge)
  {
    if (!edge->isBoundary())
      continue;
    bool onLeft = true;
    bool onRight = true;
    for (auto vertex = edge->getVertex(); !vertex.end(); ++vertex)
    {
      const Real x = vertex->getCoordinates().x();
      onLeft = onLeft && std::abs(x) < tolerance;
      onRight = onRight && std::abs(x - 1) < tolerance;
    }
    mesh.setAttribute({1, edge->getIndex()},
      onLeft      ? fixed
        : onRight ? loaded
                  : insulated);
  }
  //! [mixed-boundary-mesh]

  //! [mixed-boundary-form]
  P1 vh(mesh);
  TrialFunction u(vh);
  TestFunction v(vh);
  RealFunction flux = 1.0;
  Problem diffusion(u, v);
  diffusion = Integral(Grad(u), Grad(v)) - BoundaryIntegral(flux, v).over(loaded) +
    DirichletBC(u, Zero()).on(fixed);
  //! [mixed-boundary-form]

  //! [mixed-boundary-solve]
  Solver::CG solver(diffusion);
  solver.setTolerance(1e-10).solve();
  if (!solver.success())
    return EXIT_FAILURE;
  IO::XDMF xdmf("MixedBoundary");
  xdmf.grid().setMesh(mesh).add("u", u.getSolution());
  xdmf.write().flush();
  //! [mixed-boundary-solve]
  return 0;
}
//! [mixed-boundary-source]
