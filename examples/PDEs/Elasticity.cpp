/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [elasticity-source]
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
  //! [elasticity-mesh]
  constexpr Real spacing = 0.1;
  constexpr Real length = 4;
  constexpr Real boundaryTolerance = 1e-12; // Coordinate tolerance for this grid.
  constexpr Attribute fixed = 2;
  constexpr Attribute loaded = 3;
  constexpr Attribute free = 4;

  Mesh mesh = Mesh().UniformGrid(Polytope::Type::Triangle, {41, 11});
  mesh.scale(spacing);
  mesh.getConnectivity().compute(1, 2);

  for (auto face = mesh.getFace(); !face.end(); ++face)
  {
    if (!face->isBoundary())
      continue;
    bool onLeft = true;
    bool onRight = true;
    for (auto vertex = face->getVertex(); !vertex.end(); ++vertex)
    {
      const Real x = vertex->getCoordinates().x();
      onLeft = onLeft && std::abs(x) < boundaryTolerance;
      onRight = onRight && std::abs(x - length) < boundaryTolerance;
    }
    mesh.setAttribute({1, face->getIndex()}, onLeft ? fixed : onRight ? loaded : free);
  }
  //! [elasticity-mesh]

  //! [elasticity-form]
  P1 vh(mesh, mesh.getSpaceDimension());
  TrialFunction u(vh);
  TestFunction v(vh);

  constexpr Real lambda = 0.5769;
  constexpr Real mu = 0.3846;
  VectorFunction traction{0, -1e-4};

  Problem elasticity(u, v);
  elasticity = LinearElasticityIntegral(u, v)(lambda, mu) -
    BoundaryIntegral(traction, v).over(loaded) + DirichletBC(u, Zero(2)).on(fixed);
  //! [elasticity-form]

  //! [elasticity-solve]
  Solver::CG solver(elasticity);
  solver.setTolerance(1e-10).solve();
  if (!solver.success())
    return EXIT_FAILURE;
  //! [elasticity-solve]

  //! [elasticity-export]
  IO::XDMF xdmf("Elasticity");
  xdmf.grid().setMesh(mesh).add("displacement", u.getSolution());
  xdmf.write().flush();
  //! [elasticity-export]
  return 0;
}
//! [elasticity-source]
