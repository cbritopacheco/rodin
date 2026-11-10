/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [periodic-source]
#include <map>
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
  //! [periodic-mesh]
  constexpr size_t nodes = 33;
  const Real length = 2 * Math::Constants::pi();
  Mesh mesh = Mesh().UniformGrid(Polytope::Type::Triangle, {nodes, nodes});
  mesh.scale(length / (nodes - 1));
  mesh.getConnectivity().compute(1, 2);
  P1 vh(mesh);
  //! [periodic-mesh]

  //! [periodic-map]
  std::map<std::pair<size_t, size_t>, Index> vertexDOFs;
  for (auto vertex = mesh.getVertex(); !vertex.end(); ++vertex)
  {
    const auto& point = vertex->getCoordinates();
    const size_t i = std::llround(point.x() * (nodes - 1) / length);
    const size_t j = std::llround(point.y() * (nodes - 1) / length);
    vertexDOFs[{i, j}] = vh.getDOFs(0, vertex->getIndex())[0];
  }
  IndexMap<IndexSet> periodicDOFs;
  for (const auto& [position, dof] : vertexDOFs)
  {
    const auto [i, j] = position;
    if (i == nodes - 1 || j == nodes - 1)
      periodicDOFs[dof].insert(vertexDOFs.at({i % (nodes - 1), j % (nodes - 1)}));
  }
  //! [periodic-map]

  //! [periodic-form]
  P0g meanSpace(mesh);
  TrialFunction u(vh);
  TrialFunction multiplier(meanSpace);
  TestFunction v(vh);
  TestFunction meanTest(meanSpace);
  RealFunction source = [](const Point& point) {
    return 2 * std::sin(point.x()) * std::cos(point.y());
  };
  Problem poisson(u, multiplier, v, meanTest);
  poisson = Integral(Grad(u), Grad(v)) + Integral(multiplier, v) + Integral(u, meanTest) -
    Integral(source, v) + PeriodicBC(u, periodicDOFs);
  //! [periodic-form]

  //! [periodic-solve]
  Solver::SparseLU solver(poisson);
  solver.solve();
  if (!solver.success())
    return EXIT_FAILURE;
  IO::XDMF xdmf("Periodic");
  xdmf.grid().setMesh(mesh).add("u", u.getSolution());
  xdmf.write().flush();
  //! [periodic-solve]
  return 0;
}
//! [periodic-source]
