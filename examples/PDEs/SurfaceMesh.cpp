/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [surfacemesh-source]
#include <cmath>
#include <cstdlib>
#include <charconv>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string_view>

#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <Rodin/IO/XDMF.h>
#include <Rodin/QF/QuadratureFormula.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int argc, char** argv)
{
  size_t nodes = 11;
  if (argc == 2)
  {
    const std::string_view argument(argv[1]);
    const auto parsed =
      std::from_chars(argument.data(), argument.data() + argument.size(), nodes);
    if (parsed.ec != std::errc{} || parsed.ptr != argument.data() + argument.size() ||
      nodes < 3)
      return EXIT_FAILURE;
  }
  else if (argc != 1)
    return EXIT_FAILURE;
  //! [surface-mesh]
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
  //! [surface-solve]
  //! [surface-error]
  const auto& system = reactionDiffusion.getLinearSystem();
  const Real residual =
    (system.getOperator() * system.getSolution() - system.getVector()).norm() /
    std::max(Real(1), system.getVector().norm());
  constexpr size_t ErrorQuadratureOrder = 6;
  Real squaredError = 0, squaredReference = 0;
  for (auto cell = mesh.getCell(); cell; ++cell)
  {
    const auto& formula =
      QF::QuadratureFormula::get(ErrorQuadratureOrder, cell->getGeometry());
    const auto& quadrature = cell->getQuadrature(formula);
    for (size_t q = 0; q < quadrature.getSize(); ++q)
    {
      const auto& point = quadrature.getPoint(q);
      const Real reference = point.x() * point.y();
      const Real difference = u.getSolution().getValue(point) - reference;
      const Real weight = formula.getWeight(q) * point.getDistortion();
      squaredError += weight * difference * difference;
      squaredReference += weight * reference * reference;
    }
  }
  const Real relativeError = std::sqrt(squaredError / squaredReference);
  //! [surface-error]
  std::ofstream statistics("SurfaceMesh.csv");
  statistics << std::setprecision(17)
             << "nodes,vertices,triangles,area,relative_l2_error,relative_residual\n"
             << nodes << ',' << mesh.getVertexCount() << ',' << mesh.getCellCount() << ','
             << mesh.getArea() << ',' << relativeError << ',' << residual << '\n';
  statistics.flush();
  constexpr Real ResidualTolerance = Real(1e-9);
  if (!std::isfinite(relativeError) || !std::isfinite(residual) ||
    residual > ResidualTolerance || !statistics)
    return EXIT_FAILURE;
  std::cout << "nodes=" << nodes << " relative_l2_error=" << relativeError
            << " relative_residual=" << residual << '\n';
  //! [surface-export]
  IO::XDMF xdmf("SurfaceMesh");
  xdmf.grid().setMesh(mesh).add("u", u.getSolution());
  xdmf.write().close();
  //! [surface-export]
  return 0;
}
//! [surfacemesh-source]
