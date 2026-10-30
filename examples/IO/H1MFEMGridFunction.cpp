/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <iostream>

#include <Rodin/Geometry.h>
#include <Rodin/IO.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  Mesh mesh = LocalMesh::UniformGrid(Polytope::Type::Wedge, { 8, 8, 8 });
  mesh.scale(1.0 / 7.0);
  mesh.getConnectivity().compute(3, 1);
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 0);

  constexpr size_t order = 5;

  H1 scalarFES(std::integral_constant<size_t, order>{}, mesh);
  GridFunction scalarField(scalarFES);
  scalarField.project(RealFunction([](const Point& p)
  {
    return 1.0 + p.x() * p.x() - 0.5 * p.y() + p.x() * p.y() + p.z();
  }));

  using VectorFES = H1<order, Math::SpatialVector<Real>>;
  VectorFES vectorFES(std::integral_constant<size_t, order>{}, mesh, 3);
  GridFunction<VectorFES, Math::Vector<Real>> vectorField(vectorFES);
  vectorField.project(VectorFunction{
    [](const Point& p) { return p.x() * p.x() - p.y() * p.z(); },
      [](const Point& p) { return p.x() + p.y() * p.y() * p.z(); },
      [](const Point& p) { return p.z() - p.x() * p.y(); }});

  //! [h1-mfem-mesh]
  {
    std::ofstream output("H1MFEMGridFunction.mesh");
    output.precision(std::numeric_limits<Real>::max_digits10);
    IO::MeshPrinter<IO::FileFormat::MFEM, Context::Local> printer(mesh);
    printer.print(output);
    if (!output)
      return EXIT_FAILURE;
  }
  //! [h1-mfem-mesh]
  scalarField.save("H1MFEMGridFunction.scalar.gf", IO::FileFormat::MFEM);
  vectorField.save("H1MFEMGridFunction.vector.gf", IO::FileFormat::MFEM);

  Mesh loaded;
  loaded.load("H1MFEMGridFunction.mesh", IO::FileFormat::MFEM);
  loaded.getConnectivity().compute(3, 1);
  loaded.getConnectivity().compute(2, 1);
  loaded.getConnectivity().compute(1, 0);
  H1 restoredScalarSpace(std::integral_constant<size_t, order>{}, loaded);
  VectorFES restoredVectorSpace(std::integral_constant<size_t, order>{}, loaded, 3);
  GridFunction restoredScalar(restoredScalarSpace);
  GridFunction<VectorFES, Math::Vector<Real>> restoredVector(restoredVectorSpace);
  restoredScalar.load("H1MFEMGridFunction.scalar.gf", IO::FileFormat::MFEM);
  restoredVector.load("H1MFEMGridFunction.vector.gf", IO::FileFormat::MFEM);
  Real maximumError = 0;
  for (auto cell = loaded.getCell(); !cell.end(); ++cell)
  {
    const Point point(*cell, Math::SpatialPoint{1.0 / 3, 1.0 / 3, 0.5});
    const Real scalar =
      1 + point.x() * point.x() - 0.5 * point.y() + point.x() * point.y() + point.z();
    const Math::SpatialVector<Real> vector{point.x() * point.x() - point.y() * point.z(),
      point.x() + point.y() * point.y() * point.z(), point.z() - point.x() * point.y()};
    if (!std::isfinite(restoredScalar(point)) ||
      !std::isfinite(restoredVector(point).norm()))
      return EXIT_FAILURE;
    maximumError = std::max(maximumError, std::abs(restoredScalar(point) - scalar));
    maximumError = std::max(maximumError, (restoredVector(point) - vector).norm());
  }
  std::cout << "Maximum restored field error: " << maximumError << std::endl;
  constexpr Real tolerance = 1e-10;
  return maximumError < tolerance ? EXIT_SUCCESS : EXIT_FAILURE;
}
