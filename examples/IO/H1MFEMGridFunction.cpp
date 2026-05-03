/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
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
    [](const Point& p) { return p.x() + p.y() * p.y() * p.z(); }
  });

  mesh.save("H1MFEMGridFunction.mesh", IO::FileFormat::MFEM);
  scalarField.save("H1MFEMGridFunction.scalar.gf", IO::FileFormat::MFEM);
  vectorField.save("H1MFEMGridFunction.vector.gf", IO::FileFormat::MFEM);

  std::cout
    << "Wrote H1MFEMGridFunction.mesh\n"
    << "Wrote H1MFEMGridFunction.scalar.gf\n"
    << "Wrote H1MFEMGridFunction.vector.gf\n";
}
