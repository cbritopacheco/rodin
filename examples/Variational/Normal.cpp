/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  constexpr size_t n = 16;
  const size_t vdim = 2;

  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, n, n);
  mesh.scale(1.0 / (n - 1));
  mesh.getConnectivity().compute(1, 2);

  mesh.displace(VectorFunction{ 0, 0.2 * sin(2 * M_PI * F::x) });

  P1 fes(mesh, vdim);
  GridFunction gf(fes);
  gf.projectOnBoundary(BoundaryNormal(mesh));

  mesh.save("Normal.mesh");
  gf.save("Normal.gf");
}


