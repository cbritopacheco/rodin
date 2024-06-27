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
  const size_t n = 32;
  const size_t vdim = 2;
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { n, n });
  mesh.scale(1.0 / (n - 1.0));

  P1 fes(mesh, vdim);
  GridFunction gf(fes);

  gf = [](const Geometry::Point& p) { return Math::Vector<Scalar>{{ p.y() - 0.5, -p.x() + 0.5 }}; };

  mesh.save("Projection.mesh");
  gf.save("Projection.gf");
}

