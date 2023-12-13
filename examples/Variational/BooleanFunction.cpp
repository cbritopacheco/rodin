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
using namespace Rodin::Math::Constants;

int main(int, char**)
{
  constexpr size_t n = 16;
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, n, n);
  mesh.scale(1.0 / (n - 1));
  mesh.save("Square.mesh");

  P1 fes(mesh);
  GridFunction gf(fes);
  const auto f = cos(2 * pi() * F::x + pi() / 2) * sin(2 * pi() * F::y);
  const auto g = sin(2 * pi() * F::x) * sin(2 * pi() * F::y);

  gf = f > g;

  gf.save("GT.gf");

  gf = f < g;

  gf.save("LT.gf");

  gf = (f > g) || (f < g);

  gf.save("OR.gf");

  gf = (f > g) && (f < g);

  gf.save("AND.gf");

  return 0;
}



