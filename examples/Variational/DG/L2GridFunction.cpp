/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  const char* meshFile = "../resources/mfem/elasticity-example.mesh";
  Mesh Omega;
  Omega.load(meshFile);
  L2 Vh(Omega);
  GridFunction u(Vh);
  u = ScalarFunction(
      [](const Point& v) -> double
      {
        return cos(v.x());
      });
  Omega.save("l2.mesh");
  u.save("l2.gf");
  return 0;
}
