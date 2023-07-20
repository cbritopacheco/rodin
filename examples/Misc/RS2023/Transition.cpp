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

static constexpr Attribute dQ = 2;
static constexpr Attribute dR = 3;

int main(int, char**)
{
  // Build a mesh
  Mesh mesh;
  mesh.load("Q.medit.mesh", IO::FileFormat::MEDIT);
  mesh.getConnectivity().compute(1, 2);

  // Functions
  P1 vh(mesh);

  const Scalar epsilon = 0.01;
  const Scalar R = 0.5;
  assert(R > epsilon);
  const Math::Vector x0{{0.5, 0.5}};

  // Define problem
  const auto f =
    [&](const Scalar& x)
    {
      if (x > 0)
        return std::exp(-1.0 / x);
      else
        return 0.0;
    };
  const auto g =
    [&](const Scalar& x)
    {
      return f(x) / (f(x) + f(1 - x));
    };
  const auto h =
    [&](const Geometry::Point& p)
    {
      const Scalar r = (p.getCoordinates() - x0).norm();
      return g((r - epsilon) / (R - epsilon));
    };

  GridFunction transition(vh);
  transition = h;
  transition.save("F.gf");
  mesh.save("F.mesh");

  return 0;
}


