/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
//! [xdmf-source]
#include <cmath>
#include <cstdlib>

#include <Rodin/Geometry.h>
#include <Rodin/Math/Constants.h>
#include <Rodin/Variational.h>
#include <Rodin/IO/XDMF.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main()
{
  //! [xdmf-setup]
  constexpr size_t nodes = 32;
  constexpr size_t steps = 40;
  constexpr Real dt = 0.05;
  Mesh mesh = Mesh().UniformGrid(Polytope::Type::Triangle, {nodes, nodes});
  mesh.scale(1.0 / (nodes - 1));
  P1 vh(mesh, 2);
  GridFunction u(vh);
  //! [xdmf-setup]

  //! [xdmf-writer]
  IO::XDMF xdmf("TimeEvolution");
  xdmf.grid().setMesh(mesh).add("u", u);
  //! [xdmf-writer]

  //! [xdmf-loop]
  for (size_t k = 0; k <= steps; ++k)
  {
    const Real t = k * dt;
    u = [t](const Point& point) {
      return Math::Vector<Real>{{std::sin(2.0 * Math::Constants::pi() * (point.x() + t)),
        std::cos(2.0 * Math::Constants::pi() * (point.y() - t))}};
    };
    xdmf.write(t).flush();
  }
  //! [xdmf-loop]
  return EXIT_SUCCESS;
}
//! [xdmf-source]
