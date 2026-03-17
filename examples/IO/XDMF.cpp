/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Math/Constants.h"
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/IO/XDMF.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  const Real dt = 0.05;
  const size_t Nt = 40;
  const size_t Nc = 32;

  // ---- mesh ---------------------------------------------------------------
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { Nc, Nc });
  mesh.scale(1.0 / (Nc - 1));

  // ---- finite element space ----------------------------------------------
  P1 fes(mesh, 2);

  // evolving field
  GridFunction u(fes);
  u.setName("u");

  // ---- XDMF writer --------------------------------------------------------
  IO::XDMF xdmf("TimeEvolution");

  auto grid = xdmf.grid();
  grid.setMesh(mesh);
  grid.add(u); // vector field

  // ---- time loop ----------------------------------------------------------

  for (size_t k = 0; k < Nt; ++k)
  {
    const Real t = k * dt;

    // update field
    u = [t](const Geometry::Point& p)
    {
      const Real x = p.x();
      const Real y = p.y();

      return Math::Vector<Real>{{
        std::sin(2.0 * Math::Constants::pi() * (x + t)),
        std::cos(2.0 * Math::Constants::pi() * (y - t))
      }};
    };

    // export snapshot
    xdmf.write(t);
  }

  return 0;
}
