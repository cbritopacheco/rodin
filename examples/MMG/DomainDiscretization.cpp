/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Geometry.h>
#include <RodinExternal/MMG.h>

using namespace Rodin;
using namespace Rodin::External;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

static constexpr Attribute interior = 1;
static constexpr Attribute exterior = 2;
static constexpr Scalar hmax = 0.1;
static constexpr Scalar radius = 0.1;

int main(int, char**)
{
  const size_t n = 16;
  MMG::Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Geometry::Triangle, n, n);
  mesh.scale(1. / (n - 1));

  P1 fes(mesh);
  MMG::ScalarGridFunction gf(fes);
  gf = [](const Geometry::Point& p)
       {
         return (p.x() - 0.5) * (p.x() - 0.5) + (p.y() - 0.5) * (p.y() - 0.5) - radius;
       };

  gf.save("LevelSet.gf");
  mesh.save("Domain.mesh");

  MMG::ImplicitDomainMesher().split(RODIN_DEFAULT_POLYTOPE_ATTRIBUTE, { interior, exterior })
                             .setBoundaryReference(RODIN_DEFAULT_POLYTOPE_ATTRIBUTE)
                             .setHMax(hmax)
                             .discretize(gf)
                             .save("Discretized.mesh", IO::FileFormat::MEDIT);
  return 0;
}
