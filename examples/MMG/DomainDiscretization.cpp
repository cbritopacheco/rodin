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

int main(int, char**)
{
  const size_t n = 16;
  MMG::Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Geometry::Triangle, n, n);
  mesh.getConnectivity().compute(1, 2);

  P1 fes(mesh);
  MMG::ScalarGridFunction gf(fes);
  gf = [](const Geometry::Point& p)
       {
         if (p.x() < 8)
           return -1;
         else return 1;
       };

  gf.save("LevelSet.gf");
  mesh.save("Domain.mesh");

  MMG::ImplicitDomainMesher().split(1, { 1, 2 })
                             .setHMax(1)
                             .discretize(gf)
                             .save("Discretized.mesh");

  return 0;
}
