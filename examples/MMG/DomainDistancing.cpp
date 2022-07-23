/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <iostream>

#include <Rodin/Mesh.h>
#include <RodinExternal/MMG.h>

using namespace Rodin;
using namespace Rodin::Variational;
using namespace Rodin::External;

int main(int, char**)
{
  Mesh mesh;
  mesh.load("thks.mesh", IO::FileFormat::MEDIT);

  FiniteElementSpace<H1> fes(mesh);

  auto dist = MMG::Distancer(fes).distance(mesh);

  mesh = MMG::ImplicitDomainMesher().setHMax(0.05)
                                    .setRMC(1e-3)
                                    .setBoundaryReference(666)
                                    .discretize(dist);

  mesh.save("mfem.mesh");
  mesh.save("medit.mesh", IO::FileFormat::MEDIT);
  // dist.save("dist.gf");


  // GridFunction ls(fes);
  // ls.load("ls.gf");

  // auto implicitDomain = MMG::ImplicitDomainMesher().setHMax(0.02).discretize(ls);

  // implicitDomain.save("implicit.mesh");

  return 0;
}

