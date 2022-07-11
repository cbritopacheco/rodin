/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Mesh.h>
#include <Rodin/Solver.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Variational;

int main(int, char**)
{
  Mesh mesh;
  mesh.load("dmmg.mesh", IO::MeshFormat::MEDIT);

  FiniteElementSpace<H1> fes(mesh, 3);
  GridFunction gf(fes);
  gf.load("dmmg.sol", IO::GridFunctionFormat::MEDIT);

  gf.save("miaow.gf");

  mesh.save("OmegaMedit.mesh");

  return 0;
}
