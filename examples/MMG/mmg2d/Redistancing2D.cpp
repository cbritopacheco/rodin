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
using namespace Rodin::External::MMG;

int main(int argc, char** argv)
{
  auto box = Mesh2D::load(argv[1]);
  auto ls  = ScalarSolution2D(box).load(argv[2]);

  Distancer2D().redistance(ls);

  ls.save("ls.sol");
  box.save("ls.sol.mesh");

  return 0;
}
