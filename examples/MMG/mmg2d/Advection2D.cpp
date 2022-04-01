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
  if (argc == 5)
  {
    auto mesh = Mesh2D::load(argv[1]);
    auto ls = ScalarSolution2D::load(argv[2]).setMesh(mesh);
    auto displacement = VectorSolution2D::load(argv[3]).setMesh(mesh);
    double dt = std::atof(argv[4]);

    Advect2D(ls, displacement).step(dt);

    ls.save("advected.sol");
    mesh.save("advected.sol.mesh");
  }
  else
  {
    Alert::Exception() << "Incorrect number of parameters.\n"
                       << "Usage: " << argv[0] << " omega.mesh ls.sol disp.sol dt"
                       << Alert::Raise;
  }
  return 0;
}
