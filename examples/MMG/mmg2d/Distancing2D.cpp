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
  if (argc == 2)
  {
    // Generation of distance function to a subdomain
    int Interior = 1;

    Mesh2D box;
    box.load(argv[1]);
    auto ls = Distancer2D().setInteriorDomain({ Interior }).distance(box);

    ls.save("ls.sol");
    box.save("ls.sol.mesh");
  }
  else if (argc == 3)
  {
    // Generation of distance function with explicit contour
    Mesh2D box;
    box.load(argv[1]);

    Mesh2D contour;
    contour.load(argv[2]);
    auto ls = Distancer2D().enableScaling(false).distance(box, contour);

    ls.save("ls.sol");
    box.save("ls.sol.mesh");
  }
  else
  {
    Alert::Exception() << "Incorrect number of parameters.\n"
                       << "Usage: " << argv[0] << " box.mesh [contour.mesh]"
                       << Alert::Raise;
  }

  return 0;
}
