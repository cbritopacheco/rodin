#include <iostream>

#include "RodinExternal/MMG.h"

using namespace Rodin::External;

int main(int argc, char** argv)
{
  if (argc != 2)
  {
    std::cout << "Usage:\t" << argv[0] << " <filename>.mesh" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  MMG::Mesh2D mesh = MMG::Mesh2D::load(argv[1]);

  MMG::MeshOptimizer2D optimizer;
  optimizer.optimize(mesh);

  mesh.save("optimized.mesh");

  return 0;
}
