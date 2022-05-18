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

  MMG::Mesh2D mesh;
  mesh.load(argv[1]);

  MMG::MeshAdaptor2D adaptor;
  adaptor.setGradation(1.23).adapt(mesh);

  mesh.save("graded.mesh");

  return 0;
}
