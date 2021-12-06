#include <iostream>

#include <RodinExternal/MMG.h>

using namespace Rodin::Cast;
using namespace Rodin::External;

int main(int argc, char** argv)
{
  if (argc != 3)
  {
    std::cout << "Usage:\t" << argv[0]
              << " <filename>.mesh" << " <filename>.sol"
              << std::endl;
    std::exit(EXIT_FAILURE);
  }
  MMG::Mesh2D mmgMesh = MMG::Mesh2D::load(argv[1]);

  return 0;
}
