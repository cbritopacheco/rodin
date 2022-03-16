#include <Rodin/Mesh.h>
#include <RodinExternal/MMG.h>

using namespace Rodin;
using namespace Rodin::External;

int main(int argc, char** argv)
{
  if (argc != 2)
  {
    std::cout << "Usage:\t" << argv[0] << " <filename>.mesh" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  Rodin::Mesh rodinMesh = Rodin::Mesh::load(argv[1]);

  MMG::SurfaceMesh mmgMesh = Cast(rodinMesh).to<MMG::SurfaceMesh>();

  mmgMesh.save("mmg.mesh");

  return 0;
}
