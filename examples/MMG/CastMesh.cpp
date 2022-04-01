#include <iostream>

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

  Mesh rodinMesh = Mesh::load(argv[1]);
  switch (rodinMesh.getDimension())
  {
    case 2:
    {
      switch (rodinMesh.getSpaceDimension())
      {
        case 2:
        {
          MMG::Mesh2D mmgMesh = Cast(rodinMesh).to<MMG::Mesh2D>();
          mmgMesh.save("mmg.mesh");
          break;
        }
        case 3:
        {
          MMG::MeshS mmgMesh = Cast(rodinMesh).to<MMG::MeshS>();
          mmgMesh.save("mmg.mesh");
          break;
        }
        default:
          Alert::Exception("Bad space dimension").raise();
      }
      break;
    }
    case 3:
    {
      MMG::Mesh3D mmgMesh = Cast(rodinMesh).to<MMG::Mesh3D>();
      mmgMesh.save("mmg.mesh");
      break;
    }
    default:
      Alert::Exception("Bad mesh dimension").raise();
  }

  return 0;
}


