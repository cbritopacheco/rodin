#include <iostream>

#include "RodinExternal/MMG.h"

using namespace Rodin;
using namespace Rodin::External;

int main(int argc, char** argv)
{
  const char* meshFile = "../resources/mfem/meshes/hole.mesh";

  mfem::OptionsParser args(argc, argv);
  args.AddOption(&meshFile, "-m", "--mesh", "Mesh file to use.");
  args.ParseCheck();

  Rodin::Mesh rodinMesh = Rodin::Mesh::load(meshFile);

  MMG::Mesh2D mmgMesh = Cast(rodinMesh).to<MMG::Mesh2D>();

  mmgMesh.save("mmg.mesh");

  return 0;
}
