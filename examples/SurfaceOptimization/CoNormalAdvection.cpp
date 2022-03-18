#include <Rodin/Mesh.h>
#include <Rodin/Variational.h>

using namespace Rodin;

int main(int, char**)
{
  const char* meshFile = "rodin.mesh";

  Mesh Omega = Mesh::load(meshFile);
  SubMesh skin = Omega.skin();
  skin.save("skin.mesh");

  return 0;
}
