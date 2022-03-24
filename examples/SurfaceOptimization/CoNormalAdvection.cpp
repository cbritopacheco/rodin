#include <Rodin/Mesh.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Variational;

int main(int, char**)
{
  const char* meshFile = "../resources/mfem/meshes/poisson-example.mesh";

  Mesh Omega = Mesh::load(meshFile);

  return 0;
}
