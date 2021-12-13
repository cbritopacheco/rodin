#include <Rodin/Variational.h>
#include <Rodin/Mesh.h>

using namespace Rodin;
using namespace Rodin::Variational;

int main(int, char**)
{
  const char* meshFile = "../resources/mfem/meshes/star.mesh";

  // Load mesh
  Mesh Omega = Mesh::load(meshFile);

  // Build finite element space
  H1 Vh(Omega);

  // Build a grid function
  GridFunction u(Vh);

  return 0;
}
