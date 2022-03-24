#include <Rodin/Mesh.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Variational;

int main(int, char**)
{
  const char* meshFile = "rodin.mesh";

  Mesh Omega = Mesh::load(meshFile);
  H1 Vh(Omega, 2);
  H1 Ph(Omega);

  GridFunction u(Vh);
  GridFunction der(Ph);

  Derivative dx(0, 0, u);

  der.project(dx);

  return 0;
}
