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

  mfem::FunctionCoefficient f([](const mfem::Vector& v){ return v[0] * v[0]; });
  u.getHandle().ProjectCoefficient(f);

  auto tmp = Dy(u);
  tmp.buildMFEMCoefficient();

  GridFunction du(Vh);
  du.getHandle().ProjectCoefficient(tmp.getMFEMCoefficient());

  Omega.save("Omega.mesh");
  du.save("test.gf");


  return 0;
}
