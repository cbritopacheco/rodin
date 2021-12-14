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
  H1 Vh(Omega, 2);

  // Build a grid function
  GridFunction u(Vh);

  // mfem::FunctionCoefficient f([](const mfem::Vector& v){ return v[0] * v[0]; });
  // u.getHandle().ProjectCoefficient(f);

  // auto tmp = Dy(u);
  // tmp.buildMFEMCoefficient();

  // GridFunction du(Vh);
  // du.getHandle().ProjectCoefficient(tmp.getMFEMCoefficient());
  auto sum = ScalarCoefficient(1.0) + ScalarCoefficient(2.);
  sum.buildMFEMCoefficient();
  auto e = Gradient(u) + Gradient(u).T();
  e.buildMFEMMatrixCoefficient();

  u.getHandle().ProjectCoefficient(sum.getMFEMCoefficient());


  Omega.save("Omega.mesh");
  u.save("test.gf");


  return 0;
}
