#include <Rodin/Variational.h>
#include <Rodin/Mesh.h>

using namespace Rodin;
using namespace Rodin::Variational;

int main(int, char**)
{
  const char* meshFile = "../resources/mfem/meshes/square-disc.mesh";

  // Load mesh
  Mesh Omega = Mesh::load(meshFile);
  Omega.getHandle().UniformRefinement();

  // Build finite element space
  H1 Vh(Omega);
  H1 Gh(Omega, 2);

  // Build a grid function
  GridFunction u(Vh);
  auto c = mfem::FunctionCoefficient(
      [](const mfem::Vector& v)
      {
        return v[0] * v[0] * v[0] + v[1] * v[1];
      });
  u.getHandle().ProjectCoefficient(c);

  GridFunction g(Gh);
  Gradient grad(u);
  grad.buildMFEMVectorCoefficient();

  g.getHandle().ProjectCoefficient(grad.getMFEMVectorCoefficient());

  Omega.save("Omega.mesh");
  g.save("grad.gf");
  u.save("u.gf");

  return 0;
}
