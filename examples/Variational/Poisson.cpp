#include "Rodin/Mesh.h"
#include "Rodin/Solver.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Variational;

int main(int, char**)
{
  const char* meshFile = "../resources/mfem/meshes/star.mesh";

  // Define boundary attributes
  int Gamma = 1;

  // Load mesh
  Mesh Omega = Mesh::load(meshFile);

  // Functions
  H1 Vh(Omega);
  GridFunction u(Vh);

  // Define problem
  auto f = ScalarCoefficient(1.0);
  auto g = ScalarCoefficient(0.0);

  Problem poisson(u);
  poisson = DiffusionIntegrator()
          - DomainLFIntegrator(f)
          + DirichletBC(Gamma, g);

  // Solve problem
  Solver::PCG().setMaxIterations(200)
               .setRelativeTolerance(1e-12)
               .printIterations(true)
               .solve(poisson);

  // Save solution
  u.save("u.gf");
  Omega.save("Omega.mesh");

  return 0;
}
