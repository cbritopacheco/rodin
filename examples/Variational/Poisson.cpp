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
  GridFunction u(Vh), v(Vh);

  // Define problem
  Problem poisson(u, v);
  poisson = DiffusionIntegrator()
          - DomainLFIntegrator(ScalarCoefficient(1.0))
          + DirichletBC(Gamma, 0.0);

  // Solve problem
  Solver::PCG().setMaxIterations(200)
               .setRelativeTolerance(1e-12)
               .printIterations(true)
               .solve(poisson);

  // Save solution
  u.save("sol.gf");
  Omega.save("mesh.mesh");

  return 0;
}
