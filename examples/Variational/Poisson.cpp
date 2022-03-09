#include "Rodin/Mesh.h"
#include "Rodin/Solver.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Variational;

int main(int, char**)
{
  const char* meshFile = "../resources/mfem/meshes/poisson-example.mesh";

  // Define boundary attributes
  int Gamma = 1;

  // Load mesh
  Mesh Omega = Mesh::load(meshFile);

  // Functions
  H1 Vh(Omega);
  TrialFunction u(Vh);
  TestFunction  v(Vh);

  // Define problem
  auto f = ScalarCoefficient(1.0);
  auto g = ScalarCoefficient(0.0);

  Problem poisson(u, v);
  poisson = Inner(Grad(u), Grad(v))
          - Inner(f, v)
          + DirichletBC(u, g).on(Gamma);

  // Solve problem
  Solver::CG().setMaxIterations(200)
              .setRelativeTolerance(1e-12)
              .printIterations(true)
              .solve(poisson);

  // Save solution
  u.getGridFunction().save("u.gf");
  Omega.save("Omega.mesh");

  return 0;
}
