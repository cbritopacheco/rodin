#include "Rodin/Solver.h"
#include "Rodin/Variational.h"

using namespace Rodin;
using namespace Rodin::Variational;

int main(int argc, char** argv)
{
  const char* meshFile = "../resources/mfem/meshes/holes.mesh";

  // Define boundary attributes
  int Gamma0 = 1, GammaD = 2, GammaN = 3;

  // Load mesh
  Mesh Omega = Mesh::load(meshFile);

  // Functions
  H1 Vh(Omega, 2);
  GridFunction u(Vh), v(Vh);

  // Lamé coefficients
  double mu = 0.3846,
         lambda = 0.5769;

  // Define problem
  Problem elasticity(u, v);
  elasticity = ElasticityIntegrator(mu, lambda)
             + DirichletBC(GammaD, VectorCoefficient{0, 0})
             + NeumannBC(GammaN, VectorCoefficient{0, -1});

  // Solve problem
  Solver::PCG().setMaxIterations(200)
               .setRelativeTolerance(1e-12)
               .printIterations(true)
               .solve(elasticity);

  // Save solution
  u.save("sol.gf");
  Omega.save("mesh.mesh");

  return 0;
}
