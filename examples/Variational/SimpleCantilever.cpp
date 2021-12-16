#include <Rodin/Mesh.h>
#include <Rodin/Solver.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Variational;

int main(int, char**)
{
  const char* meshFile = "../resources/mfem/meshes/holes.mesh";

  // Define boundary attributes
  int Gamma0 = 1, GammaD = 2, GammaN = 3;

  // Load mesh
  Mesh Omega = Mesh::load(meshFile);

  // Build finite element space
  int d = 2;
  H1 Vh(Omega, d);

  // Lamé coefficients
  auto mu     = ScalarCoefficient(0.3846),
       lambda = ScalarCoefficient(0.5769);

  // Displacement
  GridFunction u(Vh);

  // Shape gradient
  GridFunction g(Vh);

  // Elasticity equation
  Problem elasticity(u);
  elasticity = ElasticityIntegrator(mu, lambda)
             + DirichletBC(GammaD, VectorCoefficient{0, 0})
             + NeumannBC(GammaN, VectorCoefficient{0, -1});

  // Hilbert regularization
  Problem hilbert(g);
  auto e = ScalarCoefficient(0.5) * (Jacobian(u) + Jacobian(u).T());
  auto Ae = ScalarCoefficient(2.0) * mu * e + lambda * Trace(e) * IdentityMatrix(d);
  auto Compliance = Dot(Ae, e);

  // Solve problem
  Solver::PCG().setMaxIterations(200)
               .setRelativeTolerance(1e-12)
               .printIterations(true)
               .solve(elasticity);

  Omega.save("Omega.mesh");
  u.save("u.gf");

  return 0;
}
