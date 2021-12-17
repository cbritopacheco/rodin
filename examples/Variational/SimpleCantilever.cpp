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
  elasticity = ElasticityIntegrator(lambda, mu)
             + DirichletBC(GammaD, VectorCoefficient{0, 0})
             + NeumannBC(GammaN, VectorCoefficient{0, -1});

  // Hilbert regularization
  auto e = ScalarCoefficient(0.5) * (Jacobian(u) + Jacobian(u).T());
  auto Ae = ScalarCoefficient(2.0) * mu * e + lambda * Trace(e) * IdentityMatrix(d);

  H1 Sh(Omega);
  GridFunction one(Sh);
  GridFunction c(Sh);
  one.getHandle() = 1.0;


  // Solve problem
  Solver::PCG().setMaxIterations(200)
               .setRelativeTolerance(1e-12)
               .printIterations(true)
               .solve(elasticity);

  LinearForm lf(Sh);
  DomainLFIntegrator Compliance(Dot(Ae, e));
  Compliance.setLinearForm(lf);
  Compliance.eval();

  BilinearForm bf(Vh);
  ElasticityIntegrator elas(lambda, mu);
  elas.setBilinearForm(bf);
  elas.eval();

  std::cout << lf(one) << std::endl;
  std::cout << bf(u, u) << std::endl;


  Omega.save("Omega.mesh");
  u.save("u.gf");
  c.save("c.gf");

  return 0;
}
