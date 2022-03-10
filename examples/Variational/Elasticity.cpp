/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Mesh.h>
#include <Rodin/Solver.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Variational;

int main(int argc, char** argv)
{
  const char* meshFile = "../resources/mfem/meshes/elasticity-example.mesh";

  // Define boundary attributes
  int Gamma = 1, GammaD = 2, GammaN = 3, Gamma0 = 4;

  // Load mesh
  Mesh Omega = Mesh::load(meshFile);

  // Functions
  int d = 2;
  H1 Vh(Omega, d);

  TrialFunction u(Vh);
  TestFunction  v(Vh);

  // Lamé coefficients
  auto mu     = ScalarCoefficient(0.3846),
       lambda = ScalarCoefficient(0.5769);

  auto f = VectorCoefficient{0, -1};

  // Define problem
  Problem elasticity(u, v);
  elasticity = ElasticityIntegrator(lambda, mu)
             - BoundaryIntegral(f, v).over(GammaN)
             + DirichletBC(u, VectorCoefficient{0, 0}).on(GammaD);

  // Solve problem
  Solver::CG().setMaxIterations(200)
              .setRelativeTolerance(1e-12)
              .printIterations(true)
              .solve(elasticity);

  // Save solution
  u.getGridFunction().save("u.gf");
  Omega.save("Omega.mesh");

  return 0;
}
