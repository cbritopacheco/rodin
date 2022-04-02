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
  const char* meshFile = "../resources/mfem/elasticity-example.mesh";

  // Define boundary attributes
  int Gamma = 1, GammaD = 2, GammaN = 3, Gamma0 = 4;

  // Load mesh
  Mesh Omega = Mesh::load(meshFile);

  // Functions
  int d = 2;
  H1 Vh(Omega, d);

  // Lamé coefficients
  auto mu     = ScalarFunction(0.3846),
       lambda = ScalarFunction(0.5769);

  // Pull force
  auto f = VectorCoefficient{0, -1};

  // Define problem
  TrialFunction u(Vh);
  TestFunction  v(Vh);
  Problem elasticity(u, v);
  elasticity = Integral(lambda * Div(u), Div(v))
             + Integral(
                 mu * (Jacobian(u) + Jacobian(u).T()), 0.5 * (Jacobian(v) + Jacobian(v).T()))
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
