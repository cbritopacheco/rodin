/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Variational/LinearElasticity.h>>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int argc, char** argv)
{
  const char* meshFile = "/Users/carlos/Projects/rodin/resources/mfem/elasticity-example.mesh";

  // Define boundary attributes
  int Gamma = 1, GammaD = 2, GammaN = 3, Gamma0 = 4;

  // Load mesh
  Mesh Omega;
  Omega.load(meshFile);

  // Functions
  size_t d = Omega.getSpaceDimension();
  H1 vh(Omega, d);

  // Lamé coefficients
  Scalar lambda = 0.5769, mu = 0.3846;

  // Pull force
  VectorFunction f{0, -1};

  // Solver object
  Solver::CG solver;
  mfem::GSSmoother smooth;
  solver.printIterations(false).setPreconditioner(smooth);

  // Define problem
  TrialFunction u(vh);
  TestFunction  v(vh);

  Problem elasticity(u, v);
  // elasticity = Integral(lambda * Div(u), Div(v))
  //            + Integral(mu * (Jacobian(u) + Jacobian(u).T()), 0.5 * (Jacobian(v) + Jacobian(v).T()))
  //            - BoundaryIntegral(f, v).over(GammaN)
  //            + DirichletBC(u, VectorFunction{0, 0}).on(GammaD);


  elasticity = LinearElasticityIntegral(u, v)(lambda, mu)
             - BoundaryIntegral(f, v).over(GammaN)
             + DirichletBC(u, VectorFunction{0, 0}).on(GammaD);

  for (int i = 0; i < 50; i++)
    elasticity.assemble();

  // elasticity.solve(solver);


  // Save solution
  // u.getSolution().save("u.gf");
  // Omega.save("Omega.mesh");

  return 0;
}
