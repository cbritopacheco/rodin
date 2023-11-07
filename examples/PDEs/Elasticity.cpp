/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Variational/LinearElasticity.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int argc, char** argv)
{
  const char* meshFile = "../resources/examples/PDEs/Elasticity.mfem.mesh";

  // Define boundary attributes
  Attribute Gamma = 1, GammaD = 2, GammaN = 3, Gamma0 = 4;

  // Load mesh
  Mesh mesh;
  mesh.load(meshFile);
  mesh.getConnectivity().compute(1, 2);

  // Functions
  size_t d = mesh.getSpaceDimension();
  P1 fes(mesh, d);

  // Lamé coefficients
  Scalar lambda = 0.5769, mu = 0.3846;

  // Pull force
  VectorFunction f{0, -1};

  // Define problem
  TrialFunction u(fes);
  TestFunction  v(fes);

  Problem elasticity(u, v);
  elasticity = Integral(lambda * Div(u), Div(v))
             + Integral(
                 mu * (Jacobian(u) + Jacobian(u).T()), 0.5 * (Jacobian(v) + Jacobian(v).T()))
             - BoundaryIntegral(f, v).over(GammaN)
             + DirichletBC(u, VectorFunction{0, 0}).on(GammaD);

  // Solver object
  Solver::CG solver;
  elasticity.solve(solver);

  // Save solution
  u.getSolution().save("Elasticity.gf");
  mesh.save("Elasticity.mesh");

  return 0;
}
