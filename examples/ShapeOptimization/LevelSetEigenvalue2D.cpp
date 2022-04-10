#include <Rodin/Mesh.h>
#include <Rodin/Solver.h>
#include <Rodin/Variational.h>
#include <RodinExternal/MMG.h>

// Dependencies required for eigenvalue computation
#include <Eigen/Core>
#include <Spectra/GenEigsSolver.h>
#include <Spectra/MatOp/SparseGenMatProd.h>

using namespace Rodin;
using namespace Rodin::Variational;
using namespace Rodin::External;

int main(int argc, char** argv)
{
  const char* meshFile = "../resources/mfem/poisson-example.mesh";

  // Define boundary attributes
  int Gamma = 1;

  // Load mesh
  Mesh Omega;
  Omega.load(meshFile);

  // You can do a uniform refinement if you want
  Omega.refine();

  // Functions
  FiniteElementSpace<H1> Vh(Omega);
  TrialFunction u(Vh);
  TestFunction  v(Vh);

  // Define problem
  auto sigma = ScalarFunction(1.0);

  // A - sigma * B
  Problem op(u, v);
  op = Integral(Grad(u), Grad(v))
       - Integral(sigma * u, v)
       + DirichletBC(u, ScalarFunction(0.0)).on(Gamma);
  op.update().assemble();

  // B
  Problem b(u, v);
  b = Integral(u, v);
  b.update().assemble();


  auto& m1 = op.getStiffnessMatrix();
  auto& m2 = b.getStiffnessMatrix();

  // 1. Convert m1 and m2 to Eigen::SparseMatrix
  // 2. Use Spectra with the shift-inverse method to compute eigenvalue
  // 3. ???
  // 4. Profit

  return 0;
}
