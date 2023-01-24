#include <Rodin/Geometry.h>
#include <Rodin/Solver.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  const char* meshFile = "../resources/mfem/surface-meshes-example.mesh";

  // Load mesh
  Mesh Omega;
  Omega.load(meshFile);

  // Functions
  H1 Vh(Omega);
  TrialFunction u(Vh);
  TestFunction  v(Vh);

  // Right hand side
  auto f = ScalarFunction(
    [](const Point& x)
    {
      double l2 = x(0) * x(0) + x(1) * x(1) + x(2) * x(2);
      return 7 * x(0) * x(1) / l2;
    });

  Solver::CG cg;
  cg.setMaxIterations(200).setRelativeTolerance(1e-12).printIterations(true);

  // Reaction-diffusion steady state
  Problem rd(u, v);
  rd = Integral(Grad(u), Grad(v))
    + Integral(u, v)
    - Integral(f, v);
  rd.solve(cg);

  // Save solution
  u.getSolution().save("u.gf");
  Omega.save("Omega.mesh");
}
