/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  // Build a mesh
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 16, 16 });

  // Functions
  P1<Complex> vh(mesh);

  GridFunction miaow(vh);
  miaow = [](const Point& p){ return Complex(0, 1) * Complex(p.x(), 2 * p.y()); };
  miaow.setWeights();

  Complex s = Integral(miaow).compute();
  std::cout << s << std::endl;

  // TrialFunction u(vh);
  // TestFunction  v(vh);

  // // Define problem
  // Problem poisson(u, v);
  // poisson = Integral(Grad(u), Grad(v))
  //         - Integral(v)
  //         + DirichletBC(u, Zero());

  // poisson.assemble();
  // Solver::CG solver;
  // poisson.solve(solver);

  // // Save solution
  // u.getSolution().save("Poisson.gf");
  // mesh.save("Poisson.mesh");

  return 0;
}
