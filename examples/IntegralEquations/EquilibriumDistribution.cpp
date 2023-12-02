/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

inline
Scalar K(const Point& x, const Point& y)
{
  return 1. / (4 * M_PI * (x - y).norm());
}

inline
Scalar exact(const Point& x)
{
  return 4. / (M_PI * std::sqrt(1 - x.squaredNorm()));
}

int main(int, char**)
{
  Mesh mesh;
  mesh.load("D1.mesh");
  // mesh.save("miaow.medit.mesh", IO::FileFormat::MEDIT);
  mesh.getConnectivity().compute(1, 2);

  P1 fes(mesh);

  TrialFunction u(fes);
  TestFunction  v(fes);
  DenseProblem eq(u, v);
  eq = Integral(Potential(K, u), v)
     - Integral(v);

  std::cout << "assemblage\n";
  eq.assemble();

  std::cout << "resolution\n";
  Solver::LDLT solver;
  eq.solve(solver);

  u.getSolution().save("u.gf");
  mesh.save("u.mesh");

  // GridFunction phi(fes);
  // phi = [](const Point& p) { return 4. / (M_PI * std::sqrt(1 - p.squaredNorm())); };
  // phi.projectOnBoundary(ScalarFunction(0));
  // phi.save("phi.gf");

  return 0;
}
