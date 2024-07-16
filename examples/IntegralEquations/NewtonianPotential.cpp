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
Real K(const Point& x, const Point& y)
{
  return 1. / (4 * M_PI * ((x - y).norm()));
}

inline
Real exact(const Point& x)
{
  if (abs(1 - x.squaredNorm()) < 0.00001)
    return 4. / (M_PI * std::sqrt(abs(0.00001)));
  else
    return 4. / (M_PI * std::sqrt(abs(1 - x.squaredNorm())));
}

int main(int, char**)
{
  Mesh mesh;
  mesh.load("D1.o.mesh", IO::FileFormat::MEDIT);
  // mesh.load("miaow.medit.o.mesh", IO::FileFormat::MEDIT);
  mesh.getConnectivity().compute(1, 2);

  P1 fes(mesh);
  TrialFunction u(fes);
  TestFunction  v(fes);

  DenseProblem eq(u, v);
  eq = Integral(0.001 * Grad(u), Grad(v))
     + Integral(Potential(K, u), v)
     - Integral(v)
     ;

  std::cout << "assemblage\n";
  eq.assemble();

  std::cout << "save\n";
  std::ofstream file("matrix.csv");
  file << eq.getStiffnessOperator().format(
      Eigen::IOFormat(Eigen::StreamPrecision, 0, ", ", "\n"));
  file.close();

  std::cout << "resolution\n";
  Solver::LDLT(eq).solve();

  u.getSolution().save("u.gf");
  mesh.save("u.mesh");

  std::cout << "average:\n";
  u.getSolution().setWeights();
  std::cout << Integral(u.getSolution()) << std::endl;

  GridFunction one(fes);
  one = RealFunction(1);

  std::cout << "phi\n";
  GridFunction phi(fes);
  phi = Potential(K, u.getSolution());
  phi.save("phi.gf");

  std::cout << "potential\n";
  phi.setWeights();
  one.setWeights();
  std::cout << Integral(phi).compute() / Integral(one) << std::endl;

  return 0;
}
