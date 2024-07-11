/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Math.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

const Real E = 100;
const Real nu = 0.48;
const Real lambda = (E * nu) / ((1.0 + nu) * (1.0 - 2 * nu));
const Real mu = E / (1.0 + nu);

inline
void K(Math::Matrix<Real>& res, const Point& x, const Point& y)
{
  const auto norm = (x - y).norm();
  res.resize(3, 3);
  const Real L00 = (1 - nu) / (2 * M_PI * mu * norm) + nu * (x(0) - y(0)) * (x(0) - y(0)) / (2 * M_PI * mu * norm * norm * norm);
  const Real L01 = nu * (x(0) - y(0)) * (x(1) - y(1)) / (2 * M_PI * mu * norm * norm * norm);

  const Real L10 = nu * (x(1) - y(1)) * (x(0) - y(0)) / (2 * M_PI * mu * norm * norm * norm);
  const Real L11 = (1 - nu) / (2 * M_PI * mu * norm) + nu * (x(1) - y(1)) * (x(1) - y(1)) / (2 * M_PI * mu * norm * norm * norm);

  const Real L20 = -(1 - 2 * nu) * (x(0) - y(0)) / (4 * M_PI * mu * norm * norm);
  const Real L21 = -(1 - 2 * nu) * (x(1) - y(1)) / (4 * M_PI * mu * norm * norm);
  const Real L22 = (1 - nu) / (2 * M_PI * mu * norm);

  const Real L02 = -L20;
  const Real L12 = -L21;

  res << L00, L01, L02,
         L10, L11, L12,
         L20, L21, L22;
}

int main(int, char**)
{
  // Threads::getGlobalThreadPool().reset(6);
  Mesh mesh;
  mesh.load("D1.o.mesh", IO::FileFormat::MEDIT);
  mesh.save("D1.mfem.mesh");
  mesh.getConnectivity().compute(1, 2);

  P1 fes(mesh, 3);

  VectorFunction e1{1, 1, 0};

  TrialFunction u(fes);
  TestFunction  v(fes);
  DenseProblem eq(u, v);
  eq = 0.001 * LinearElasticityIntegral(u, v)(lambda, mu)
     + Integral(Potential(K, u), v)
     - Integral(e1, v);

  std::cout << "assemblage\n";
  eq.assemble();

  // std::cout << "save\n";
  // std::ofstream file("matrix.csv");
  // file << eq.getStiffnessOperator().format(
  //     Eigen::IOFormat(Eigen::StreamPrecision, 0, ", ", "\n"));
  // file.close();

  std::cout << "resolution\n";
  Solver::CG(eq).solve();

  u.getSolution().save("u.gf");
  mesh.save("u.mesh");

  // std::cout << "average:\n";
  // u.getSolution().setWeights();
  // std::cout << Integral(u.getSolution()) << std::endl;

  std::cout << "phi\n";
  GridFunction phi(fes);
  phi = Potential(K, u.getSolution());
  phi.save("phi.gf");
  mesh.save("phi.mesh");

  // std::cout << "potential\n";
  // phi.setWeights();
  // std::cout << Integral(phi).compute() / M_PI << std::endl;

  return 0;
}

