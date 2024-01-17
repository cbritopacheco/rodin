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

const Scalar lambda = 0.5769, mu = 0.3846;
const Scalar nu = lambda / (2 * (lambda + mu));

inline
void K(Math::Matrix& res, const Point& x, const Point& y)
{
  const auto norm = (x - y).norm();
  res.resize(3, 3);
  const Scalar L00 = (1 - nu) / (2 * M_PI * mu * norm) + nu * (x(0) - y(0)) * (x(0) - y(0)) / (2 * M_PI * mu * norm * norm * norm);
  const Scalar L01 = nu * (x(0) - y(0)) * (x(1) - y(1)) / (2 * M_PI * mu * norm * norm * norm);
  const Scalar L02 = nu * (x(0) - y(0)) * (x(2) - y(2)) / (2 * M_PI * mu * norm * norm * norm);

  const Scalar L10 = nu * (x(1) - y(1)) * (x(0) - y(0)) / (2 * M_PI * mu * norm * norm * norm);
  const Scalar L11 = (1 - nu) / (2 * M_PI * mu * norm) + nu * (x(1) - y(1)) * (x(1) - y(1)) / (2 * M_PI * mu * norm * norm * norm);
  const Scalar L12 = nu * (x(1) - y(1)) * (x(2) - y(2)) / (2 * M_PI * mu * norm * norm * norm);

  const Scalar L20 = (1 - 2 * nu) * (x(0) - y(0)) / (4 * M_PI * mu * norm * norm);
  const Scalar L21 = (1 - 2 * nu) * (x(1) - y(1)) / (4 * M_PI * mu * norm * norm);
  const Scalar L22 = (1 - nu) / (2 * M_PI * mu * norm);

  res << L00, L01, L02,
         L10, L11, L12,
         L20, L21, L22;
}

int main(int, char**)
{
  Mesh mesh;
  mesh.load("D1.mesh");
  mesh.getConnectivity().compute(1, 2);

  P1 fes(mesh, 3);

  VectorFunction e1{0, 0, 1};

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
  Solver::CG<Math::Matrix, Math::Vector> cg;
  eq.solve(cg);

  u.getSolution().save("u.gf");
  mesh.save("u.mesh");

  // std::cout << "average:\n";
  // u.getSolution().setWeights();
  // std::cout << Integral(u.getSolution()) << std::endl;

  std::cout << "phi\n";
  GridFunction phi(fes);
  phi = Potential(K, u.getSolution());
  phi.save("phi.gf");

  // std::cout << "potential\n";
  // phi.setWeights();
  // std::cout << Integral(phi).compute() / M_PI << std::endl;

  return 0;
}

