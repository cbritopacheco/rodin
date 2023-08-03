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

static constexpr Geometry::Attribute Boundary = 2;

static constexpr Scalar m = 1;

int main(int, char**)
{
  const size_t n = 200;

  // Build a mesh
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Geometry::Triangle, n, n);
  mesh.getConnectivity().compute(1, 2);
  mesh.scale(1.0 / (n - 1));
  mesh.save("Cell.mesh");

  // Functions
  P1 vh(mesh);

  TrialFunction psi1(vh);
  TrialFunction psi2(vh);
  TestFunction  v(vh);

  // Define problem
  ScalarFunction gamma =
    [](const Point& p)
    {
      return 2 + sin(2 * M_PI * m * p.x()) * cos(2 * M_PI * m * p.y());
    };

  ScalarFunction dxgamma =
    [](const Point& p)
    {
      return 2 * M_PI * m * cos(2 * M_PI * m * p.x()) * cos(2 * M_PI * m * p.y());
    };

  ScalarFunction dygamma =
    [](const Point& p)
    {
      return -2 * M_PI * m * sin(2 * M_PI * m * p.x()) * sin(2 * M_PI * m * p.y());
    };

  GridFunction conductivity(vh);
  conductivity = gamma;
  conductivity.save("Conductivity.gf");

  ScalarFunction g = 0.0;

  IndexMap<IndexSet> dofs;
  for (Index i = 0; i < vh.getSize(); i += n)
    dofs[i].insert(i + n - 1);

  for (Index i = 0; i < n; i++)
    dofs[i].insert(i + n * (n - 1));

  Problem poisson1(psi1, v);
  poisson1 = Integral(gamma * Grad(psi1), Grad(v))
           + Integral(dxgamma, v)
           + PeriodicBC(psi1, dofs);
  poisson1.assemble();

  Problem poisson2(psi2, v);
  poisson2 = Integral(gamma * Grad(psi2), Grad(v))
           + Integral(dygamma, v)
           + PeriodicBC(psi2, dofs);
  poisson2.assemble();

  // Solve the problem
  Solver::SparseLU solver;

  poisson1.solve(solver);
  psi1.getSolution().save("Psi1.gf");

  poisson2.solve(solver);
  psi2.getSolution().save("Psi2.gf");

  P1 gh(mesh, 2);

  GridFunction psi(gh);
  psi = VectorFunction{ psi1.getSolution(), psi2.getSolution() };
  psi.save("Psi.gf");

  GridFunction gpsi_1(vh);
  gpsi_1 = gamma * Pow(Frobenius(Grad(psi1.getSolution())), 2);
  gpsi_1.setWeights();

  GridFunction gpsi_2(vh);
  gpsi_2 = gamma * Pow(Frobenius(Grad(psi2.getSolution())), 2);
  gpsi_2.setWeights();

  GridFunction Ah00(vh);
  Ah00 = [&](const Point& p) { return gamma(p) - gamma(p) * Jacobian(psi)(p).coeff(0, 0); };
  Ah00.setWeights();
  Ah00.save("dx00.gf");

  GridFunction Ah11(vh);
  Ah11 = [&](const Point& p) { return gamma(p) - gamma(p) * Jacobian(psi)(p).coeff(1, 1); };
  Ah11.setWeights();
  Ah11.save("dx11.gf");

  const Scalar value_1 = Integral(gpsi_1);
  const Scalar value_2 = Integral(gpsi_2);
  const Scalar value_00 = Integral(Ah00);
  const Scalar value_11 = Integral(Ah11);
  std::cout << value_1 << std::endl;
  std::cout << value_2 << std::endl;
  std::cout << value_00 << std::endl;
  std::cout << value_11 << std::endl;


  return 0;
}


