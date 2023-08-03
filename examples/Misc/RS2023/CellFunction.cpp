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
  const size_t n = 32;

  // Build a mesh
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Geometry::Triangle, n, n);
  mesh.getConnectivity().compute(1, 2);
  mesh.scale(1.0 / (n - 1));
  mesh.save("Cell.mesh");

  Alert::Info() << "Number of mesh elements: " << n << "."
                << Alert::NewLine
                << "Saved mesh to Cell.mesh "
                << Alert::Raise;

  // Functions
  P1 vh(mesh);

  TrialFunction psi1(vh);
  TrialFunction psi2(vh);
  TestFunction  v(vh);

  // Define problem
  ScalarFunction gamma =
    [](const Point& p)
    {
      return 2 + sin(2 * M_PI * m * p.x()) * sin(2 * M_PI * m * p.y());
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

  Alert::Info() << "Saved conductivity coefficient to Conductivity.gf" << Alert::Raise;

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

  Problem poisson2(psi2, v);
  poisson2 = Integral(gamma * Grad(psi2), Grad(v))
           + Integral(dygamma, v)
           + PeriodicBC(psi2, dofs);

  // Solve the problem
  Solver::SparseLU solver;

  Alert::Info() << "Solving for first component..." << Alert::Raise;
  poisson1.solve(solver);
  psi1.getSolution().save("CellFunction1.gf");
  Alert::Success() << "Done! Saved to CellFunction1.gf" << Alert::Raise;

  Alert::Info() << "Solving for second component..." << Alert::Raise;
  poisson2.solve(solver);
  psi2.getSolution().save("CellFunction2.gf");
  Alert::Success() << "Done! Saved to CellFunction2.gf" << Alert::Raise;

  P1 gh(mesh, 2);

  GridFunction psi(gh);
  psi = VectorFunction{ psi1.getSolution(), psi2.getSolution() };
  psi.save("CellFunction.gf");
  Alert::Success() << "Saved cell function to CellFunction.gf" << Alert::Raise;

  Alert::Info() << "Computing homogenized coefficient..." << Alert::Raise;

  GridFunction Ah00(vh);
  Ah00 = [&](const Point& p) { return gamma(p) - gamma(p) * Jacobian(psi)(p).coeff(0, 0); };
  Ah00.setWeights();

  GridFunction Ah11(vh);
  Ah11 = [&](const Point& p) { return gamma(p) - gamma(p) * Jacobian(psi)(p).coeff(1, 1); };
  Ah11.setWeights();

  Alert::Success() << "Homogenized coefficient: "
                   << Alert::NewLine
                   << "[ " << std::setprecision(16) << Integral(Ah00).compute() << " 0 ]"
                   << Alert::NewLine
                   << "[ 0 " << Integral(Ah11).compute() << " ]"
                   << Alert::Raise;

  return 0;
}


