/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 16, 16 });
  mesh.getConnectivity().compute(1, 2);

  H1 vh(std::integral_constant<size_t, 12>{}, mesh);
  // P1 vh(mesh);
  GridFunction u(vh);

  u = [](const Geometry::Point& p)
  {
    return std::sin(M_PI * p.x()) * std::cos(M_PI * p.y());
    // return p.x() * (1 - p.x()) * p.y() * (1 - p.y());
  };
  std::cout << "Miaow\n";
  u.save("Poisson.gf");
  // u.load("Poisson.gf");
  // u.save("Poisson2.gf");
  mesh.save("Poisson.mesh");

  // P1 vh(mesh);

  // TrialFunction u(vh);
  // TestFunction  v(vh);

  // RealFunction f = 1;

  // // Apply Dirichlet conditions on the entire boundary.
  // Problem poisson(u, v);
  // poisson = Integral(Grad(u), Grad(v))
  //         - Integral(f, v)
  //         + DirichletBC(u, Zero());
  // CG(poisson).solve();

  // // Save solution
  // u.getSolution().save("Poisson.gf");
  // mesh.save("Poisson.mesh");

  return 0;
}
