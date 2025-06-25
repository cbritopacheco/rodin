/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Assembly.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  auto pi = Rodin::Math::Constants::pi();
  size_t M = 4;

  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { M, M });

      mesh.scale(1.0 / (M - 1));
      mesh.getConnectivity().compute(1, 2);

  P1 vh(mesh, mesh.getSpaceDimension());
  VectorFunction f{
    (pi * pi - 1) * sin(pi * F::x) * exp(F::y),
    (pi * pi - 1) * sin(pi * F::y) * exp(F::x)
  };

  TrialFunction u(vh);
  TestFunction  v(vh);

  // Apply Dirichlet conditions on the entire boundary.
  Problem poisson(u, v);
  poisson = Integral(Jacobian(u), Jacobian(v))
          - Integral(f, v)
          + DirichletBC(u, VectorFunction{
              sin(pi * F::x) * exp(F::y),
              sin(pi * F::y) * exp(F::x)
            });
  poisson.assemble();
  // CG(poisson).solve();

  std::cout << poisson.getLinearSystem().getOperator() << std::endl;
  std::cout << "==============" << std::endl;
  std::cout << poisson.getLinearSystem().getVector() << std::endl;

  // VectorFunction solution{
  //   sin(pi * F::x) * exp(F::y),
  //   sin(pi * F::y) * exp(F::x)
  // };

  // // Save solution
  // u.getSolution().save("Poisson.gf");
  // mesh.save("Poisson.mesh");

  return 0;
}
