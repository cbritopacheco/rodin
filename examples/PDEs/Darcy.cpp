/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <Rodin/Variational/LinearElasticity.h>

#include <Rodin/Assembly/Sequential.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int argc, char** argv)
{
  // Load mesh
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 16, 16 });
  mesh.getConnectivity().compute(1, 2);

  // Functions
  P1 vh(mesh);
  P0 ph(mesh);

  TrialFunction u(vh);
  TestFunction  v(vh);

  TrialFunction p(ph);
  TestFunction  q(ph);

  Problem darcy(u, p, v, q);
  darcy = Integral(Grad(u), Grad(v))
        - Integral(v);
        //+ DirichletBC(u, Zero());
  darcy.assemble();

  Problem poisson(u, v);
  poisson = Integral(Grad(u), Grad(v))
          - Integral(v);
          //+ DirichletBC(u, Zero());
  poisson.assemble();

  std::cout << darcy.getStiffnessOperator().norm() << std::endl;
  std::cout << poisson.getStiffnessOperator().norm() << std::endl;

  // Assembly::Sequential<std::vector<Eigen::Triplet<Scalar>>, decltype(t)> assembly;

  // auto res = is.reduce([](const auto& l, const auto& r){return l +r ;});
  // std::cout << res << std::endl;


  // trw.apply([](auto&& v){ std::cout << v << " "; });

  return 0;
}


