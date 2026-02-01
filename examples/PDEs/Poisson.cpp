/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include "Rodin/Geometry/Polytope.h"
#include "Rodin/IO/ForwardDecls.h"
#include "Rodin/Solver/ForwardDecls.h"
#include "Rodin/Variational/TrialFunction.h"
#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Assembly.h>
#include <Rodin/Variational.h>
#include <type_traits>

using namespace Rodin;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  Mesh mesh;
  mesh = mesh.Build().initialize(2).nodes(3)
                     .vertex({0.0, 0.0})
                     .vertex({1.0, 0.0})
                     .vertex({0.0, 1.0})
                     .polytope(Polytope::Type::Triangle, {0, 1, 2})
                     .finalize();

  P1 vh(mesh, 2);

  TrialFunction u(vh);
  TestFunction v(vh);

  const Real lambda = 1.0, mu = 1.0;
  Problem prob(u, v);
      prob = Integral(lambda * Div(u), Div(v))
                 + Integral(
                     mu * (Jacobian(u) + Jacobian(u).T()), 0.5 * (Jacobian(v) + Jacobian(v).T()));
  prob.assemble();
  std::cout << "Operator:\n";
  std::cout << prob.getLinearSystem().getOperator() << std::endl;


  // P0 p0h(mesh);
  // P1 p1h(mesh);

  // TrialFunction u(p1h);
  // TestFunction v(p1h);

  // TrialFunction p(p0h);
  // TestFunction q(p0h);

  // RealFunction f = [](const Geometry::Point& p) -> double
  // {
  //   return 2.0 * (p.x() * (1.0 - p.x()) + p.y() * p.y());
  // };

  // Problem p_l2(p, q);
  // p_l2 = Integral(p, q) - Integral(f, q);
  // CG(p_l2).solve();

  // p.getSolution().save("p0.mesh");

  // Problem u_l2(u, v);
  // u_l2 = Integral(u, v) - Integral(p.getSolution(), v);
  // CG(u_l2).solve();

  // u.getSolution().save("u0.mesh");

  // p.getSolution() = 0.0;
  // u.getSolution() = 0.0;

  // Problem mixed(u, v, p, q);
  // mixed = Integral(u, v)
  //       - Integral(p, v)
  //       + Integral(p, q)
  //       - Integral(f, q);

  // BiCGSTAB(mixed).solve();

  // std::cout << "Operator:\n";
  // std::cout << mixed.getLinearSystem().getOperator() << std::endl;

  // std::cout << "Vector:\n";
  // std::cout << mixed.getLinearSystem().getVector() << std::endl;

  // u.getSolution().save("u.mesh");
  // p.getSolution().save("p.mesh");
  // mesh.save("mesh.mesh");

  // // mesh.load("CoronaryArtery.mesh", IO::FileFormat::MEDIT);
  // // mesh.getConnectivity().compute(2, 3);

  // // auto skin = mesh.skin();
  // // skin.save("CoronaryArtery_Skin.mesh", IO::FileFormat::MEDIT);

  return 0;
}
