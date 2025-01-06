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
#include <Rodin/Solver/UMFPack.h>

#include <Rodin/Assembly/Sequential.h>

using namespace Rodin;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int argc, char** argv)
{
  // Load mesh
  Mesh mesh;
  mesh.load("Star.mesh");
  // mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 16, 16 });
  // mesh.scale(1. / (15));
  // mesh.displace(VectorFunction{-1, -1});
  mesh.getConnectivity().compute(1, 2);

  // Functions
  P1 vh(mesh, 2);
  P1 ph(mesh);

  TrialFunction u(vh);
  TestFunction  v(vh);

  TrialFunction p(ph);
  TestFunction  q(ph);

  RealFunction f = [](const Geometry::Point& p) { return exp(p.x()) * sin(p.y()); };
  VectorFunction g = {
    [](const Geometry::Point& p) { return exp(p.x()) * cos(p.y()); },
    [](const Geometry::Point& p) { return exp(p.x()) * sin(p.y()); } };

  auto n = BoundaryNormal(mesh);

  Problem darcy(u, p, v, q);
  darcy = //Integral(u, v)
        - Integral(p, Div(v))
        // + BoundaryIntegral(f, Dot(v, n))
        - Integral(Div(u), q);


  // darcy = Integral(Jacobian(u), Jacobian(v))
  //       - Integral(VectorFunction{1, 0}, v)
  //       + Integral(Grad(p), Grad(q))
  //       - Integral(Div(u), q)
  //       + DirichletBC(u, VectorFunction{2, 0})
  //       + DirichletBC(p, Zero());

  darcy.assemble();

  //std::cout << darcy.getStiffnessOperator().block(0, 0, vh.getSize(), vh.getSize()) << std::endl;
  // std::cout << darcy.getStiffnessOperator().block(0, vh.getSize(), vh.getSize(), ph.getSize()) << std::endl;
  // std::cout << darcy.getStiffnessOperator().block(vh.getSize(), 0, ph.getSize(), vh.getSize()) << std::endl;
  // auto a = darcy.getStiffnessOperator().block(vh.getSize(), 0, ph.getSize(), vh.getSize()).toDense();
  // auto b = darcy.getStiffnessOperator().block(vh.getSize(), 0, ph.getSize(), vh.getSize()).toDense();
  // std::cout << (a - b).norm() << std::endl;
  std::exit(1);

  std::cout << "miaow\n";

  // UMFPack(darcy).solve();

  u.getSolution().save("u.gf");
  mesh.save("u.mesh");

  p.getSolution().save("p.gf");
  mesh.save("q.mesh");

  return 0;
}


