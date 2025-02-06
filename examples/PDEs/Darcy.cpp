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


constexpr Real pi = Math::Constants::pi();

int main(int argc, char** argv)
{
  // Load mesh
  Mesh mesh;
  // mesh.load("Star.mesh");
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, { 32, 32 });
  mesh.scale(1. / (31));
  mesh.getConnectivity().compute(1, 2);

  // Functions
  P1 vh(mesh);
  P1 ph(mesh, 2);

  TrialFunction u(vh);
  TestFunction  w1(vh);

  TrialFunction v(ph);
  TestFunction  w2(ph);

  RealFunction f1 =
    [](const Geometry::Point& p)
    {
      return 2 * pi * pi * sin(pi * p.x()) * sin(pi * p.y());
    };

  VectorFunction f2 = VectorFunction{
    [](const Geometry::Point& p)
    {
      return 2 * pi * (
          cos(2 * pi * p.x()) * sin(pi * p.y()) * sin(pi * p.y())
          + cos(2 * pi * p.y()) * sin(pi * p.x()) * sin(pi * p.x()))
        - 4 * pi * pi * sin(pi * p.x()) * sin(pi * p.x()) * sin(pi * p.y()) * sin(pi * p.y());
    },
    [](const Geometry::Point& p)
    {
      return -2 * pi * (
          cos(2 * pi * p.x()) * sin(pi * p.y()) * sin(pi * p.y())
          + cos(2 * pi * p.y()) * sin(pi * p.x()) * sin(pi * p.x()))
        + 4 * pi * pi * sin(pi * p.x()) * sin(pi * p.x()) * sin(pi * p.y()) * sin(pi * p.y());
    }
  };

  // RealFunction f1 =
  //   [](const Geometry::Point& p)
  //   {
  //     return 2 * pi * pi * sin(pi * p.x()) * sin(pi * p.y())
  //       + 2 * pi * sin(pi * p.x()) * cos(pi * p.x()) * sin(pi * p.y()) * sin(pi * p.y())
  //       - 2 * pi * sin(pi * p.y()) * cos(pi * p.y()) * sin(pi * p.x()) * sin(pi * p.x());
  //   };

  // VectorFunction f2 = VectorFunction{
  //   [](const Geometry::Point& p)
  //   {
  //     return -2 * pi * (
  //         cos(2 * pi * p.x()) * sin(pi * p.y()) * sin(pi * p.y())
  //         + cos(2 * pi * p.y()) * sin(pi * p.x()) * sin(pi * p.x()))
  //       + 4 * pi * pi * sin(pi * p.x()) * sin(pi * p.x()) * sin(pi * p.y()) * sin(pi * p.y())
  //       + pi * cos(pi * p.x()) * sin(pi * p.y());
  //   },
  //   [](const Geometry::Point& p)
  //   {
  //     return 2 * pi * (
  //         cos(2 * pi * p.x()) * sin(pi * p.y()) * sin(pi * p.y())
  //         + cos(2 * pi * p.y()) * sin(pi * p.x()) * sin(pi * p.x()))
  //       + 4 * pi * pi * sin(pi * p.x()) * sin(pi * p.x()) * sin(pi * p.y()) * sin(pi * p.y())
  //       + pi * sin(pi * p.x()) * cos(pi * p.y());
  //   }
  // };

  auto n = BoundaryNormal(mesh);

  Problem darcy(u, v, w1, w2);
  darcy = Integral(Grad(u), Grad(w1))
        // - Integral(Div(v), w1)
        - Integral(f1, w1)
        + Integral(Jacobian(v), Jacobian(w2))
        // + Integral(u, Div(w2))
        - Integral(f2, w2)
        + DirichletBC(u, Zero())
        + DirichletBC(v, Zero(2))
        ;

  UMFPack(darcy).solve();

  mesh.save("th.mesh");
  u.getSolution().save("u.gf");
  v.getSolution().save("v.gf");

  GridFunction uEx(vh);
  uEx = [](const Geometry::Point& p) { return sin(M_PI * p.x()) * sin(M_PI * p.y()); };
  uEx.save("uEx.gf");

  GridFunction vEx(ph);
  vEx = VectorFunction{
    [](const Geometry::Point& p)
    {
      return sin(M_PI * p.x()) * sin(M_PI * p.y());
    },
    [](const Geometry::Point& p)
    {
      return sin(M_PI * p.x()) * sin(M_PI * p.y());
    }
  };
  vEx.save("vEx.gf");


  VectorFunction f3 = VectorFunction{
    [](const Geometry::Point& p)
    {
      return 2 * pi * pi * sin(pi * p.x()) * sin(pi * p.y());
    },
    [](const Geometry::Point& p)
    {
      return 2 * pi * pi * sin(pi * p.x()) * sin(pi * p.y());
    }
  };

  TrialFunction v2(ph);
  Problem darcy2(v2, w2);
  darcy2 = Integral(Jacobian(v2), Jacobian(w2))
         - Integral(f3, w2)
         + DirichletBC(v2, Zero(2))
         ;


  return 0;
}


