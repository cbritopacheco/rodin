/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <RodinExternal/MMG.h>
#include <Rodin/Variational.h>

#include <Rodin/Geometry/RippleSimplexification.h>

using namespace Rodin;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::External;

int main(int, char**)
{
  Mesh mesh;
  mesh = mesh.Build().initialize(3)
                     .nodes(6)
                     .vertex({0, 0, 0})
                     .vertex({1, 0, 0})
                     .vertex({0, 1, 0})

                     .vertex({0, 0, 1})
                     .vertex({1, 0, 1})
                     .vertex({0, 1, 1})

                     // .vertex({0, 0, 2})
                     // .vertex({1, 0, 2})
                     // .vertex({0, 1, 2})

                     .polytope(Polytope::Type::TriangularPrism, { 0, 1, 2, 3, 4, 5 })
                     // .polytope(Polytope::Type::TriangularPrism, { 3, 4, 5, 6, 7, 8 })
                     .finalize();
  mesh.getConnectivity().compute(2, 3);

  RippleSimplexification(mesh).simplexify();
  std::exit(1);

  mesh.save("prism.mesh");


  Mesh tet;
  auto build = tet.Build();
  build.initialize(3).nodes(mesh.getVertexCount());
  build.setVertices(mesh.getVertices());
  build.tetrahedron(IndexArray{{ 1, 2, 5, 3 }});
  build.tetrahedron(IndexArray{{ 2, 3, 0, 1 }});
  build.tetrahedron(IndexArray{{ 3, 4, 1, 5  }});
  // build.tetrahedron(IndexArray{ 0, 2, 5, 3, 1 });
  tet = build.finalize();
  tet.save("tet.mesh", IO::FileFormat::MEDIT);



  std::exit(1);

  const auto& conn = mesh.getConnectivity();
  std::vector<Connectivity<Context::Local>::SubPolytope> subpolytopes;
  for (auto it = mesh.getCell(); it; ++it)
  {
    conn.getSubPolytopes(subpolytopes, it->getIndex(), 2);
    const auto& vs = it->getVertices();
    for (auto& [geometry, svs] : subpolytopes)
    {
      // Reorient quad
      Index idx;
      svs.minCoeff(&idx);
      if (geometry == Polytope::Type::Quadrilateral)
      {
        if (idx == 0)
        {
        }
        else if (idx == 1)
        {
          std::swap(svs(0), svs(1));
          std::swap(svs(2), svs(3));
        }
        else if (idx == 2)
        {
          std::swap(svs(0), svs(2));
          std::swap(svs(1), svs(3));
        }
        else if (idx == 3)
        {
          std::swap(svs(0), svs(3));
          std::swap(svs(1), svs(2));
        }
        else
        {
          std::cout << "miaow\n";
        }
      }
    }

    const auto& f1 = subpolytopes[1].vertices;
    const auto& f2 = subpolytopes[2].vertices;
    const auto& f3 = subpolytopes[3].vertices;

    build.polytope(Polytope::Type::Tetrahedron, { vs(0), vs(1), vs(2), vs(3) });
    build.polytope(Polytope::Type::Tetrahedron, { f2(0), f2(1), f2(2), f1(0) });
    build.polytope(Polytope::Type::Tetrahedron, { f2(1), f3(2), f3(3), f3(3) });
  }
  tet = build.finalize();
  tet.getConnectivity().compute(2, 3);
  tet.save("tet.mesh", IO::FileFormat::MEDIT);
  mesh.getConnectivity().compute(2, 3);

  auto skin2 = mesh.skin();
  skin2.save("skin.mesh", IO::FileFormat::MEDIT);
  mesh.save("prism.mesh");
  std::exit(1);




  for (auto it = mesh.getFace(); it; ++it)
  {
    mesh.setAttribute({2, it->getIndex()}, it->getIndex() + 1);
  }
  std::cout << mesh.getConnectivity().getIncidence({3, 2}, 0).size() << std::endl;
  mesh.save("prism.mesh");
  mesh.load("prism.mesh");
  mesh.save("prism2.mesh");
  std::exit(1);

  mesh.load("atria_fluid.mesh", IO::FileFormat::MEDIT);
  mesh.getConnectivity().compute(2, 3);
  auto skin = mesh.skin();
  skin.save("atria_fluid_skin.mesh", IO::FileFormat::MEDIT);
  std::exit(1);

  P1 vh(mesh);
  RealFunction center = 0.01;
  RealFunction one = 1e10;

  TrialFunction u(vh);
  TestFunction  v(vh);
  Problem poisson(u, v);
  poisson = 0.00001 * Integral(Grad(u), Grad(v))
          + Integral(one, v)
          + DirichletBC(u, center);
  CG(poisson).solve();

  // u.getSolution().save("atria_fluid2.o.sol", IO::FileFormat::MEDIT);

  // mesh.save("atria_fluid3.mesh");
  // std::exit(1);

  // P1 vh(mesh);
  // GridFunction gf(vh);
  // // gf.load("atria_fluid.o.sol");

  // Math::Vector<Real> c{{43, 1340, -3}};
  // gf = [&](const Point& p)
  // {
  //   auto r = (p - c).norm();
  //   return r;
  // };
  Real max = u.getSolution().max(), min = u.getSolution().min();
  Real hmax = 1, hmin = -1;
  GridFunction newgf(vh);
  newgf = hmin + (u.getSolution() - min) * (hmax - hmin) / (max - min);
  newgf.save("atria_fluid.o.sol", IO::FileFormat::MEDIT);

  // mesh.getConnectivity().compute(1, 2);
  // // std::set<size_t> attrs;
  // // size_t count54 = 0;
  // // size_t count120 = 0;
  // // for (auto edge = mesh.getFace(); edge; ++edge)
  // // {
  // //   if (edge->getAttribute() == 54)
  // //     count54++;
  // //   else if (edge->getAttribute() == 120)
  // //     count120++;
  // // }
  // // std::cout << count54 << " " << count120 << std::endl;

  // mesh.trace({{6, 666}});
  // mesh.trace({{{6, 7}, 666}});

  // mesh.save("test.mesh", IO::FileFormat::MEDIT);

  // // Build a mesh
  // Mesh mesh;
  // mesh = mesh.UniformGrid(Polytope::Type::Quadrilateral, { 32, 32 });
  // mesh.scale(1.0 / 31);
  // mesh.getConnectivity().compute(1, 2);

  // // Functions
  // P1 vh(mesh);

  // // auto f = cos(2 * M_PI * F::x) * sin(2 * M_PI * F::y);
  // ScalarFunction f(1.0);

  // TrialFunction u(vh);
  // TestFunction  v(vh);

  // // Define problem
  // Problem poisson(u, v);
  // poisson = Integral(Grad(u), Grad(v))
  //         - Integral(f, v)
  //         + DirichletBC(u, Zero());

  // // Solve
  // CG(poisson).solve();

  // // Save solution
  // u.getSolution().save("Poisson.gf");
  // mesh.save("Poisson.mesh");

  return 0;
}
