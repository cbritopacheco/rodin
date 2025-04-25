/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <chrono>
#include <Rodin/Types.h>
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>

using namespace Rodin;
using namespace Rodin::Solver;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  Mesh mesh;
  mesh.load("Sphere.mesh", IO::FileFormat::MEDIT);

  mesh.scale(1.125 / 1.25);
  mesh.getConnectivity().compute(2, 3);
  auto skin = mesh.skin();

  auto keep = skin.keep({ 2, 3 });

  keep.save("Sphere1_125.mesh", IO::FileFormat::MEDIT);


  // // Build a mesh
  // Mesh mesh;
  // //mesh = mesh.UniformGrid(Polytope::Type::Quadrilateral, { 2048, 2048 });
  // // mesh.save("big.mesh");
  // mesh.load("big.mesh");

  // ScalarFunction f = [](const Point& p)
  // {
  //   return p.x() + p.y();
  // };



  // // mesh.scale(1.0 / 31);
  // // mesh.getConnectivity().compute(1, 2);

  // // Functions
  // P1 vh(mesh);
  // GridFunction gf(vh);

  // using std::chrono::high_resolution_clock;
  // using std::chrono::duration_cast;
  // using std::chrono::duration;
  // using std::chrono::milliseconds;

  // auto t1 = high_resolution_clock::now();
  // gf = f;
  // auto t2 = high_resolution_clock::now();

  // /* Getting number of milliseconds as an integer. */
  // auto ms_int = duration_cast<milliseconds>(t2 - t1);
/*//  Getting number of milliseconds as a double. */
  //   duration<double, std::milli> ms_double = t2 - t1;

  //   std::cout << ms_int.count() << "ms\n";
  //   std::cout << ms_double.count() << "ms\n";

  // gf.save("big.gf");

  // // ScalarFunction f(1.0);

  // // TrialFunction u(vh);
  // // TestFunction  v(vh);

  // // // Define problem
  // // Problem poisson(u, v);
  // // poisson = Integral(Grad(u), Grad(v))
  // //         - Integral(f, v)
  // //         + DirichletBC(u, Zero());

  // // // Solve
  // // CG(poisson).solve();

  // // // Save solution
  // // u.getSolution().save("Poisson.gf");
  // // mesh.save("Poisson.mesh");

  return 0;
}
