/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Solver.h>
#include <Rodin/Geometry.h>
#include <Rodin/Variational.h>
#include <RodinExternal/MMG.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;
using namespace Rodin::External;

inline
Scalar K(const Point& x, const Point& y)
{
  return 1. / (4 * M_PI * ((x.asVector() - y.asVector()).stableNorm()));
}

int main(int, char**)
{
  size_t n = 32;
  MMG::Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, n, n);
  mesh.scale(2. / (n - 1));
  mesh.displace(VectorFunction{-1, -1});

  auto it = mesh.getElement();
  // for (auto it = mesh.getVertex(); it; ++it)
  // {
  //   const auto& coords = it->getCoordinates();
  //   const Scalar x = coords.x();
  //   const Scalar y = coords.y();
  //   const Scalar u = x * sqrt(x * x + y * y - x * x * y * y) / sqrt(x * x + y * y);
  //   const Scalar v = y * sqrt(x * x + y * y - x * x * y * y) / sqrt(x * x + y * y);
  //   mesh.setVertexCoordinates(it->getIndex(), u, 0);
  //   mesh.setVertexCoordinates(it->getIndex(), v, 1);
  // }
  MMG::Optimizer().setHMax(0.04).setHMin(0.002).optimize(mesh);

  for (const Index& i : mesh.getRidges())
  {
    mesh.setAttribute({ 1, i }, 1);
  }
  mesh.flush();

  mesh.save("D1.mesh");

  P1 fes(mesh);
  // GridFunction f(fes);
  ScalarFunction f =
    [](const Geometry::Point& p)
    {
      return 4.0 / (M_PI * std::sqrt(std::abs(1 - p.asVector().squaredNorm())));
    };

  Potential sl(K, f);
  GridFunction u(fes);
  u = sl;
  u.setWeights();
  u.save("u.gf");

  Alert::Info() << Integral(u).compute() << Alert::Raise;
  Alert::Info() << mesh.getVolume() << Alert::Raise;

  GridFunction fn(fes);
  fn = f;
  fn.save("f.gf");
  fn.setWeights();
  Alert::Info() << Integral(fn).compute() << Alert::Raise;


  return 0;
}

