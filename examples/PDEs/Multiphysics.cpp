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

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

template <class T>
struct IsFloatingPoint
{
  static constexpr Boolean Value = std::is_floating_point_v<std::decay_t<T>>;
};

template <class T>
auto foo(const T& v)
{
  return 2 * v;
}

int main(int argc, char** argv)
{
  // Load mesh
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, 16, 16);
  mesh.getConnectivity().compute(1, 2);

  // Functions
  P1 vh(mesh);

  TrialFunction u(vh);
  TestFunction  v(vh);

  Problem problem(u, v, u, v, v);

  // trw.apply([](auto&& v){ std::cout << v << " "; });
  std::cout << std::endl;


  return 0;
}

