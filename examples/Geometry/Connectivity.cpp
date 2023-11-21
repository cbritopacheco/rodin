/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <chrono>
#include <Rodin/Geometry.h>
#include <Rodin/Alert.h>

using namespace Rodin;
using namespace Geometry;


int main(int, char**)
{
  Mesh mesh;
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, 1024, 1024);

  auto t0 = std::chrono::high_resolution_clock::now();
  mesh.getConnectivity().compute(1, 2);
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 1);
  auto t1 = std::chrono::high_resolution_clock::now();
  auto d = t1 - t0;
  std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(d).count() << std::endl;
  return 0;
}

