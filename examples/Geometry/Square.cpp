/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Geometry.h>
#include <Rodin/Alert.h>

using namespace Rodin;
using namespace Geometry;

int main(int, char**)
{
  Mesh mesh;
  mesh.load("../resources/mfem/StarSquare.mfem.mesh", IO::FileFormat::MFEM);
  // mesh = mesh.Build().initialize(2).nodes(4)
  //                    .vertex({ 0, 0 }).vertex({ 1, 0 }).vertex({ 0, 1 }).vertex({ 1, 1 })
  //                    .polytope(Polytope::Type::Quadrilateral, { 0, 1, 2, 3 }).finalize();

  mesh.getConnectivity().compute(1, 2);
  // const auto& inc = mesh.getConnectivity().getIncidence(1, 0);
  // for (auto& c : inc)
  // {
  //   for (auto& v : c)
  //     std::cout << v << ",";
  //   std::cout << std::endl;
  // }

  mesh.save("miaow.mesh", IO::FileFormat::MFEM);

  return 0;
}


