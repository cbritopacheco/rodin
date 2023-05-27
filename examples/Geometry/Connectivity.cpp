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
  // Mesh<Context::Serial>::Builder build;
  // const size_t dim = 2;

  // Mesh mesh;
  // mesh.initialize(dim, dim)
  //     .vertex({0, 0})
  //     .vertex({1, 0})
  //     .vertex({0, 1})
  //     .vertex({1, 1})
  //     .polytope(Polytope::Geometry::Triangle, {0, 1, 2})
  //     .polytope(Polytope::Geometry::Triangle, {1, 2, 3})
  //     .finalize();

  // mesh.setAttribute(dim, 0, 666);
  // mesh.setAttribute(dim, 1, 999);

  // std::cout << mesh.getConnectivity().getCount(1) << std::endl;

  // size_t d = 1, dp = 2;
  // mesh.getConnectivity().compute(d, dp);
  // for (const auto& s : mesh.getConnectivity().getIncidence(d, dp))
  // {
  //   std::cout << "\n-----\n";
  //   for (Index i : s)
  //     std::cout << i << ", ";
  //   std::cout << "\n-----\n";
  // }

  // // mesh.getConnectivity().compute(1, 0);
  // // mesh.getConnectivity().intersection(mesh.getDimension(), mesh.getDimension(), 0);

  // // mesh.getConnectivity().transpose(0, mesh.getDimension());
  // // std::cout << "count: " << mesh.getConnectivity().count(1) << std::endl;
  // // mesh.getConnectivity().allocate();
  // // const auto& s = mesh.getConnectivity().candidate(1, 1);
  // // for (auto& v : s)
  // // {
  // //   std::cout << "\n---\n";
  // //   for (auto& i : v)
  // //     std::cout << i << ", ";
  // // }
  // // std::cout << "\n-----\n";

  // mesh.save("square.mesh", IO::FileFormat::MEDIT);


  Mesh miaow;
  miaow.load("../resources/mfem/surface-meshes-example.mesh", IO::FileFormat::MFEM);
  // miaow.load(std::string(RODIN_RESOURCES_DIR) + "/mmg/Box.medit.mesh", IO::FileFormat::MEDIT);
  miaow.save("test.mesh", IO::FileFormat::MEDIT);

  return 0;
}

