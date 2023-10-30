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
  mesh = mesh.UniformGrid(Polytope::Type::Triangle, 6, 6);
  mesh.getConnectivity().compute(1, 2);
  mesh.getConnectivity().compute(2, 1);
  mesh.getConnectivity().compute(1, 1);

  mesh.save("miaow.mesh", IO::FileFormat::MEDIT);

  auto skin = mesh.skin();
  const auto& inc = skin.getConnectivity().getIncidence(1, 1);
  for (size_t i = 0; i < inc.size(); i++)
  {
    std::cout << i << ": ";
    for (const Index s : inc[i])
      std::cout << s << " ";
    std::cout << std::endl;
  }
  skin.save("skin.mesh", IO::FileFormat::MEDIT);
  return 0;
}

