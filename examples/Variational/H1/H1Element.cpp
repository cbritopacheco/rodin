/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/Variational/H1/H1Element.h>

using namespace Rodin;
using namespace Rodin::Geometry;
using namespace Rodin::Variational;

int main(int, char**)
{
  H1Element<2, Real> e1(Polytope::Type::Triangle);

  std::cout << "Number of DOFs: " << e1.getCount() << std::endl;
  std::cout << "Nodes:" << std::endl;

  for (size_t i = 0; i < e1.getCount(); ++i)
  {
    std::cout << e1.getNode(i) << std::endl;
    std::cout << "----------" << std::endl;
  }
}



