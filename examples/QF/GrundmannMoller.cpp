/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <Rodin/QF/GrundmannMoller.h>
#include <Rodin/QF/QF1P1.h>

using namespace Rodin;
using namespace Rodin::Geometry;

int main(int, char**)
{
  QF::GrundmannMoller qfe(1, Polytope::Type::Segment);
  std::cout << "Order: " << qfe.getOrder() << std::endl;
  std::cout << "Size: " << qfe.getSize() << std::endl;
  for (size_t i = 0; i < qfe.getSize(); i++)
  {
    std::cout << "Weight:\n";
    std::cout << qfe.getWeight(i) << std::endl;
    std::cout << "Point:\n" << qfe.getPoint(i) << std::endl;
  }

  // QF::GrundmannMoller qf(1, Polytope::Type::Triangle);
  // std::cout << "Order: " << qf.getOrder() << std::endl;
  // std::cout << "Size: " << qf.getSize() << std::endl;
  // for (size_t i = 0; i < qf.getSize(); i++)
  // {
  //   std::cout << "Weight:\n";
  //   std::cout << qf.getWeight(i) << std::endl;
  //   std::cout << "Point:\n" << qf.getPoint(i) << std::endl;
  // }

  return 0;
}
