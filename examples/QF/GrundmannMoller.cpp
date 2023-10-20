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
  QF::GrundmannMoller qf(0, Polytope::Type::Triangle);
  std::cout << qf.getOrder() << std::endl;
  for (size_t i = 0; i < qf.getSize(); i++)
  {
    std::cout << "Weight:\n";
    std::cout << qf.getWeight(i) << std::endl;
    std::cout << "Point:\n" << qf.getPoint(i) << std::endl;
  }

  std::cout << "miaow ?\n";
  QF::QF1P1 qf1p1(Polytope::Type::Triangle);
  for (size_t i = 0; i < qf.getSize(); i++)
  {
    std::cout << "Weight:\n";
    std::cout << qf1p1.getWeight(i) << std::endl;
    std::cout << "Point:\n" << qf1p1.getPoint(i) << std::endl;
  }
  return 0;
}
