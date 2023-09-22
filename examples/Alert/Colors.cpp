/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <iostream>
#include <Rodin/Alert.h>

using namespace Rodin::Alert;

int main(int, char**)
{
  std::stringstream ss;
  ss << Color<RGB<255, 0, 0>>() << "miaow" << Reset << " quack";
}
