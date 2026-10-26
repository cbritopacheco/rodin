/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/// @brief Executable for the transient active atrial solid component.
#include <exception>
#include <iostream>
#include "Solid.h"

int main(int argc, char** argv)
{
  try
  {
    namespace AtriumSolid = Rodin::Examples::Heart::AtriumContraction;
    AtriumSolid::Solid atrium(AtriumSolid::readOptions(argc, argv));
    atrium.run();
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << std::endl;
    return 1;
  }
  return 0;
}
