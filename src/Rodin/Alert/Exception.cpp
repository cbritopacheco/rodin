/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#include <iostream>

#include <rang.hpp>

#include "Exception.h"

namespace Rodin::Alert
{
   Exception::Exception(const std::string& what)
      : Alert(what)
   {}

   void Exception::raise()
   {
      std::cerr << rang::fg::red
                << "Error: "
                << rang::fg::reset
                << what()
                << std::endl;
      std::abort();
   }
}

