/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_LOADER_H
#define RODIN_IO_LOADER_H

#include <istream>

#include "ForwardDecls.h"

namespace Rodin::IO
{
   template <class T>
   class Loader
   {
      public:
         virtual IO::Status load(std::istream& is) = 0;
   };
}

#endif
