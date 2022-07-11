/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_PRINTER_H
#define RODIN_IO_PRINTER_H

#include <ostream>

#include "ForwardDecls.h"

namespace Rodin::IO
{
   template <class T>
   class Printer
   {
      public:
         virtual void print(std::ostream& os) = 0;

      protected:
         virtual const T& getObject() const = 0;
   };
}

#endif

