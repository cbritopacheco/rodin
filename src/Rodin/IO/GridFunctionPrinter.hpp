/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRIDFUNCTIONPRINTER_HPP
#define RODIN_VARIATIONAL_GRIDFUNCTIONPRINTER_HPP

#include "GridFunctionPrinter.h"

#include "MEDIT.h"

namespace Rodin::IO
{
   template <class FES>
   void GridFunctionPrinter<FileFormat::MEDIT, FES>::print(std::ostream& os)
   {
      assert(false);
   }
}

#endif
