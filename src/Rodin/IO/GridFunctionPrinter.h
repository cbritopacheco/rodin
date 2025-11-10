/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_GRIDFUNCTIONPRINTER_H
#define RODIN_VARIATIONAL_GRIDFUNCTIONPRINTER_H

#include <boost/filesystem.hpp>

#include "ForwardDecls.h"

namespace Rodin::IO
{
  /**
   * @brief Forward declaration of base class for printing grid functions.
   *
   * GridFunctionPrinterBase provides the foundation for writing finite element
   * solution data in different file formats. Specialized implementations exist
   * for each supported file format.
   *
   * @tparam Fmt File format to use for printing
   * @tparam FES Finite element space type
   * @tparam Data Data storage type (typically a vector type)
   *
   * @see GridFunctionPrinter
   */
  template <FileFormat Fmt, class FES, class Data>
  class GridFunctionPrinterBase;
}

#endif
