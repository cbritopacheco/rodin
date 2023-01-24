/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_FORWARDDECLS_H
#define RODIN_IO_FORWARDDECLS_H

#include <string>
#include <ostream>
#include <optional>
#include <boost/filesystem.hpp>

namespace Rodin::IO
{
  template <class T>
  class Loader;

  template <class T>
  class Printer;

  /**
   * @brief Enumeration class for different file formats.
   */
  enum class FileFormat
  {
    MFEM, ///< MFEM file format
    GMSH, ///< GMSH file format
    MEDIT ///< MEDIT file format
  };

  template <FileFormat fmt, class Trait>
  class MeshLoader;

  template <FileFormat fmt, class Trait>
  class MeshPrinter;

  template <FileFormat fmt, class FES>
  class GridFunctionLoader;

  template <FileFormat fmt, class FES>
  class GridFunctionPrinter;

  std::ostream& operator<<(std::ostream& os, FileFormat fmt);
}

#endif

