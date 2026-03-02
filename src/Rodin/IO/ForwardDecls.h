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
  /**
   * @brief Base template for loading objects from files or streams.
   * @tparam T Type of object to load
   * @see Loader
   */
  template <class T>
  class Loader;

  /**
   * @brief Base template for printing objects to streams.
   * @tparam T Type of object to print
   * @see Printer
   */
  template <class T>
  class Printer;

  /**
   * @brief Enumeration of supported file formats for mesh and grid function I/O.
   *
   * This enum identifies the different file formats that can be used for
   * reading and writing meshes and grid functions in Rodin.
   */
  enum class FileFormat
  {
    MFEM,     ///< MFEM mesh format - native format for MFEM library
    MEDIT,    ///< MEDIT mesh format - used by MMG remeshing software
  };

  /**
   * @brief Loader template for meshes with specific file format and context.
   * @tparam Fmt File format to use for loading
   * @tparam Trait Context trait (e.g., sequential or parallel)
   * @see MeshLoader
   */
  template <FileFormat Fmt, class Trait>
  class MeshLoader;

  /**
   * @brief Printer template for meshes with specific file format and context.
   * @tparam Fmt File format to use for printing
   * @tparam Trait Context trait (e.g., sequential or parallel)
   * @see MeshPrinter
   */
  template <FileFormat Fmt, class Trait>
  class MeshPrinter;

  /**
   * @brief Loader template for grid functions with specific file format.
   * @tparam Fmt File format to use for loading
   * @tparam FES Finite element space type
   * @tparam Data Data storage type
   * @see GridFunctionLoader
   */
  template <FileFormat Fmt, class FES, class Data>
  class GridFunctionLoader;

  /**
   * @brief Printer template for grid functions with specific file format.
   * @tparam Fmt File format to use for printing
   * @tparam FES Finite element space type
   * @tparam Data Data storage type
   * @see GridFunctionPrinter
   */
  template <FileFormat Fmt, class FES, class Data>
  class GridFunctionPrinter;

  /**
   * @brief Converts a FileFormat enum value to its string representation.
   * @param[in] fmt File format enum value
   * @returns C-style string name of the format, or nullptr if invalid
   *
   * ## Example
   * ```cpp
   * const char* name = toCharString(FileFormat::MFEM);  // Returns "MFEM"
   * ```
   */
  inline
  constexpr
  const char* toCharString(FileFormat fmt)
  {
    switch (fmt)
    {
      case FileFormat::MFEM:
        return "MFEM";
      case FileFormat::MEDIT:
        return "MEDIT";
    }
    return nullptr;
  }

  /**
   * @brief Stream output operator for FileFormat enum.
   * @param[in,out] os Output stream
   * @param[in] fmt File format to output
   * @returns Reference to the output stream
   *
   * Writes the string representation of the file format to the stream.
   *
   * ## Example
   * ```cpp
   * std::cout << FileFormat::MEDIT;  // Outputs "MEDIT"
   * ```
   */
  inline
  std::ostream& operator<<(std::ostream& os, FileFormat fmt)
  {
    os << toCharString(fmt);
    return os;
  }
}

#endif

