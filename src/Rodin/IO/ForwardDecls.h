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
   * @see <a href="class_rodin_1_1_i_o_1_1_loader.html">Loader</a>
   */
  template <class T>
  class Loader;

  /**
   * @brief Base template for printing objects to streams.
   * @tparam T Type of object to print
   * @see <a href="class_rodin_1_1_i_o_1_1_printer.html">Printer</a>
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
    HDF5,     ///< HDF5 hierarchical binary format for raw mesh/function storage
  };

  /**
   * @brief Loader template for meshes with specific file format and context.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | [MeshLoader<IO::FileFormat::MFEM, Context::Local>](class_rodin_1_1_i_o_1_1_mesh_loader_3_01_i_o_1_1_file_format_1_1_m_f_e_m_00_01_context_1_1_local_01_4.html) | Loads local meshes from MFEM mesh files. |
   * | [MeshLoader<IO::FileFormat::MEDIT, Context::Local>](class_rodin_1_1_i_o_1_1_mesh_loader_3_01_i_o_1_1_file_format_1_1_m_e_d_i_t_00_01_context_1_1_local_01_4.html) | Loads local meshes from MEDIT `.mesh` files. |
   * | [MeshLoader<IO::FileFormat::HDF5, Context::Local>](class_rodin_1_1_i_o_1_1_mesh_loader_3_01_file_format_1_1_h_d_f5_00_01_context_1_1_local_01_4.html) | Loads complete local meshes from Rodin HDF5 datasets. |
   * | [MeshLoader<IO::FileFormat::HDF5, Context::MPI>](class_rodin_1_1_i_o_1_1_mesh_loader_3_01_file_format_1_1_h_d_f5_00_01_context_1_1_m_p_i_01_4.html) | Loads distributed MPI meshes from rank-local Rodin HDF5 @ref Rodin::Geometry::Shard "shards". |
   *
   * @tparam Fmt File format to use for loading
   * @tparam Trait Context trait (e.g., sequential or parallel)
   * @see <a href="class_rodin_1_1_i_o_1_1_mesh_loader.html">MeshLoader</a>
   */
  template <FileFormat Fmt, class Trait>
  class MeshLoader;

  /**
   * @brief Printer template for meshes with specific file format and context.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | [MeshPrinter<IO::FileFormat::MFEM, Context::Local>](class_rodin_1_1_i_o_1_1_mesh_printer_3_01_file_format_1_1_m_f_e_m_00_01_context_1_1_local_01_4.html) | Writes local meshes in MFEM mesh format. |
   * | [MeshPrinter<IO::FileFormat::MEDIT, Context::Local>](class_rodin_1_1_i_o_1_1_mesh_printer_3_01_file_format_1_1_m_e_d_i_t_00_01_context_1_1_local_01_4.html) | Writes local meshes in MEDIT `.mesh` format. |
   * | [MeshPrinter<IO::FileFormat::HDF5, Context::Local>](class_rodin_1_1_i_o_1_1_mesh_printer_3_01_file_format_1_1_h_d_f5_00_01_context_1_1_local_01_4.html) | Writes complete local meshes to Rodin HDF5 datasets. |
   * | [MeshPrinter<IO::FileFormat::HDF5, Context::MPI>](class_rodin_1_1_i_o_1_1_mesh_printer_3_01_file_format_1_1_h_d_f5_00_01_context_1_1_m_p_i_01_4.html) | Writes distributed MPI @ref Rodin::Geometry::Shard "mesh shards" to rank-local Rodin HDF5 datasets. |
   *
   * @tparam Fmt File format to use for printing
   * @tparam Trait Context trait (e.g., sequential or parallel)
   * @see <a href="class_rodin_1_1_i_o_1_1_mesh_printer.html">MeshPrinter</a>
   */
  template <FileFormat Fmt, class Trait>
  class MeshPrinter;

  /**
   * @brief Loader template for grid functions with specific file format.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref GridFunctionLoader "GridFunctionLoader<IO::FileFormat::MFEM, FES, Data>" | Loads MFEM grid functions into Rodin-owned vector storage. |
   * | @ref GridFunctionLoader "GridFunctionLoader<IO::FileFormat::MFEM, FES, Math::Vector<Scalar>&>" | Loads MFEM grid functions into referenced vector storage. |
   * | @ref GridFunctionLoader "GridFunctionLoader<IO::FileFormat::MEDIT, FES, Data>" | Loads MEDIT grid functions into Rodin-owned vector storage. |
   * | @ref GridFunctionLoader "GridFunctionLoader<IO::FileFormat::MEDIT, FES, Math::Vector<Scalar>&>" | Loads MEDIT grid functions into referenced vector storage. |
   * | @ref GridFunctionLoader "GridFunctionLoader<IO::FileFormat::HDF5, FES, Math::Vector<Scalar>>" | Loads HDF5 grid functions into Rodin vector storage. |
   * | @ref GridFunctionLoader "GridFunctionLoader<IO::FileFormat::HDF5, FES, Vec>" | Loads HDF5 grid functions into PETSc vector storage. |
   *
   * @tparam Fmt File format to use for loading
   * @tparam FES Finite element space type
   * @tparam Data Data storage type
   * @see <a href="class_rodin_1_1_i_o_1_1_grid_function_loader.html">GridFunctionLoader</a>
   */
  template <FileFormat Fmt, class FES, class Data>
  class GridFunctionLoader;

  /**
   * @brief Printer template for grid functions with specific file format.
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | @ref GridFunctionPrinter "GridFunctionPrinter<IO::FileFormat::MFEM, FES, Math::Vector<Scalar>>" | Writes Rodin vector-backed grid functions in MFEM format. |
   * | @ref GridFunctionPrinter "GridFunctionPrinter<IO::FileFormat::MFEM, FES, Vec>" | Writes PETSc vector-backed grid functions in MFEM format. |
   * | @ref GridFunctionPrinter "GridFunctionPrinter<IO::FileFormat::MEDIT, FES, Math::Vector<Scalar>>" | Writes Rodin vector-backed grid functions in MEDIT format. |
   * | @ref GridFunctionPrinter "GridFunctionPrinter<IO::FileFormat::MEDIT, FES, Vec>" | Writes PETSc vector-backed grid functions in MEDIT format. |
   * | @ref GridFunctionPrinter "GridFunctionPrinter<IO::FileFormat::HDF5, FES, Math::Vector<Scalar>>" | Writes Rodin vector-backed grid functions to HDF5 datasets. |
   * | @ref GridFunctionPrinter "GridFunctionPrinter<IO::FileFormat::HDF5, FES, Vec>" | Writes PETSc vector-backed grid functions to HDF5 datasets. |
   *
   * @tparam Fmt File format to use for printing
   * @tparam FES Finite element space type
   * @tparam Data Data storage type
   * @see <a href="class_rodin_1_1_i_o_1_1_grid_function_printer.html">GridFunctionPrinter</a>
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
      case FileFormat::HDF5:
        return "HDF5";
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
