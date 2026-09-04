/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_IO_PRINTER_H
#define RODIN_IO_PRINTER_H

#include <ostream>
#include <fstream>

#include <boost/filesystem/path.hpp>

#include "ForwardDecls.h"

namespace Rodin::IO
{
  /**
   * @defgroup PrinterSpecializations Printer Template Specializations
   * @brief Template specializations of the Printer class.
   * @see <a href="class_rodin_1_1_i_o_1_1_printer.html">Printer</a>
   *
   * | Specialization | Description |
   * |----------------|-------------|
   * | <a href="class_rodin_1_1_i_o_1_1_mesh_printer_3_01_file_format_1_1_m_f_e_m_00_01_context_1_1_local_01_4.html">MeshPrinter<IO::FileFormat::MFEM, Context::Local></a> | Writes sequential meshes in MFEM format. |
   * | <a href="class_rodin_1_1_i_o_1_1_mesh_printer_3_01_file_format_1_1_m_e_d_i_t_00_01_context_1_1_local_01_4.html">MeshPrinter<IO::FileFormat::MEDIT, Context::Local></a> | Writes sequential meshes in MEDIT format. |
   * | <a href="class_rodin_1_1_i_o_1_1_mesh_printer_3_01_file_format_1_1_h_d_f5_00_01_context_1_1_local_01_4.html">MeshPrinter<IO::FileFormat::HDF5, Context::Local></a> | Writes sequential meshes in HDF5 format. |
   * | <a href="class_rodin_1_1_i_o_1_1_mesh_printer_3_01_file_format_1_1_h_d_f5_00_01_context_1_1_m_p_i_01_4.html">MeshPrinter<IO::FileFormat::HDF5, Context::MPI></a> | Writes distributed mesh shards and metadata to per-rank HDF5 files. |
   */

  /**
   * @brief Abstract base class template for printing objects to streams.
   *
   * The Printer class provides a common interface for writing objects of type
   * @p T to output streams. This class serves as the foundation for specialized
   * printers that handle different file formats and object types.
   *
   * @tparam T Type of object to print
   *
   * Derived classes must implement:
   * - print(std::ostream&) - Write object to an output stream
   * - getObject() - Access the object being printed
   *
   * ## Usage Example
   * ```cpp
   * class MyObjectPrinter : public Printer<MyObject>
   * {
   * public:
   *   MyObjectPrinter(const MyObject& obj) : m_object(obj) {}
   *
   *   void print(std::ostream& os) override
   *   {
   *     // Write object to stream
   *     os << getObject();
   *   }
   *
   *   const ObjectType& getObject() const override
   *   {
   *     return m_object;
   *   }
   *
   * private:
   *   const MyObject& m_object;
   * };
   *
   * // Usage
   * MyObject obj;
   * MyObjectPrinter printer(obj);
   * printer.print(std::cout);
   * ```
   *
   * @see <a href="class_rodin_1_1_i_o_1_1_mesh_printer.html">MeshPrinter</a>
   * @see <a href="class_rodin_1_1_i_o_1_1_grid_function_printer.html">GridFunctionPrinter</a>
   */
  template <class T>
  class Printer
  {
    public:
      /**
       * @brief Type of object being printed.
       */
      using ObjectType = T;

      /**
       * @brief Prints object to an output stream.
       * @param[in,out] os Output stream to write to
       *
       * This pure virtual method must be implemented by derived classes to
       * write the object data to the provided stream.
       */
      virtual void print(std::ostream& os) = 0;

      /**
       * @brief Prints object to a file path.
       * @param[in] filename Output file path.
       */
      virtual void print(const boost::filesystem::path& filename)
      {
        std::ofstream os(filename.string());
        this->print(os);
      }

      /**
       * @brief Gets a const reference to the object being printed.
       * @returns Const reference to the object
       *
       * This method provides read-only access to the object being written.
       * Derived classes must implement this to return their internal object reference.
       */
      virtual const ObjectType& getObject() const = 0;
  };
}

#endif
