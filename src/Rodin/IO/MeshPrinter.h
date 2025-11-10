/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MESH_MESHPRINTER_H
#define RODIN_MESH_MESHPRINTER_H

#include <utility>

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/IO/Printer.h"

#include "ForwardDecls.h"

namespace Rodin::IO
{
  /**
   * @brief Base class for printing meshes to streams or files.
   *
   * MeshPrinterBase provides the foundation for writing mesh data in different
   * file formats. It extends the generic Printer class to handle mesh-specific
   * operations and manages the mesh object being output.
   *
   * @tparam Context Context type (e.g., Context::Local for sequential execution)
   *
   * Specialized printers for specific file formats should derive from this class
   * and implement the print(std::ostream&) method to generate their respective formats.
   *
   * ## Usage Example
   * ```cpp
   * const Mesh<Context::Local>& mesh = getMesh();
   * MeshPrinter<FileFormat::MEDIT, Context::Local> printer(mesh);
   * std::ofstream file("output.mesh");
   * printer.print(file);
   * ```
   *
   * @see Printer, MeshLoader
   */
  template <class Context>
  class MeshPrinterBase : public IO::Printer<Geometry::Mesh<Context>>
  {
    public:
      /**
       * @brief Context type for the mesh (e.g., sequential or parallel).
       */
      using ContextType = Context;

      /**
       * @brief Type of mesh object being printed.
       */
      using ObjectType = Geometry::Mesh<ContextType>;

      /**
       * @brief Constructs a mesh printer for the given mesh object.
       * @param[in] mesh Mesh object to be written to output
       *
       * The mesh is stored by const reference and remains unchanged during printing.
       */
      MeshPrinterBase(const ObjectType& mesh)
        : m_mesh(mesh)
      {}

      /**
       * @brief Gets a const reference to the mesh being printed.
       * @returns Const reference to the mesh object
       *
       * Provides read-only access to the mesh for derived classes during output.
       */
      const ObjectType& getObject() const override
      {
        return m_mesh.get();
      }

    private:
      std::reference_wrapper<const ObjectType> m_mesh;
  };
}

#endif

