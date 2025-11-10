/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MESH_MESHLOADER_H
#define RODIN_MESH_MESHLOADER_H

#include <utility>
#include <fstream>

#include "Rodin/IO/Loader.h"
#include "Rodin/Geometry/Mesh.h"

#include "ForwardDecls.h"

namespace Rodin::IO
{
  /**
   * @defgroup MeshLoaderSpecializations MeshLoader Template Specializations
   * @brief Template specializations of the MeshLoader class.
   * @see MeshLoader
   */

  /**
   * @brief Base class for loading meshes from files or streams.
   *
   * MeshLoaderBase provides the foundation for loading mesh data in different
   * file formats. It extends the generic Loader class to handle mesh-specific
   * operations and manages the mesh object being populated.
   *
   * @tparam Context Context type (e.g., Context::Local for sequential execution)
   *
   * Specialized loaders for specific file formats should derive from this class
   * and implement the load(std::istream&) method to parse their respective formats.
   *
   * ## Usage Example
   * ```cpp
   * Mesh<Context::Local> mesh;
   * MeshLoader<FileFormat::MFEM, Context::Local> loader(mesh);
   * loader.load("mesh.mfem");
   * ```
   *
   * @see Loader, MeshPrinter
   */
  template <class Context>
  class MeshLoaderBase : public IO::Loader<Rodin::Geometry::Mesh<Context>>
  {
    public:
      /**
       * @brief Context type for the mesh (e.g., sequential or parallel).
       */
      using ContextType = Context;

      /**
       * @brief Type of mesh object being loaded.
       */
      using ObjectType = Geometry::Mesh<ContextType>;

      /**
       * @brief Parent Loader class type.
       */
      using Parent = IO::Loader<ObjectType>;

      /**
       * @brief Constructs a mesh loader for the given mesh object.
       * @param[in,out] mesh Mesh object to be populated with loaded data
       *
       * The mesh is stored by reference and will be modified during loading.
       */
      MeshLoaderBase(ObjectType& mesh)
        : m_mesh(mesh)
      {}

    protected:
      /**
       * @brief Gets a reference to the mesh being loaded.
       * @returns Reference to the mesh object
       *
       * Provides access to the mesh for derived classes to populate during loading.
       */
      ObjectType& getObject() override
      {
        return m_mesh.get();
      }

    private:
      std::reference_wrapper<ObjectType> m_mesh;
  };
}

#endif
