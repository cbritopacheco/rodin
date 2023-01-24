/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MESH_MESHLOADER_H
#define RODIN_MESH_MESHLOADER_H

#include <map>
#include <optional>
#include <mfem.hpp>
#include <boost/filesystem.hpp>

#include "Rodin/Geometry/Mesh.h"
#include "Rodin/IO/Loader.h"

#include "ForwardDecls.h"
#include "Helpers.h"

namespace Rodin::IO
{
  /**
   * @defgroup MeshLoaderSpecializations MeshLoader Template Specializations
   * @brief Template specializations of the MeshLoader class.
   * @see MeshLoader
   */

  /**
   * @brief Base class for mesh loader objects.
   */
  template <class Trait>
  class MeshLoaderBase : public IO::Loader<Rodin::Geometry::Mesh<Trait>>
  {
    public:
      MeshLoaderBase(Rodin::Geometry::Mesh<Trait>& mesh)
        :  m_mesh(mesh),
          m_fixOrientation(true)
      {}

      MeshLoaderBase& fixOrientation(bool fixOrientation)
      {
        m_fixOrientation = fixOrientation;
        return *this;
      }

    protected:
      Rodin::Geometry::Mesh<Trait>& getObject() override
      {
        return m_mesh;
      }

      bool getFixOrientation() const
      {
        return m_fixOrientation;
      }

    private:
      Rodin::Geometry::Mesh<Context::Serial>& m_mesh;
      bool m_fixOrientation;
  };

  template <>
  class MeshLoader<IO::FileFormat::MFEM, Context::Serial>
    : public MeshLoaderBase<Context::Serial>
  {
    public:
      MeshLoader(Rodin::Geometry::Mesh<Context::Serial>& mesh)
        : MeshLoaderBase<Context::Serial>(mesh)
      {}

      void load(std::istream& is) override;

    private:
      std::string m_parseTag;
  };

  template <>
  class MeshLoader<IO::FileFormat::GMSH, Context::Serial>
    : public MeshLoaderBase<Context::Serial>
  {
    public:
      MeshLoader(Rodin::Geometry::Mesh<Context::Serial>& mesh)
        : MeshLoaderBase<Context::Serial>(mesh)
      {}

      void load(std::istream& is) override;
  };

  /**
   * @ingroup MeshLoaderSpecializations
   * @brief Specialization for loading Serial meshes in the MEDIT file format.
   *
   * The MEDIT file format specification can be found by visiting
   * <a href="https://www.ljll.math.upmc.fr/frey/logiciels/Docmedit.dir/index.html">this
   * link</a>.
   */
  template <>
  class MeshLoader<IO::FileFormat::MEDIT, Context::Serial>
    : public MeshLoaderBase<Context::Serial>
  {
    public:
      using Object = Rodin::Geometry::Mesh<Context::Serial>;

      struct IOContext
      {
        std::map<Medit::Keyword, std::optional<std::istream::pos_type>> pos;
        std::map<Medit::Keyword, size_t> count;
      };

      MeshLoader(Object& mesh)
        : MeshLoaderBase<Context::Serial>(mesh)
      {}

      void load(std::istream& is) override;

      IOContext& getContext()
      {
        return m_context;
      }

      const IOContext& getContext() const
      {
        return m_context;
      }

    private:
      IOContext m_context;
  };
}

#endif
