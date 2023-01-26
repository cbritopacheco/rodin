/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2023.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MESH_SUBMESH_H
#define RODIN_MESH_SUBMESH_H

#include <map>
#include <optional>
#include <functional>
#include <boost/bimap.hpp>

#include "ForwardDecls.h"
#include "Mesh.h"

namespace Rodin::Geometry
{
  /**
   * @brief A SubMesh object represents a subregion of a Mesh object.
   *
   * A SubMesh object contains a reference to the parent Mesh object. It also
   * contains information regarding the mapping of elements and vertices
   * between the child and parent Mesh.
   *
   * A Mesh which is also a SubMesh may be casted into down to access
   * the SubMesh functionality. For example:
   * @code{.cpp}
   * if (mesh.isSubMesh())
   * {
   *   // Cast is well defined
   *   auto& submesh = static_cast<SubMesh&>(mesh);
   * }
   * @endcode
   *
   */
  template <>
  class SubMesh<Context::Serial> : public Mesh<Context::Serial>
  {
    public:
      class Builder : public BuilderBase
      {
        public:
          Builder();

          Builder& setMesh(SubMesh<Context::Serial>& mesh)
          {
            m_mesh.emplace(mesh);
            m_indices.resize(mesh.getDimension() + 1);
            return *this;
          }

          Builder& include(size_t dim, std::set<Index> indices);

          void finalize() override;

        private:
          std::optional<std::reference_wrapper<SubMesh<Context::Serial>>> m_mesh;
          std::vector<std::set<Index>> m_indices;
      };

      SubMesh(const MeshBase& parent);

      SubMesh(const SubMesh& other);

      SubMesh(SubMesh&& other);

      SubMesh& operator=(const SubMesh&) = delete;

      SubMesh& operator=(SubMesh&& other)
      {
        MeshBase::operator=(std::move(other));
        m_mesh = std::move(other.m_mesh);
        m_parent = std::move(other.m_parent);
        m_s2ps = std::move(other.m_s2ps);
        return *this;
      }

      bool isSubMesh() const override
      {
        return true;
      }

      /**
       * @returns Reference to the parent Mesh object
       */
      const MeshBase& getParent() const;

      const boost::bimap<Index, Index>& getSimplexMap(size_t d) const
      {
        return m_s2ps.at(d);
      }

      [[deprecated]] const boost::bimap<Index, Index>& getElementMap() const
      {
        return m_s2ps.at(getDimension());
      }

      SubMesh<Context::Serial>::Builder initialize(size_t dim);

      Mesh<Context::Serial>::Builder initialize(size_t dim, size_t sdim) = delete;

    private:
      Mesh<Context::Serial> m_mesh;
      std::reference_wrapper<const MeshBase> m_parent;
      std::vector<boost::bimap<Index, Index>> m_s2ps;
  };
}

#endif

