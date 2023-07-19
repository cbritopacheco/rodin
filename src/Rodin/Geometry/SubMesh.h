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
   * @defgroup DotSpecializations Dot Template Specializations
   * @brief Template specializations of the Dot class.
   * @see Dot
   */

  /**
   * @ingroup DotSpecializations
   */

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
      using Parent = Mesh<Context::Serial>;

      /**
       * @brief Class used to build SubMesh<Context::Serial> instances.
       */
      class Builder
      {
        public:
          Builder() = default;

          Builder& initialize(const Mesh<Context::Serial>& parent);

          Builder& include(size_t d, Index parentIdx);

          Builder& include(size_t d, const IndexSet& indices);

          SubMesh finalize();

        private:
          std::optional<std::reference_wrapper<const Mesh<Context::Serial>>> m_parent;
          Mesh<Context::Serial>::Builder m_build;
          std::vector<Index> m_sidx;
          std::vector<boost::bimap<Index, Index>> m_s2ps;
      };

      explicit
      SubMesh(std::reference_wrapper<const Mesh<Context::Serial>> parent);

      SubMesh(const SubMesh& other);

      SubMesh(SubMesh&& other);

      SubMesh& operator=(const SubMesh&) = delete;

      SubMesh& operator=(SubMesh&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_parent = std::move(other.m_parent);
          m_s2ps = std::move(other.m_s2ps);
        }
        return *this;
      }

      bool isSubMesh() const override
      {
        return true;
      }

      /**
       * @returns Reference to the parent Mesh object
       */
      const Mesh<Context::Serial>& getParent() const;

      /**
       * @brief Gets the map of simplex indices from the submesh to the parent
       * mesh.
       */
      const boost::bimap<Index, Index>& getPolytopeMap(size_t d) const
      {
        return m_s2ps.at(d);
      }

    private:
      std::reference_wrapper<const Mesh<Context::Serial>> m_parent;
      std::vector<boost::bimap<Index, Index>> m_s2ps;
  };
}

#endif

