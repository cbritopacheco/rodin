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
   * @defgroup SubMeshSpecializations SubMesh Template Specializations
   * @brief Template specializations of the SubMesh class.
   * @see SubMesh
   */

  class SubMeshBase
  {
    public:
      /**
       * @brief Represents the inclusion of a Point of a SubMesh @f$ C @f$ into
       * a Point of a Mesh @f$ P @f$.
       */
      virtual Point inclusion(const Point& p) const = 0;

      /**
       * @brief Represents the restriction of a Point of a Mesh @f$ P @f$ into
       * a Point of a SubMesh @f$ C @f$.
       */
      virtual Point restriction(const Point& p) const = 0;

      /**
       * @returns Reference to the parent Mesh object
       */
      virtual const MeshBase& getParent() const = 0;

      /**
       * @brief Gets the map of polytope indices from the SubMesh to the parent
       * Mesh.
       */
      virtual const boost::bimap<Index, Index>& getPolytopeMap(size_t d) const = 0;
  };

  /**
   * @ingroup SubMeshSpecializations
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
   *   auto& submesh = static_cast<SubMesh<Context::Sequential>&>(mesh);
   * }
   * @endcode
   *
   * Alternatively, you can use the asSubMesh() to access the SubMeshBase
   * interface:
   * @code{.cpp}
   * if (mesh.isSubMesh())
   * {
   *   // Cast is well defined
   *   auto& submesh = mesh.asSubMesh();
   * }
   * @endcode
   *
   */
  template <>
  class SubMesh<Context::Sequential> final : public SubMeshBase, public Mesh<Context::Sequential>
  {
    public:
      using Parent = Mesh<Context::Sequential>;

      /**
       * @brief Class used to build SubMesh<Context::Sequential> instances.
       */
      class Builder
      {
        public:
          Builder() = default;

          Builder& initialize(const Mesh<Context::Sequential>& parent);

          Builder& include(size_t d, Index parentIdx);

          Builder& include(size_t d, const IndexSet& indices);

          SubMesh finalize();

        private:
          std::optional<std::reference_wrapper<const Mesh<Context::Sequential>>> m_parent;
          Mesh<Context::Sequential>::Builder m_build;
          std::vector<Index> m_sidx;
          std::vector<boost::bimap<Index, Index>> m_s2ps;
          size_t m_dimension;
      };

      explicit
      SubMesh(std::reference_wrapper<const Mesh<Context::Sequential>> parent);

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

      Point inclusion(const Point& p) const override;

      Point restriction(const Point& p) const override;

      bool isSubMesh() const override
      {
        return true;
      }

      /**
       * @returns Reference to the parent Mesh object
       */
      const Mesh<Context::Sequential>& getParent() const override;

      /**
       * @brief Gets the map of polytope indices from the SubMesh to the parent
       * Mesh.
       */
      inline
      const boost::bimap<Index, Index>& getPolytopeMap(size_t d) const override
      {
        return m_s2ps.at(d);
      }

    private:
      std::reference_wrapper<const Mesh<Context::Sequential>> m_parent;
      std::vector<boost::bimap<Index, Index>> m_s2ps;
  };
}

#endif

