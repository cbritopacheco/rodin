/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2023.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_SUBMESH_H
#define RODIN_GEOMETRY_SUBMESH_H

/**
 * @file
 * @brief SubMesh class for representing subregions of meshes with parent relationships.
 */

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

  /**
   * @brief Abstract base class for SubMesh implementations.
   *
   * Defines the interface for submeshes, which represent subregions of
   * parent meshes while maintaining bidirectional index mappings.
   */
  class SubMeshBase
  {
    public:
      /**
       * @brief Polytope index mapping structure.
       *
       * Maintains bidirectional mapping between submesh and parent mesh indices.
       */
      struct PolytopeMap
      {
        std::vector<Index> left;        ///< Submesh index to parent index
        FlatMap<Index, Index> right;    ///< Parent index to submesh index
      };

      /**
       * @brief Type alias for ancestor mesh references.
       */
      using Ancestor = std::reference_wrapper<const MeshBase>;

      /**
       * @brief Restricts a point from parent mesh to submesh coordinates.
       * @param[in] p Point in parent mesh
       * @returns Optional point in submesh (empty if point not in submesh)
       *
       * Projects a point @f$ p @f$ from the parent mesh @f$ \mathcal{T}_h @f$
       * to the submesh @f$ \mathcal{S} \subset \mathcal{T}_h @f$.
       */
      virtual Optional<Point> restriction(const Point& p) const = 0;

      /**
       * @brief Gets the immediate parent mesh.
       * @returns Reference to the parent mesh
       */
      virtual const MeshBase& getParent() const = 0;

      /**
       * @brief Gets all ancestor meshes.
       * @returns Deque of ancestor mesh references (parent, grandparent, etc.)
       *
       * Returns the chain of meshes from immediate parent to root mesh.
       */
      virtual const Deque<Ancestor>& getAncestors() const = 0;

      /**
       * @brief Gets the polytope index mapping for a dimension.
       * @param[in] d Dimension
       * @returns Bidirectional mapping between submesh and parent indices
       */
      virtual const PolytopeMap& getPolytopeMap(size_t d) const = 0;

      /**
       * @brief Equality comparison.
       * @param[in] other SubMesh to compare with
       * @returns True if this is the same object (pointer equality)
       */
      constexpr
      bool operator==(const SubMeshBase& other) const
      {
        return this == &other;
      }
  };

  /**
   * @ingroup SubMeshSpecializations
   * @brief SubMesh representing a subregion of a parent mesh.
   *
   * A SubMesh is a mesh that represents a subset of polytopes from a parent
   * mesh, maintaining bidirectional index mappings between child and parent.
   * This is useful for:
   * - Extracting boundary meshes
   * - Defining material regions
   * - Mesh refinement
   * - Domain decomposition
   *
   * # Mathematical Foundation
   *
   * For a parent mesh @f$ \mathcal{T}_h @f$ and submesh @f$ \mathcal{S} @f$:
   * @f[
   *   \mathcal{S} \subseteq \mathcal{T}_h
   * @f]
   * The submesh maintains mappings @f$ \phi: I_{\mathcal{S}} \to I_{\mathcal{T}_h} @f$
   * where @f$ I_{\mathcal{S}} @f$ and @f$ I_{\mathcal{T}_h} @f$ are index sets.
   *
   * # Usage Examples
   *
   * ## Creating a boundary submesh:
   * @code{.cpp}
   * SubMesh<Context::Local>::Builder builder;
   * builder.initialize(mesh);
   *
   * // Include all boundary faces
   * for (auto it = mesh.getFace(); it; ++it)
   *   if (it->isBoundary())
   *     builder.include(mesh.getDimension() - 1, it->getIndex());
   *
   * SubMesh<Context::Local> boundary = builder.finalize();
   * @endcode
   *
   * ## Accessing as submesh:
   * @code{.cpp}
   * if (mesh.isSubMesh())
   * {
   *   // Downcast to submesh
   *   auto& submesh = static_cast<SubMesh<Context::Local>&>(mesh);
   *   const auto& parent = submesh.getParent();
   * }
   *
   * // Or use base interface
   * if (mesh.isSubMesh())
   * {
   *   auto& submesh = mesh.asSubMesh();
   *   const auto& parent = submesh.getParent();
   * }
   * @endcode
   *
   * # Thread Safety
   * SubMesh objects are not thread-safe during construction. Once finalized,
   * read-only operations are thread-safe.
   *
   * @see SubMeshBase, Mesh, Shard
   */
  template <>
  class SubMesh<Context::Local> final : public SubMeshBase, public Mesh<Context::Local>
  {
    public:
      /**
       * @brief Parent mesh type.
       */
      using Parent = Mesh<Rodin::Context::Local>;
      
      /**
       * @brief Context type.
       */
      using Context = typename Parent::Context;

      /**
       * @brief Builder for constructing SubMesh instances.
       *
       * Provides a fluent interface for incrementally building submeshes
       * from parent meshes by selecting specific polytopes.
       */
      class Builder
      {
        public:
          /**
           * @brief Default constructor.
           */
          Builder() = default;

          /**
           * @brief Initializes builder with parent mesh.
           * @param[in] parent Parent mesh to extract submesh from
           * @returns Reference to this builder
           */
          Builder& initialize(const Mesh<Context>& parent);

          /**
           * @brief Includes a single polytope in the submesh.
           * @param[in] d Dimension of the polytope
           * @param[in] parentIdx Index in parent mesh
           * @returns Reference to this builder for method chaining
           */
          Builder& include(size_t d, Index parentIdx);

          /**
           * @brief Includes multiple polytopes in the submesh.
           * @param[in] d Dimension of the polytopes
           * @param[in] indices Vector of parent mesh indices
           * @returns Reference to this builder for method chaining
           */
          Builder& include(size_t d, const IndexVector& indices);

          /**
           * @brief Finalizes construction and returns the submesh.
           * @returns Newly constructed SubMesh
           */
          SubMesh finalize();

        private:
          Optional<std::reference_wrapper<const Mesh<Context>>> m_parent;
          Mesh<Context>::Builder m_build;
          std::vector<Index> m_sidx;
          std::vector<PolytopeMap> m_s2ps;
          size_t m_dimension;
      };

      /**
       * @brief Constructs a submesh with a parent.
       * @param[in] parent Reference to the parent mesh
       */
      explicit
      SubMesh(std::reference_wrapper<const Mesh<Context>> parent);

      /**
       * @brief Copy constructor.
       * @param[in] other SubMesh to copy from
       */
      SubMesh(const SubMesh& other);

      /**
       * @brief Move constructor.
       * @param[in] other SubMesh to move from
       */
      SubMesh(SubMesh&& other);

      /**
       * @brief Copy assignment (deleted).
       */
      SubMesh& operator=(const SubMesh&) = delete;

      /**
       * @brief Move assignment operator.
       * @param[in] other SubMesh to move from
       * @returns Reference to this object
       */
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

      /**
       * @brief Restricts a point to this submesh.
       * @param[in] p Point in parent mesh coordinates
       * @returns Optional point in submesh (empty if not in submesh)
       *
       * @note Implements virtual method from SubMeshBase.
       */
      Optional<Point> restriction(const Point& p) const override;

      /**
       * @brief Checks if this mesh is a submesh.
       * @returns Always true for SubMesh instances
       */
      bool isSubMesh() const override
      {
        return true;
      }

      /**
       * @brief Gets the immediate parent mesh.
       * @returns Reference to the parent mesh
       */
      const Mesh<Context>& getParent() const override;

      /**
       * @brief Gets the polytope index mapping.
       * @param[in] d Dimension
       * @returns Bidirectional mapping between submesh and parent indices
       *
       * The mapping allows efficient queries:
       * - `left[i]`: parent index for submesh index i
       * - `right[j]`: submesh index for parent index j
       */
      const PolytopeMap& getPolytopeMap(size_t d) const override
      {
        return m_s2ps.at(d);
      }

      /**
       * @brief Gets all ancestor meshes.
       * @returns Deque of ancestor mesh references
       *
       * Returns the full chain of parent meshes up to the root mesh.
       */
      const Deque<Ancestor>& getAncestors() const override
      {
        return m_ancestors;
      }

    private:
      std::reference_wrapper<const Mesh<Context>> m_parent;  ///< Parent mesh reference
      std::vector<PolytopeMap> m_s2ps;                      ///< Index mappings per dimension
      Deque<Ancestor> m_ancestors;                           ///< Ancestor mesh chain
  };
}

#endif

