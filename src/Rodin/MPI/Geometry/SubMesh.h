/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MPI_GEOMETRY_SUBMESH_H
#define RODIN_MPI_GEOMETRY_SUBMESH_H

/**
 * @file
 * @brief Distributed SubMesh specialization for MPI contexts.
 *
 * This file defines @ref Rodin::Geometry::SubMesh<Rodin::Context::MPI>, a
 * submesh that stores a rank-local shard subset of a distributed parent mesh
 * and maintains bidirectional index mappings between submesh-local and
 * parent-local shard indices.
 */

#include "Rodin/Geometry/SubMesh.h"
#include "Rodin/Geometry/Shard.h"

#include "Mesh.h"

namespace Rodin::Geometry
{
  /**
   * @ingroup SubMeshSpecializations
   * @brief SubMesh representing a distributed subregion of a parent MPI mesh.
   *
   * A `SubMesh<Context::MPI>` is an MPI mesh that represents a subset of
   * polytopes from a parent `Mesh<Context::MPI>`, maintaining bidirectional
   * index mappings between submesh-local shard indices and parent-local shard
   * indices.
   *
   * Each rank independently builds its portion of the submesh by including
   * entities from its parent shard. The submesh shard inherits the ownership
   * state (Owned, Shared, Ghost) from the parent shard for each included
   * entity, so global counting and reductions remain correct.
   *
   * The distributed vertex indices of the submesh shard are taken directly
   * from the parent shard, so no additional reconciliation is needed for
   * vertices. Only the inclusion call order matters.
   *
   * # Usage
   *
   * ## Creating a distributed boundary submesh:
   * @code{.cpp}
   * SubMesh<Context::MPI>::Builder builder;
   * builder.initialize(mpiMesh);
   *
   * // Include all owned boundary faces on this rank
   * for (auto it = mpiMesh.getBoundary(); it; ++it)
   *   builder.include(mpiMesh.getDimension() - 1, it->getIndex());
   *
   * SubMesh<Context::MPI> boundary = builder.finalize();
   * @endcode
   *
   * @see SubMeshBase, Mesh<Context::MPI>, SubMesh<Context::Local>
   */
  template <>
  class SubMesh<Context::MPI> final : public SubMeshBase, public Mesh<Context::MPI>
  {
    public:
      /**
       * @brief Parent mesh type.
       */
      using Parent = Mesh<Context::MPI>;

      /**
       * @brief Builder for constructing distributed SubMesh instances.
       *
       * Provides a fluent interface for incrementally building submeshes from
       * a parent `Mesh<Context::MPI>` by selecting specific polytopes. Each
       * rank includes entities independently using parent-shard-local indices.
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
           * @param[in] parent Parent distributed mesh to extract submesh from.
           * @returns Reference to this builder.
           */
          Builder& initialize(const Mesh<Context::MPI>& parent);

          /**
           * @brief Includes a single polytope in the submesh.
           *
           * The index @p parentLocalIdx is the rank-local shard index of the
           * polytope in the parent mesh. The ownership state (Owned, Shared,
           * Ghost) is inherited from the parent shard.
           *
           * For polytopes of dimension @p d > 0, all associated vertices are
           * automatically included with their parent-shard ownership state.
           *
           * @param[in] d Topological dimension of the polytope.
           * @param[in] parentLocalIdx Local shard index in the parent mesh.
           * @returns Reference to this builder for method chaining.
           */
          Builder& include(size_t d, Index parentLocalIdx);

          /**
           * @brief Finalizes construction and returns the distributed submesh.
           *
           * After finalization:
           * - The submesh shard uses the parent shard's distributed vertex
           *   indices for consistent global operations.
           * - Ownership and halo metadata are inherited from the parent shard.
           * - `getPolytopeMap(d)` provides the submesh-local to parent-local
           *   mapping for restriction operations.
           *
           * @returns Newly constructed SubMesh<Context::MPI>.
           */
          SubMesh finalize();

        private:
          Optional<std::reference_wrapper<const Mesh<Context::MPI>>> m_parent;
          Shard::Builder m_shardBuilder;
          std::vector<SubMeshBase::PolytopeMap> m_s2ps;
          std::vector<Index> m_sidx;
          size_t m_dimension = 0;
      };

      /**
       * @brief Constructs a submesh with a parent distributed mesh.
       * @param[in] parent Reference to the parent MPI mesh.
       */
      explicit
      SubMesh(std::reference_wrapper<const Mesh<Context::MPI>> parent);

      /**
       * @brief Copy constructor.
       * @param[in] other SubMesh to copy from.
       */
      SubMesh(const SubMesh& other);

      /**
       * @brief Move constructor.
       * @param[in] other SubMesh to move from.
       */
      SubMesh(SubMesh&& other);

      /**
       * @brief Copy assignment (deleted).
       */
      SubMesh& operator=(const SubMesh&) = delete;

      /**
       * @brief Move assignment operator.
       * @param[in] other SubMesh to move from.
       * @returns Reference to this object.
       */
      SubMesh& operator=(SubMesh&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_parent = std::move(other.m_parent);
          m_s2ps = std::move(other.m_s2ps);
          m_ancestors = std::move(other.m_ancestors);
        }
        return *this;
      }

      /**
       * @brief Restricts a point from the parent mesh to this submesh.
       * @param[in] p Point in parent mesh coordinates.
       * @returns Optional point in submesh (empty if point not in submesh).
       */
      Optional<Point> restriction(const Point& p) const override;

      /**
       * @brief Checks if this mesh is a submesh.
       * @returns Always true for SubMesh instances.
       */
      bool isSubMesh() const override
      {
        return true;
      }

      /**
       * @brief Checks whether a point is local to this submesh or its shard.
       * @param[in] p Point to test.
       * @returns True if the point belongs to this submesh or its local shard.
       */
      bool isLocalPoint(const Point& p) const override;

      /**
       * @brief Gets the immediate parent mesh.
       * @returns Reference to the parent distributed mesh.
       */
      const Mesh<Context::MPI>& getParent() const override;

      /**
       * @brief Gets the polytope index mapping for a dimension.
       * @param[in] d Topological dimension.
       * @returns Bidirectional mapping between submesh-local and parent-local indices.
       *
       * The mapping allows efficient queries:
       * - `left[i]`: parent shard-local index for submesh shard-local index i
       * - `right[j]`: submesh shard-local index for parent shard-local index j
       */
      const PolytopeMap& getPolytopeMap(size_t d) const override
      {
        return m_s2ps.at(d);
      }

      /**
       * @brief Gets all ancestor meshes.
       * @returns Deque of ancestor mesh references (parent, grandparent, etc.)
       */
      const Deque<Ancestor>& getAncestors() const override
      {
        return m_ancestors;
      }

      /**
       * @brief Downcasts to SubMeshBase.
       * @returns Reference to this as SubMeshBase.
       */
      SubMeshBase& asSubMesh() override
      {
        return *this;
      }

      /**
       * @brief Downcasts to SubMeshBase (const).
       * @returns Const reference to this as SubMeshBase.
       */
      const SubMeshBase& asSubMesh() const override
      {
        return *this;
      }

    private:
      std::reference_wrapper<const Mesh<Context::MPI>> m_parent;  ///< Parent mesh reference
      std::vector<PolytopeMap> m_s2ps;                            ///< Submesh-to-parent index maps
      Deque<Ancestor> m_ancestors;                                 ///< Ancestor mesh chain
  };
}

#endif
