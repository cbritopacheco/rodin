/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2026.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_GEOMETRY_SHARD_H
#define RODIN_GEOMETRY_SHARD_H

/**
 * @file
 * @brief Local mesh shard with partition, ownership, and overlap metadata.
 */

#include "Rodin/Types.h"
#include "Rodin/Geometry/Mesh.h"

namespace Rodin::Geometry
{
  /**
   * @brief Local view of one partition of a distributed mesh.
   *
   * A `Shard` is a local mesh enriched with metadata describing how its local
   * entities relate to a distributed partitioning of a larger mesh.
   *
   * A shard contains:
   * - entities that belong to the local partition, either as locally owned
   *   entities or as shared entities owned by another partition, and
   * - ghost entities imported as overlap from neighboring partitions.
   *
   * The class derives from `Mesh<Context::Local>` because, once constructed,
   * a shard behaves as an ordinary local mesh for geometric and topological
   * queries. What makes it special is the additional distributed metadata:
   * - local-to-distributed index correspondences,
   * - per-entity local state (`Owned`, `Shared`, `Ghost`),
   * - owner information for non-owned entities,
   * - halo information for owned entities.
   *
   * ## Local vs distributed indices
   *
   * Every entity stored in a shard has a **local shard index**. This is the
   * index used by the inherited local mesh API.
   *
   * In addition, each shard entity corresponds to an entity of the distributed
   * mesh. For every topological dimension @f$ d @f$, `Shard` stores a
   * bidirectional map
   *
   * @f[
   *   i_{\mathrm{local}} \leftrightarrow i_{\mathrm{distributed}}
   * @f]
   *
   * through `PolytopeMap`.
   *
   * Therefore:
   * - all inherited `Mesh<Context::Local>` methods use **local shard indices**,
   * - `getPolytopeMap(d)` is used to recover the corresponding distributed index.
   *
   * ## Local state semantics
   *
   * Each local shard entity is classified by two orthogonal notions:
   * - whether it belongs to the local partition, and
   * - which partition owns it.
   *
   * The following three mutually exclusive states are used:
   *
   * - **Owned**:
   *   The entity belongs to the local partition and this partition is its
   *   designated owner.
   *
   *   Owned entities are the ones on which assembly, counting, and global
   *   reductions should typically be performed.
   *
   * - **Shared**:
   *   The entity belongs to the local partition but is owned by another
   *   partition.
   *
   *   Shared entities are part of the partition’s geometric/topological domain,
   *   but this rank is not responsible for their global contribution.
   *
   * - **Ghost**:
   *   The entity does not belong to the local partition but is present locally
   *   as overlap because it is needed by neighboring local entities, for example
   *   for stencil support, adjacency completeness, or distributed connectivity.
   *
   *   A ghost entity is always owned by another partition.
   *
   * In particular:
   * - `Owned ∪ Shared` are exactly the entities that define the local partition,
   * - `Ghost` entities form the overlap with neighboring partitions,
   * - every distributed entity present on multiple partitions has exactly one
   *   owner.
   *
   * After finalization, users normally reason in terms of:
   * - `isOwned(d, i)`
   * - `isShared(d, i)`
   * - `isGhost(d, i)`
   * - `isLocal(d, i)` (equivalent to `Owned || Shared`)
   *
   * ## Owner map
   *
   * For each dimension @f$ d @f$, `getOwner(d)` stores
   *
   * @f[
   *   \text{local non-owned index} \mapsto \text{owner rank}
   * @f]
   *
   * That is, this map is meaningful for entities in state `Shared` or `Ghost`.
   *
   * In particular:
   * - if `i` is `Owned`, then `getOwner(d)` need not contain `i`,
   * - if `i` is `Shared` or `Ghost`, then `getOwner(d).at(i)` gives the rank
   *   that owns the entity.
   *
   * ## Halo map
   *
   * For each dimension @f$ d @f$, `getHalo(d)` stores
   *
   * @f[
   *   \text{owned local index} \mapsto \{\text{remote ranks that also contain it}\}
   * @f]
   *
   * So the halo map is dual to the owner map:
   * - `owner` is attached to non-owned local entities (`Shared` or `Ghost`),
   * - `halo` is attached to owned local entities.
   *
   * If `k` is an owned local entity of dimension `d`, then `getHalo(d).at(k)`
   * is the set of remote ranks that also contain the same distributed entity,
   * either as `Shared` or as `Ghost`.
   *
   * This information is useful for:
   * - distributed numbering,
   * - synchronization,
   * - neighborhood-restricted communication,
   * - distributed connectivity discovery.
   *
   * ## Geometric and topological contents
   *
   * A shard is not just a set of cells. It generally contains lower-dimensional
   * entities as well, for example vertices, edges, and faces needed by the local
   * finite element space or by local topological queries.
   *
   * In particular, a finalized shard may contain:
   * - owned entities,
   * - shared entities,
   * - ghost entities,
   * - local incidences inherited from the distributed mesh and filtered to the
   *   shard contents.
   *
   * Thus `Shard` is intended to be a self-consistent local mesh with overlap.
   *
   * ## Construction model
   *
   * `Shard::Builder` constructs a shard incrementally from either:
   * - a parent mesh (`Mode::Parent`), or
   * - explicit direct insertion of vertices and polytopes (`Mode::Direct`).
   *
   * In parent-based mode, a typical construction pattern is:
   *
   * @code{.cpp}
   * Shard::Builder builder;
   * builder.initialize(parentMesh);
   *
   * // insert owned entities
   * builder.include({D, cellIdx}, Shard::State::Owned);
   *
   * // insert overlap entities
   * builder.include({D, neighborCellIdx}, Shard::State::Ghost);
   *
   * Shard shard = builder.finalize();
   * @endcode
   *
   * `include()` is idempotent with respect to the parent index:
   * if the entity is already present in the shard, it is not duplicated.
   * The returned pair is:
   * - the local shard index of the entity,
   * - a boolean telling whether the entity was newly inserted.
   *
   * ## Finalization guarantees
   *
   * After `Builder::finalize()`:
   * - all local entity indices are stable,
   * - inherited local mesh queries operate on shard-local indices,
   * - `getPolytopeMap(d)` maps shard-local entities back to distributed or
   *   parent indices,
   * - `isOwned()`, `isShared()`, `isGhost()`, and `isLocal()` are available
   *   for every local entity,
   * - `getOwner(d)` is defined for non-owned local entities (`Shared` or
   *   `Ghost`),
   * - `getHalo(d)` is defined for owned local entities.
   *
   * ## Important convention
   *
   * A shard does **not** itself define a distributed global numbering of mesh
   * entities or degrees of freedom. It only provides:
   * - a local mesh with overlap,
   * - partition/ownership metadata,
   * - distributed-to-local and local-to-distributed index correspondences.
   *
   * Global distributed numbering is a higher-level concern handled for example
   * by MPI finite element spaces.
   *
   * @see Mesh<Context::Local>, Sharder, MPIMesh
   */
  class Shard : public Mesh<Context::Local>
  {
    friend class boost::serialization::access;

    public:
      /**
       * @brief Bidirectional correspondence between local shard and distributed indices.
       *
       * For a fixed topological dimension @f$ d @f$:
       * - `left[local]` is the distributed index of the local shard entity,
       * - `right[distributed]` is the local shard index of that distributed entity,
       *   if present in the shard.
       *
       * This map does not encode ownership; it only encodes index correspondence.
       */
      struct PolytopeMap
      {
        friend class boost::serialization::access;

        std::vector<Index> left;             ///< Local shard index -> distributed index
        UnorderedMap<Index, Index> right;    ///< Distributed index -> local shard index

        template <class Archive>
        void serialize(Archive& ar, const unsigned int version)
        {
          ar & left;
          ar & right;
        }
      };

      /**
       * @brief Context type (always local for shards).
       */
      using ContextType = Rodin::Context::Local;

      /**
       * @brief Parent mesh type.
       */
      using Parent = Mesh<ContextType>;

      /**
       * @brief Local state of a shard entity.
       */
      enum class State : uint8_t
      {
        Shared = 0, ///< In local partition, owned remotely
        Owned  = 1, ///< In local partition, owned locally
        Ghost  = 2  ///< Not in local partition, overlap only
      };

      /**
       * @brief Incremental builder for `Shard`.
       *
       * The builder accumulates a subset of mesh entities together with state,
       * owner, and halo metadata, and finally produces a self-consistent local mesh.
       *
       * The builder maintains, per topological dimension:
       * - the distributed-to-shard and shard-to-distributed index maps,
       * - the per-entity local state,
       * - owner information for non-owned entities,
       * - halo information for owned entities.
       *
       * `include()` inserts an entity identified by its parent index into the
       * shard if it is not already present, and returns its shard-local index.
       */
      class Builder
      {
        public:
          enum class Mode
          {
            None,    ///< Uninitialized builder
            Parent,  ///< Initialize from parent mesh
            Direct   ///< Initialize with explicit dimension and explicit contents
          };

          /**
           * @brief Default constructor.
           */
          Builder();

          /**
           * @brief Initializes the builder from a parent mesh.
           * @param[in] parent Parent mesh from which shard entities are drawn.
           * @returns Reference to this builder.
           */
          Builder& initialize(const Mesh<Context>& parent);

          /**
           * @brief Initializes the builder for direct construction.
           * @param[in] dimension Mesh topological dimension.
           * @param[in] sdim Ambient space dimension.
           * @returns Reference to this builder.
           */
          Builder& initialize(size_t dimension, size_t sdim);

          /**
           * @brief Inserts a vertex in direct-construction mode.
           * @param[in] globalIdx Distributed/global vertex index.
           * @param[in] x Vertex coordinates.
           * @param[in] state Local shard state.
           * @returns Local shard vertex index.
           */
          Index vertex(
              Index globalIdx,
              const Math::SpatialPoint& x,
              const State& state);

          /**
           * @brief Inserts a non-vertex polytope in direct-construction mode.
           * @param[in] d Topological dimension.
           * @param[in] globalIdx Distributed/global polytope index.
           * @param[in] g Geometry type.
           * @param[in] vs Local shard vertex indices.
           * @param[in] state Local shard state.
           * @returns Local shard polytope index.
           */
          Index polytope(
              size_t d,
              Index globalIdx,
              Polytope::Type g,
              const IndexArray& vs,
              const State& state);

          /**
           * @brief Records the owner rank of a non-owned local entity.
           * @param[in] d Topological dimension.
           * @param[in] localIdx Local shard index.
           * @param[in] ownerRank Owning rank.
           * @returns Reference to this builder.
           */
          Builder& setOwner(size_t d, Index localIdx, Index ownerRank);

          /**
           * @brief Records that an owned local entity is also present on a remote rank.
           * @param[in] d Topological dimension.
           * @param[in] localIdx Local shard index.
           * @param[in] neighborRank Remote rank containing the same entity.
           * @returns Reference to this builder.
           */
          Builder& halo(size_t d, Index localIdx, Index neighborRank);

          /**
           * @brief Sets the attribute of a local entity.
           * @param[in] p Pair `(d, localIdx)`.
           * @param[in] attr Optional attribute.
           * @returns Reference to this builder.
           */
          Builder& attribute(const std::pair<size_t, Index>& p, const Optional<Attribute>& attr);

          /**
           * @brief Inserts a parent entity into the shard in parent-based mode.
           *
           * @param[in] p Pair `(d, parentIdx)` identifying an entity in the parent mesh.
           * @param[in] state Local state to assign if the entity is newly inserted.
           *
           * @return Pair `(localIdx, inserted)` where:
           * - `localIdx` is the shard-local index of the entity,
           * - `inserted` is `true` iff the entity was not already present.
           *
           * If the entity is already present, no duplication occurs and the existing
           * local index is returned.
           */
          std::pair<Index, Boolean> include(const std::pair<size_t, Index>& p, const State& state);

          /**
           * @brief Finalizes construction and returns the resulting shard.
           * @returns Newly constructed shard.
           */
          Shard finalize();

          /**
           * @brief Returns the local/distributed index map for dimension `d`.
           * @param[in] d Topological dimension.
           * @returns Const reference to the polytope map.
           */
          const PolytopeMap& getPolytopeMap(size_t d) const;

          /**
           * @brief Returns the number of local entities of dimension `d`.
           * @param[in] d Topological dimension.
           * @returns Local entity count.
           */
          size_t getPolytopeCount(size_t d) const;

          /**
           * @brief Returns the owner map for dimension `d`.
           *
           * The map stores:
           *
           * @f[
           *   \text{local non-owned index} \mapsto \text{owner rank}
           * @f]
           *
           * and is meaningful for entities in state `Shared` or `Ghost`.
           */
          UnorderedMap<Index, Index>& getOwner(size_t d);

          /**
           * @brief Returns the halo map for owned entities of dimension `d`.
           *
           * The map stores:
           *
           * @f[
           *   \text{owned local index} \mapsto \{\text{remote ranks containing the same entity}\}
           * @f]
           *
           * and is meaningful only for owned local entities.
           */
          UnorderedMap<Index, IndexSet>& getHalo(size_t d);

          /**
           * @brief Returns the owner map for dimension `d` (const).
           * @param[in] d Topological dimension.
           * @returns Const reference to the owner map.
           */
          const UnorderedMap<Index, Index>& getOwner(size_t d) const;

          /**
           * @brief Returns the halo map for dimension `d` (const).
           * @param[in] d Topological dimension.
           * @returns Const reference to the halo map.
           */
          const UnorderedMap<Index, IndexSet>& getHalo(size_t d) const;

        private:
          Optional<std::reference_wrapper<const Mesh<Context>>> m_parent;
          Mesh<Context>::Builder m_build;
          std::vector<Index> m_sidx;
          std::vector<PolytopeMap> m_s2ds;
          std::vector<std::vector<State>> m_state;
          std::vector<UnorderedMap<Index, IndexSet>> m_halo;
          std::vector<UnorderedMap<Index, Index>> m_owner;
          size_t m_dimension;
          size_t m_sdim;

          Mode m_mode;
          PointCloud m_vertices;
          Connectivity<Context> m_connectivity;
          AttributeIndex m_attributes;
          PolytopeTransformationIndex m_transformations;
      };

      /**
       * @brief Default constructor.
       */
      Shard() = default;

      /**
       * @brief Copy constructor.
       * @param[in] other Shard to copy from.
       */
      Shard(const Shard& other);

      /**
       * @brief Move constructor.
       * @param[in] other Shard to move from.
       */
      Shard(Shard&& other);

      /**
       * @brief Move assignment operator.
       * @param[in] other Shard to move from.
       * @returns Reference to `*this`.
       */
      Shard& operator=(Shard&& other);

      /**
       * @brief Tests whether a local entity is ghost.
       * @param[in] d Topological dimension.
       * @param[in] idx Local shard index.
       * @returns `true` iff the entity is in state `Ghost`.
       */
      bool isGhost(size_t d, Index idx) const;

      /**
       * @brief Tests whether a local entity is owned.
       * @param[in] d Topological dimension.
       * @param[in] idx Local shard index.
       * @returns `true` iff the entity is in state `Owned`.
       */
      bool isOwned(size_t d, Index idx) const;

      /**
       * @brief Tests whether a local entity is shared.
       * @param[in] d Topological dimension.
       * @param[in] idx Local shard index.
       * @returns `true` iff the entity is in state `Shared`.
       */
      bool isShared(size_t d, Index idx) const;

      /**
       * @brief Tests whether a local entity belongs to the local partition.
       * @param[in] d Topological dimension.
       * @param[in] idx Local shard index.
       * @returns `true` iff the entity is `Owned` or `Shared`.
       */
      bool isLocal(size_t d, Index idx) const;

      /**
       * @brief Returns the owner map for dimension `d`.
       * @param[in] d Topological dimension.
       * @returns Map from local non-owned index to owner rank.
       */
      const UnorderedMap<Index, Index>& getOwner(size_t d) const;

      /**
       * @brief Returns the owner map for dimension `d`.
       * @param[in] d Topological dimension.
       * @returns Mutable map from local non-owned index to owner rank.
       */
      UnorderedMap<Index, Index>& getOwner(size_t d);

      /**
       * @brief Returns the halo map for dimension `d`.
       * @param[in] d Topological dimension.
       * @returns Map from owned local index to remote ranks containing the same entity.
       */
      const UnorderedMap<Index, IndexSet>& getHalo(size_t d) const;

      /**
       * @brief Returns the halo map for dimension `d`.
       * @param[in] d Topological dimension.
       * @returns Mutable map from owned local index to remote ranks containing the same entity.
       */
      UnorderedMap<Index, IndexSet>& getHalo(size_t d);

      /**
       * @brief Returns the local state array for dimension `d`.
       * @param[in] d Topological dimension.
       * @returns Const reference to per-entity states.
       */
      const std::vector<State>& getState(size_t d) const;

      /**
       * @brief Returns the local state array for dimension `d`.
       * @param[in] d Topological dimension.
       * @returns Mutable reference to per-entity states.
       */
      std::vector<State>& getState(size_t d);

      /**
       * @brief Returns the local/distributed polytope map for dimension `d`.
       * @param[in] d Topological dimension.
       * @returns Const reference to the polytope map.
       */
      const PolytopeMap& getPolytopeMap(size_t d) const;

      /**
       * @brief Returns the local/distributed polytope map for dimension `d`.
       * @param[in] d Topological dimension.
       * @returns Mutable reference to the polytope map.
       */
      PolytopeMap& getPolytopeMap(size_t d);

      /**
       * @brief Serialization method.
       * @tparam Archive Archive type.
       * @param[in,out] ar Archive object.
       * @param[in] version Serialization version.
       */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version)
      {
        ar & boost::serialization::base_object<Mesh<Context>>(*this);
        ar & m_s2ds;
        ar & m_state;
        ar & m_owner;
        ar & m_halo;
      }

    private:
      std::vector<PolytopeMap> m_s2ds;                  ///< Local/distributed index maps per dimension
      std::vector<std::vector<State>> m_state;          ///< Local state per dimension
      std::vector<UnorderedMap<Index, Index>> m_owner;  ///< Local non-owned index -> owner rank
      std::vector<UnorderedMap<Index, IndexSet>> m_halo;///< Owned local index -> remote ranks containing the entity
  };
}

#endif
