#ifndef RODIN_GEOMETRY_SHARD_H
#define RODIN_GEOMETRY_SHARD_H

/**
 * @file
 * @brief Advanced mesh shard with ghost polytope tracking for distributed computing.
 */

#include "Rodin/Types.h"

#include "Rodin/Geometry/Mesh.h"

namespace Rodin::Geometry
{
  /**
   * @brief Local view of a distributed mesh partition.
   *
   * A `Shard` is a local mesh containing:
   * - entities owned by the current partition, and
   * - ghost entities copied from neighboring partitions.
   *
   * The class derives from `Mesh<Context::Local>` because, once constructed,
   * a shard behaves as an ordinary local mesh for geometric and topological
   * queries. What makes it special is the additional metadata describing how
   * its local entities relate to the global parent mesh and to neighboring
   * partitions.
   *
   * ## Local vs parent indices
   *
   * Every entity stored in a shard has a **local shard index**. This is the
   * index used by the inherited local mesh API.
   *
   * In addition, each shard entity corresponds to an entity of the parent mesh.
   * For every topological dimension @f$ d @f$, `Shard` stores a bidirectional map
   *
   * @f[
   *   i_{\mathrm{local}} \leftrightarrow i_{\mathrm{parent}}
   * @f]
   *
   * through `PolytopeMap`.
   *
   * Therefore:
   * - all inherited `Mesh<Context::Local>` methods use **local shard indices**,
   * - `getPolytopeMap(d)` is used to recover the corresponding **parent index**.
   *
   * ## Ownership semantics
   *
   * Each local shard entity is tagged with one of the following logical states:
   *
   * - **Owned**:
   *   this partition is the designated owner of the entity.
   *   Owned entities are the ones on which assembly, counting, and reductions
   *   should typically be performed.
   *
   * - **Ghost**:
   *   this entity is present locally only because it is needed by an owned
   *   neighboring entity. A ghost entity is owned by another partition.
   *
   * - **None**:
   *   this flag is used only during shard construction. It means that the
   *   entity is inserted in the current shard without changing the previously
   *   established ownership. In practice this occurs when an entity has already
   *   been assigned an owner and must simply be made present in another shard.
   *
   * After finalization, users normally reason only in terms of:
   * - `isOwned(d, i)`
   * - `isGhost(d, i)`
   *
   * ## Owner map
   *
   * For each dimension @f$ d @f$, `getOwner(d)` stores:
   *
   * @f[
   *   \text{ghost local index} \mapsto \text{owner rank}
   * @f]
   *
   * That is, this map is only meaningful for ghost entities. If a local entity
   * `i` of dimension `d` is ghost, then `getOwner(d).at(i)` gives the rank that
   * owns it.
   *
   * ## Halo map
   *
   * For each dimension @f$ d @f$, `getHalo(d)` stores:
   *
   * @f[
   *   \text{owned local index} \mapsto \{\text{neighbor ranks that also contain it}\}
   * @f]
   *
   * So the halo map is the dual of the owner map:
   * - `owner` is information attached to ghosts,
   * - `halo` is information attached to owned entities.
   *
   * In particular, if `k` is an owned local entity of dimension `d`, then
   * `getHalo(d).at(k)` is the set of remote ranks for which this entity appears
   * as a ghost.
   *
   * This information is useful for:
   * - distributed numbering,
   * - ghost synchronization,
   * - neighborhood-restricted communication,
   * - distributed connectivity discovery.
   *
   * ## Geometric/topological contents of a shard
   *
   * A shard is not just a set of cells. It generally contains lower-dimensional
   * entities as well, for example vertices, edges, and faces needed by the local
   * finite element space or by local topological queries.
   *
   * In particular, a finalized shard may contain:
   * - owned entities,
   * - ghost entities,
   * - local incidences inherited from the parent mesh and filtered to the local
   *   shard content.
   *
   * Thus `Shard` is intended to be a self-consistent local mesh with overlap.
   *
   * ## Construction model
   *
   * `Shard::Builder` constructs a shard incrementally from a parent mesh.
   * The parent mesh is assumed to provide the topology needed to identify the
   * entities being inserted.
   *
   * The typical construction pattern is:
   *
   * @code{.cpp}
   * Shard::Builder builder;
   * builder.initialize(parentMesh);
   *
   * // insert owned entities
   * builder.include({D, cellIdx}, Shard::Flags::Owned);
   *
   * // insert overlap entities
   * builder.include({D, neighborCellIdx}, Shard::Flags::Ghost);
   *
   * Shard shard = builder.finalize();
   * @endcode
   *
   * `include()` is idempotent with respect to the parent index:
   * if the entity was already inserted in the shard, it is not duplicated.
   * The returned pair is:
   * - the local shard index of the entity,
   * - a boolean telling whether the entity was newly inserted.
   *
   * ## Finalization guarantees
   *
   * After `Builder::finalize()`:
   * - all local entity indices are stable,
   * - inherited local mesh queries operate on shard-local indices,
   * - `getPolytopeMap(d)` maps local shard entities back to parent entities,
   * - `isOwned()` / `isGhost()` are available for every local entity,
   * - `getOwner(d)` is defined for ghost entities,
   * - `getHalo(d)` is defined for owned entities.
   *
   * ## Important convention
   *
   * A shard does **not** itself define a distributed global numbering of degrees
   * of freedom or mesh entities. It only provides:
   * - a local mesh with overlap,
   * - ownership metadata,
   * - parent/shard index correspondences.
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
       * @brief Bidirectional correspondence between shard-local and parent indices.
       *
       * For a fixed topological dimension @f$ d @f$:
       * - `left[local]` is the parent index of the local shard entity,
       * - `right[parent]` is the local shard index of that parent entity,
       *   if present in the shard.
       *
       * This map does not encode ownership; it only encodes index correspondence.
       */
      struct PolytopeMap
      {
        friend class boost::serialization::access;

        std::vector<Index> left;        ///< Shard index to parent index
        UnorderedMap<Index, Index> right;    ///< Parent index to shard index

        /**
         * @brief Serialization method.
         * @tparam Archive Archive type
         * @param[in,out] ar Archive object
         * @param[in] version Serialization version
         */
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
       * @brief Bitmask flags indicating polytope state.
       *
       * A polytope in a shard can be marked as:
       * - **Owned**: This process is responsible for computations on this element
       * - **Ghost**: Element is owned by another process but needed for boundary calculations
       */
      class Flags
      {
        friend class boost::serialization::access;

        public:
          static const Flags None;  ///< No flags set
          static const Flags Owned; ///< Polytope is owned by local process
          static const Flags Ghost; ///< Polytope is ghosted from remote process

          /**
           * @brief Default constructor (no flags).
           */
          Flags()
            : m_bits(None.m_bits)
          {}

          /**
           * @brief Constructs from bit pattern.
           * @param[in] bits Bit pattern for flags
           */
          Flags(BitSet2 bits)
            : m_bits(bits)
          {}

          /**
           * @brief Copy constructor.
           */
          Flags(const Flags& other)
            : m_bits(other.m_bits)
          {}

          /**
           * @brief Bitwise OR operator.
           * @param[in] rhs Right-hand side flags
           * @returns Combined flags
           */
          Flags operator|(const Flags& rhs) const
          {
            return m_bits | rhs.m_bits;
          }

          /**
           * @brief Checks if flags are set (bitwise AND).
           * @param[in] rhs Flags to check
           * @returns True if any of rhs flags are set
           */
          Boolean operator&(const Flags& rhs) const
          {
            return has(rhs);
          }

          /**
           * @brief Bitwise OR assignment.
           * @param[in] rhs Flags to add
           * @returns Reference to this object
           */
          Flags& operator|=(const Flags& rhs)
          {
            m_bits |= rhs.m_bits;
            return *this;
          }

          /**
           * @brief Copy assignment.
           * @param[in] other Flags to copy
           * @returns Reference to this object
           */
          Flags& operator=(const Flags& other)
          {
            m_bits = other.m_bits;
            return *this;
          }

          /**
           * @brief Checks if specific flags are set.
           * @param[in] flag Flags to check for
           * @returns True if any of the specified flags are set
           */
          Boolean has(Flags flag) const
          {
            return (m_bits & flag.m_bits).any();
          }

          /**
           * @brief Serialization method.
           * @tparam Archive Archive type
           * @param[in,out] ar Archive object
           * @param[in] version Serialization version
           */
          template <class Archive>
          void serialize(Archive& ar, const unsigned int version)
          {
            ar & m_bits;
          }

        private:
          BitSet2 m_bits;
      };

      /**
       * @brief Incremental builder for `Shard`.
       *
       * The builder accumulates a subset of the parent mesh entities together
       * with ownership metadata and finally produces a self-consistent local mesh.
       *
       * The builder maintains, per topological dimension:
       * - the parent-to-shard index maps,
       * - per-entity ownership flags,
       * - ghost owner information,
       * - halo information for owned entities.
       *
       * `include()` inserts an entity identified by its parent index into the
       * shard if it is not already present, and returns its shard-local index.
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
           * @param[in] parent Parent mesh to extract shard from
           * @returns Reference to this builder
           */
          Builder& initialize(const Mesh<Context>& parent);

          /**
           * @brief Inserts a parent entity into the shard.
           *
           * @param[in] p Pair `(d, parentIdx)` identifying an entity in the parent mesh.
           * @param[in] flags Local classification to assign if the entity is newly inserted.
           *
           * @return Pair `(localIdx, inserted)` where:
           * - `localIdx` is the shard-local index of the entity,
           * - `inserted` is `true` iff the entity was not already present.
           *
           * If the entity already exists in the shard, no duplication occurs and
           * the existing local index is returned.
           *
           * @note `Flags::None` is meaningful during construction when an entity
           *       is required locally but its ownership was already established
           *       elsewhere.
           */
          std::pair<Index, Boolean> include(const std::pair<size_t, Index>& p, const Flags& flags);

          /**
           * @brief Finalizes construction and returns the shard.
           * @returns Newly constructed Shard
           */
          Shard finalize();

          /**
           * @brief Gets polytope map for a dimension (const).
           * @param[in] d Dimension
           * @returns Polytope index mapping
           */
          const PolytopeMap& getPolytopeMap(size_t d) const;

          /**
           * @brief Gets polytope count for a dimension.
           * @param[in] d Dimension
           * @returns Number of polytopes
           */
          size_t getPolytopeCount(size_t d) const;

          /**
           * @brief Returns the owner-rank map for ghost entities of dimension `d`.
           *
           * The map stores:
           *
           * @f[
           *   \text{ghost local index} \mapsto \text{owner rank}
           * @f]
           *
           * and is only meaningful for ghost entities.
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
           * and is only meaningful for owned entities.
           */
          UnorderedMap<Index, IndexSet>& getHalo(size_t d);

          /**
           * @brief Gets owner map for a dimension (const).
           * @param[in] d Dimension
           * @returns Map from ghost index to owning process
           */
          const UnorderedMap<Index, Index>& getOwner(size_t d) const;

          /**
           * @brief Gets halo map for a dimension (const).
           * @param[in] d Dimension
           * @returns Map from owned index to processes needing it
           */
          const UnorderedMap<Index, IndexSet>& getHalo(size_t d) const;

        private:
          Optional<std::reference_wrapper<const Mesh<Context>>> m_parent;
          Mesh<Context>::Builder m_build;
          std::vector<Index> m_sidx;
          std::vector<PolytopeMap> m_s2ps;
          std::vector<std::vector<Flags>> m_flags;
          std::vector<UnorderedMap<Index, IndexSet>> m_halo;
          std::vector<UnorderedMap<Index, Index>> m_owner;
          size_t m_dimension;
      };

      /**
       * @brief Default constructor.
       */
      Shard() = default;

      /**
       * @brief Copy constructor.
       * @param[in] other Shard to copy from
       */
      Shard(const Shard& other);

      /**
       * @brief Move constructor.
       * @param[in] other Shard to move from
       */
      Shard(Shard&& other);

      /**
       * @brief Move assignment operator.
       * @param[in] other Shard to move from
       * @returns Reference to this object
       */
      Shard& operator=(Shard&& other);

      /**
       * @brief Checks if a polytope is a ghost element.
       * @param[in] d Dimension of the polytope
       * @param[in] idx Local index of the polytope
       * @returns True if the polytope is ghosted from another process
       *
       * Ghost elements are owned by another process but included in this
       * shard for boundary calculations.
       */
      bool isGhost(size_t d, Index idx) const;

      /**
       * @brief Checks if a polytope is owned by this shard.
       * @param[in] d Dimension of the polytope
       * @param[in] idx Local index of the polytope
       * @returns True if the polytope is owned by this process
       *
       * Owned elements are those this process is responsible for computing.
       */
      bool isOwned(size_t d, Index idx) const;

      /**
       * @brief Gets the owner map for a dimension.
       * @param[in] d Dimension
       * @returns Map from ghost polytope index to owning process ID
       */
      const UnorderedMap<Index, Index>& getOwner(size_t d) const;

      /**
       * @brief Gets the halo map for a dimension.
       * @param[in] d Dimension
       * @returns Map from owned polytope index to set of processes needing it
       *
       * The halo identifies which neighboring processes need ghost copies
       * of this shard's owned elements.
       */
      const UnorderedMap<Index, IndexSet>& getHalo(size_t d) const;

      /**
       * @brief Gets the polytope index map for a dimension.
       * @param[in] d Dimension
       * @returns Bidirectional mapping between shard and parent indices
       */
      const PolytopeMap& getPolytopeMap(size_t d) const;

      /**
       * @brief Serialization method for Boost.Serialization.
       * @tparam Archive Archive type
       * @param[in,out] ar Archive object
       * @param[in] version Serialization version
       */
      template <class Archive>
      void serialize(Archive& ar, const unsigned int version)
      {
        ar & boost::serialization::base_object<Mesh<Context>>(*this);
        ar & m_s2ps;
        ar & m_flags;
        ar & m_owner;
        ar & m_halo;
      }

    private:
      std::vector<PolytopeMap> m_s2ps;                  ///< Index mappings per dimension
      std::vector<std::vector<Flags>> m_flags;          ///< Ownership flags per dimension
      std::vector<UnorderedMap<Index, Index>> m_owner;       ///< Ghost to owner process map
      std::vector<UnorderedMap<Index, IndexSet>> m_halo;     ///< Owned to needing processes map
  };
}

#endif
