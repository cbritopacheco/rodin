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
   * @brief Advanced mesh shard for distributed parallel computing.
   *
   * A Shard is a specialized mesh that represents a portion of a larger mesh
   * in a distributed computing context. Unlike the simplified MeshShard, this
   * class provides full support for ghost polytopes and ownership tracking,
   * which are essential for parallel finite element computations.
   *
   * # Distributed Computing Concepts
   *
   * In parallel mesh decomposition, each process owns a subset of mesh elements:
   * - **Owned polytopes**: Elements that this process is responsible for computing
   * - **Ghost polytopes**: Elements owned by other processes but needed for boundary computations
   * - **Halo**: The set of processes that need ghost information from this shard
   *
   * # Mathematical Foundation
   *
   * For a mesh @f$ \mathcal{T}_h @f$ decomposed into @f$ N @f$ shards:
   * @f[
   *   \mathcal{T}_h = \bigcup_{i=1}^{N} \mathcal{S}_i
   * @f]
   * where each shard @f$ \mathcal{S}_i @f$ contains both owned and ghost elements.
   *
   * # Usage in MPI Applications
   *
   * @code{.cpp}
   * // Build a shard from a parent mesh
   * Shard::Builder builder;
   * builder.initialize(parentMesh);
   *
   * // Include owned polytopes
   * for (Index idx : ownedCells)
   *   builder.include({dimension, idx}, Shard::Flags::Owned);
   *
   * // Include ghost polytopes from neighboring processes
   * for (Index idx : ghostCells)
   *   builder.include({dimension, idx}, Shard::Flags::Ghost);
   *
   * Shard shard = builder.finalize();
   *
   * // Check polytope ownership
   * if (shard.isOwned(d, localIdx))
   *   // Compute on this element
   * @endcode
   *
   * # Thread Safety
   * Shard objects are not thread-safe during construction. Once finalized,
   * read-only operations are thread-safe.
   *
   * @see MeshShard, Mesh, Sharder
   */
  class Shard : public Mesh<Context::Local>
  {
    friend class boost::serialization::access;

    public:
      /**
       * @brief Bidirectional polytope index mapping.
       *
       * Maps between shard-local indices and parent mesh indices,
       * enabling efficient lookup in both directions.
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
       * @brief Builder for constructing Shard instances.
       *
       * Provides a fluent interface for incrementally building shards from
       * parent meshes with precise control over ownership and ghost polytopes.
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
           * @brief Includes a polytope with specified flags.
           * @param[in] p Pair of (dimension, parent index)
           * @param[in] flags Ownership flags (Owned or Ghost)
           * @returns Pair of (local index, was newly added)
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
           * @brief Gets owner map for a dimension (mutable).
           * @param[in] d Dimension
           * @returns Map from ghost index to owning process
           */
          UnorderedMap<Index, Index>& getOwner(size_t d);

          /**
           * @brief Gets halo map for a dimension (mutable).
           * @param[in] d Dimension
           * @returns Map from owned index to processes needing it
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
