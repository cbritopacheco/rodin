#ifndef RODIN_GEOMETRY_SHARD_H
#define RODIN_GEOMETRY_SHARD_H

#include "Rodin/Types.h"

#include "Rodin/Geometry/Mesh.h"

namespace Rodin::Geometry
{
  class Shard : public Mesh<Context::Local>
  {
    friend class boost::serialization::access;

    public:
      struct PolytopeMap
      {
        friend class boost::serialization::access;

        std::vector<Index> left;
        FlatMap<Index, Index> right;

        template <class Archive>
        void serialize(Archive& ar, const unsigned int version)
        {
          ar & left;
          ar & right;
        }
      };

      using ContextType = Rodin::Context::Local;

      using Parent = Mesh<ContextType>;

      /**
       * @brief Bitmask enum indicating the state of a polytope in the shard.
       *
       * A polytope can be:
       * - Owned locally
       * - Ghosted from a remote process
       */
      class Flags
      {
        friend class boost::serialization::access;

        public:
          static const Flags None; ///< No flags set
          static const Flags Owned; ///< Polytope is owned by the local process
          static const Flags Ghost; ///< Polytope is ghosted from a remote process

          Flags()
            : m_bits(None.m_bits)
          {}

          Flags(BitSet2 bits)
            : m_bits(bits)
          {}

          Flags(const Flags& other)
            : m_bits(other.m_bits)
          {}

          Flags operator|(const Flags& rhs) const
          {
            return m_bits | rhs.m_bits;
          }

          Boolean operator&(const Flags& rhs) const
          {
            return has(rhs);
          }

          Flags& operator|=(const Flags& rhs)
          {
            m_bits |= rhs.m_bits;
            return *this;
          }

          Flags& operator=(const Flags& other)
          {
            m_bits = other.m_bits;
            return *this;
          }

          Boolean has(Flags flag) const
          {
            return (m_bits & flag.m_bits).any();
          }

          template <class Archive>
          void serialize(Archive& ar, const unsigned int version)
          {
            ar & m_bits;
          }

        private:
          BitSet2 m_bits;
      };

      class Builder
      {
        public:
          Builder() = default;

          Builder& initialize(const Mesh<Context>& parent);

          std::pair<Index, Boolean> include(const std::pair<size_t, Index>& p, const Flags& flags);

          Shard finalize();

          const PolytopeMap& getPolytopeMap(size_t d) const;

          size_t getPolytopeCount(size_t d) const;

          FlatMap<Index, Index>& getOwner(size_t d);

          FlatMap<Index, IndexSet>& getHalo(size_t d);

          const FlatMap<Index, Index>& getOwner(size_t d) const;

          const FlatMap<Index, IndexSet>& getHalo(size_t d) const;

        private:
          Optional<std::reference_wrapper<const Mesh<Context>>> m_parent;
          Mesh<Context>::Builder m_build;
          std::vector<Index> m_sidx;
          std::vector<PolytopeMap> m_s2ps;
          std::vector<std::vector<Flags>> m_flags;
          std::vector<FlatMap<Index, IndexSet>> m_halo;
          std::vector<FlatMap<Index, Index>> m_owner;
          size_t m_dimension;
      };

      /**
       * @brief Default constructor.
       */
      Shard() = default;

      /**
       * @brief Copy constructor.
       * @param[in] other The shard to copy from
       */
      Shard(const Shard& other);

      /**
       * @brief Move constructor.
       * @param[in] other The shard to move from
       */
      Shard(Shard&& other);

      Shard& operator=(Shard&& other);

      /**
       * @brief Indicates whether the given polytope is a ghost element.
       * @param[in] d Dimension of the polytope
       * @param[in] idx Local index of the polytope
       */
      bool isGhost(size_t d, Index idx) const;

      /**
       * @brief Indicates whether the given polytope is owned by the shard.
       * @param[in] d Dimension of the polytope
       * @param[in] idx Local index of the polytope
       */
      bool isOwned(size_t d, Index idx) const;

      const FlatMap<Index, Index>& getOwner(size_t d) const;

      const FlatMap<Index, IndexSet>& getHalo(size_t d) const;

      const PolytopeMap& getPolytopeMap(size_t d) const;

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
      std::vector<PolytopeMap> m_s2ps;
      std::vector<std::vector<Flags>> m_flags;
      std::vector<FlatMap<Index, Index>> m_owner;
      std::vector<FlatMap<Index, IndexSet>> m_halo;
  };
}

#endif
