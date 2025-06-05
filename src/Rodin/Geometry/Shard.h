#ifndef RODIN_GEOMETRY_SHARD_H
#define RODIN_GEOMETRY_SHARD_H

#include <type_traits>
#include <boost/bimap/vector_of.hpp>

#include "Rodin/Pair.h"
#include "Rodin/Serialization/FlatMap.h"
#include "Rodin/Serialization/BitSet.h"

#include "SubMesh.h"

namespace Rodin::Geometry
{
  class Shard : public Mesh<Context::Local>
  {
    friend class boost::serialization::access;

    public:
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

          template<class Archive>
          void serialize(Archive& ar, const unsigned int version)
          {
            ar & m_bits;
          }

        private:
          BitSet2 m_bits;
      };

      using PolytopeMap =
        boost::bimap<
          boost::bimaps::vector_of<Index>,
          boost::bimaps::unordered_set_of<Index>>;

      using ContextType = Rodin::Context::Local;

      using Parent = Mesh<ContextType>;

      class Builder
      {
        public:
          Builder() = default;

          Builder& initialize(const Mesh<Context>& parent);

          Builder& include(size_t d, Index parentIdx, const Flags& flags = Shard::Flags::None);

          Builder& flag(size_t d, Index parentIdx, const Flags& flags);

          Shard finalize();

          const PolytopeMap& getPolytopeMap(size_t d) const;

          size_t getPolytopeCount(size_t d) const;

        private:
          std::optional<std::reference_wrapper<const Mesh<Context>>> m_parent;
          Mesh<Context>::Builder m_build;
          std::vector<Index> m_sidx;
          std::vector<PolytopeMap> m_s2ps;
          std::vector<FlatMap<Index, Flags>> m_flags;
          size_t m_dimension;
      };

      Shard() = default;

      Shard(const Shard& other);

      Shard(Shard&& other);

      Shard& operator=(Shard&& other);

      /**
       * @brief Indicates whether the given polytope is a ghost element.
       * @param[in] d Dimension of the polytope
       * @param[in] idx Local index of the polytope
       */
      bool isGhost(size_t d, Index idx) const;

      const PolytopeMap& getPolytopeMap(size_t d) const;

      template<class Archive>
      void serialize(Archive& ar, const unsigned int version)
      {
        ar & boost::serialization::base_object<Mesh<Context>>(*this);
        ar & m_s2ps;
        ar & m_flags;
      }

      const FlatMap<Index, Flags>& getFlags(size_t d) const
      {
        return m_flags[d];
      }

    private:
      std::vector<PolytopeMap> m_s2ps;
      std::vector<FlatMap<Index, Flags>> m_flags;
  };
}

#endif
