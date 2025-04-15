#ifndef RODIN_GEOMETRY_SHARD_H
#define RODIN_GEOMETRY_SHARD_H

#include "Rodin/Pair.h"

#include "SubMesh.h"

namespace Rodin::Geometry
{
  class Shard : public Mesh<Context::Local>
  {
    friend class boost::serialization::access;

    public:
      class Builder
      {
        public:
          Builder() = default;

          Builder& initialize(const Mesh<Context>& parent);

          Builder& include(size_t d, Index parentIdx);

          Builder& ghost(size_t d, Index parentIdx);

          Builder& include(size_t d, const IndexSet& indices);

          Shard finalize();

        private:
          std::optional<std::reference_wrapper<const Mesh<Context>>> m_parent;
          Mesh<Context>::Builder m_build;
          std::vector<Index> m_sidx;
          std::vector<boost::bimap<Index, Index>> m_s2ps;
          std::vector<IndexSet> m_ghosts;
          size_t m_dimension;
      };

      bool isGhost(size_t d, Index idx) const;

      const boost::bimap<Index, Index>& getPolytopeMap(size_t d) const;

      template<class Archive>
      void serialize(Archive& ar, const unsigned int version)
      {
        ar & boost::serialization::base_object<Mesh<Context>>(*this);
        ar & m_s2ps;
        ar & m_ghosts;
      }

    private:
      std::vector<boost::bimap<Index, Index>> m_s2ps;
      std::vector<IndexSet> m_ghosts;
  };
}

#endif
