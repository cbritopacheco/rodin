#ifndef RODIN_GEOMETRY_SHARD_H
#define RODIN_GEOMETRY_SHARD_H

#include <boost/bimap/vector_of.hpp>

#include "Rodin/Pair.h"

#include "SubMesh.h"

namespace Rodin::Geometry
{
  class Shard : public Mesh<Context::Local>
  {
    friend class boost::serialization::access;

    public:
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

          Builder& include(size_t d, Index parentIdx);

          Builder& ghost(size_t d, Index parentIdx);

          Builder& include(size_t d, const IndexSet& indices);

          Shard finalize();

        private:
          std::optional<std::reference_wrapper<const Mesh<Context>>> m_parent;
          Mesh<Context>::Builder m_build;
          std::vector<Index> m_sidx;
          std::vector<PolytopeMap> m_s2ps;
          std::vector<IndexSet> m_ghosts;
          size_t m_dimension;
      };

      Shard() = default;

      Shard(const Shard& other);

      Shard(Shard&& other);

      Shard& operator=(Shard&& other);

      bool isGhost(size_t d, Index idx) const;

      const PolytopeMap& getPolytopeMap(size_t d) const;

      template<class Archive>
      void serialize(Archive& ar, const unsigned int version)
      {
        ar & boost::serialization::base_object<Mesh<Context>>(*this);
        ar & m_s2ps;
        ar & m_ghosts;
      }

      const auto& getGhosts() const
      {
        return m_ghosts;
      }

      Real getVolume() const override;

      Real getVolume(Attribute attr) const override;

      Real getVolume(const FlatSet<Attribute>& attr) const override;

      Real getPerimeter() const override;

      Real getPerimeter(Attribute attr) const override;

      Real getPerimeter(const FlatSet<Attribute>& attr) const override;

      Real getArea() const override;

      Real getArea(Attribute attr) const override;

      Real getArea(const FlatSet<Attribute>& attr) const override;

      Real getMeasure(size_t d) const override;

      Real getMeasure(size_t d, Attribute attr) const override;

      Real getMeasure(size_t d, const FlatSet<Attribute>& attr) const override;

    private:
      std::vector<PolytopeMap> m_s2ps;
      std::vector<IndexSet> m_ghosts;
  };
}

#endif
