#ifndef RODIN_GEOMETRY_SHARD_H
#define RODIN_GEOMETRY_SHARD_H

#include "Mesh.h"

namespace Rodin::Geometry
{
  class Shard final : public Mesh<Context::Local>
  {
    public:
      using Parent = Mesh<Rodin::Context::Local>;
      using Context = typename Parent::Context;

      /**
       * @brief Class used to build SubMesh<Context::Local> instances.
       */
      class Builder
      {
        public:
          Builder() = default;

          Builder& initialize(const Mesh<Context::Local>& parent);

          Builder& include(size_t d, Index parentIdx);

          Builder& include(size_t d, const IndexSet& indices);

          Shard finalize();

        private:
          std::optional<std::reference_wrapper<const Mesh<Context>>> m_parent;
          Mesh<Context>::Builder m_build;
          std::vector<Index> m_sidx;
          std::vector<boost::bimap<Index, Index>> m_s2ps;
          size_t m_dimension;
      };

      Shard();

      Shard(const Shard& other);

      Shard(Shard&& other);

      Shard& operator=(const Shard&) = delete;

      Shard& operator=(Shard&& other)
      {
        if (this != &other)
        {
          Parent::operator=(std::move(other));
          m_s2ps = std::move(other.m_s2ps);
        }
        return *this;
      }

      /**
       * @brief Gets the map of polytope indices from the SubMesh to the parent
       * Mesh.
       */
      const boost::bimap<Index, Index>& getPolytopeMap(size_t d) const
      {
        return m_s2ps.at(d);
      }

    private:
      std::vector<boost::bimap<Index, Index>> m_s2ps;
  };
}

#endif

