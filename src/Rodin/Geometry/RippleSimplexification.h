#ifndef RODIN_GEOMETRY_RIPPLESIMPLEXIFICATION_H
#define RODIN_GEOMETRY_RIPPLESIMPLEXIFICATION_H

#include "Simplexification.h"

namespace Rodin::Geometry
{
  template <class T>
  class RippleSimplexification;

  template <>
  class RippleSimplexification<Mesh<Context::Local>>
    : public Simplexification<Mesh<Context::Local>>
  {
    public:
      using MeshType = Mesh<Context::Local>;
      using Parent = Simplexification<MeshType>;

      class Simplexification
      {
        public:
          constexpr
          Simplexification(Polytope::Type geometry)
            : m_geometry(geometry)
          {}

          constexpr
          size_t getCount() const;

          constexpr
          size_t getCount(size_t k, size_t d) const;

          constexpr
          const IndexArray& getSimplex(size_t k, size_t d, size_t i) const;

          constexpr
          Polytope::Type getGeometry() const;

        private:
          static GeometryIndexed<std::vector<std::vector<std::vector<IndexArray>>>> s_simplices;
          Polytope::Type m_geometry;
      };

      RippleSimplexification(const MeshType& mesh);

      MeshType simplexify();

    private:
      std::vector<size_t> m_triangulation;
      MeshType::Builder m_build;
  };

  RippleSimplexification(Mesh<Context::Local>) -> RippleSimplexification<Mesh<Context::Local>>;
}

#endif

