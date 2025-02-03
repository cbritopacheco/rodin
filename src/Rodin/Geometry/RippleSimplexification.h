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

      RippleSimplexification(const MeshType& mesh);

      MeshType simplexify();

    private:
      static GeometryIndexed<std::vector<std::vector<IndexArray>>> s_simplices;
      std::vector<size_t> m_triangulation;
      MeshType::Builder m_build;
  };

  RippleSimplexification(Mesh<Context::Local>) -> RippleSimplexification<Mesh<Context::Local>>;
}

#endif

