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

      RippleSimplexification(const MeshType& mesh)
        : Parent(mesh)
      {}

      MeshType simplexify()
      {
        const auto& mesh = getMesh();
        const size_t D = mesh.getDimension();
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, D, D - 1);
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, D - 1, D);
        RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, D - 1, D - 1);
        std::deque<Index> searchQueue;
        std::vector<size_t> visitedCount(mesh.getFaceCount(), 0);
        while (searchQueue.size() > 0)
        {
          const Index i = searchQueue.back();
          searchQueue.pop_back();
          auto it = mesh.getFace(i);
          // Check node compatibility
          bool compatible = true;
          const auto& cells = mesh.getConnectivity().getIncidence({D - 1, D}, i);
          for (const auto& j : cells)
          {
            auto cell = mesh.getCell(j);
            const auto& faces = mesh.getConnectivity().getIncidence({D, D - 1}, j);
            switch (cell->getGeometry())
            {
              case Polytope::Type::TriangularPrism:
              {
                // Check if the triangulations are compatible for the prism by
                // checking if at least two of the triangles in the
                // quadrilateral share a common vertex

                compatible = false;
                break;
              }
            }

            if (!compatible)
              break;
          }

          // If not compatible, then change the triangulation for this face
          if (!compatible)
          {
            m_triangulation[i] = visitedCount[i] % s_simplices[it->getGeometry()].size();
            for (auto adj = it->getAdjacent(); adj; ++adj)
            {
              const Index j = adj->getIndex();
              searchQueue.push_back(j);
            }
          }
        }
        return m_build.finalize();
      }

    private:
      static GeometryIndexed<std::vector<IndexArray>> s_simplices;
      std::vector<size_t> m_triangulation;
      MeshType::Builder m_build;
  };
}

#endif

