#include "GeometryIndexed.h"
#include "RippleSimplexification.h"

namespace Rodin::Geometry
{
  GeometryIndexed<std::vector<std::vector<IndexArray>>>
  RippleSimplexification<Mesh<Context::Local>>::s_simplices =
  {
    {
      Polytope::Type::Triangle,
      {
        { {{0, 1, 2}} },
      }
    },
    {
      Polytope::Type::Quadrilateral,
      {
        { {{0, 1, 2}}, {{1, 2, 3}} },
        { {{0, 1, 3}}, {{0, 2, 3}} }
      }
    },
    {
      Polytope::Type::TriangularPrism,
      {
        { {{0, 1, 2}}, {{1, 2, 3}} },
        { {{0, 1, 3}}, {{0, 2, 3}} }
      }
    },
  };

  RippleSimplexification<Mesh<Context::Local>>
  ::RippleSimplexification(const Mesh<Context::Local>& mesh)
    : Parent(mesh)
  {}

  Mesh<Context::Local> RippleSimplexification<Mesh<Context::Local>>::simplexify()
  {
    const auto& mesh = getMesh();
    const size_t D = mesh.getDimension();
    RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, D, D - 1);
    RODIN_GEOMETRY_REQUIRE_INCIDENCE(mesh, D - 1, D);
    std::deque<Index> searchQueue;
    std::vector<size_t> visitedCount(mesh.getFaceCount(), 0);
    m_triangulation.resize(mesh.getFaceCount(), 0);
    for (auto it = mesh.getFace(); it; ++it)
      searchQueue.push_back(it->getIndex());
    while (searchQueue.size() > 0)
    {
      const Index i = searchQueue.back();
      searchQueue.pop_back();
      visitedCount[i] += 1;
      auto it = mesh.getFace(i);

      // Check if the face triangulation is compatible with the adjacent cells
      bool compatible;

      if (it->getGeometry() == Polytope::Type::Triangle)
        continue;

      const auto& cells = mesh.getConnectivity().getIncidence({D - 1, D}, i);

      for (const auto& j : cells)
      {
        auto cell = mesh.getCell(j);
        IndexMap<size_t> count;
        count.reserve(6);
        switch (cell->getGeometry())
        {
          case Polytope::Type::TriangularPrism:
          {
            compatible = false;
            const auto& faces = mesh.getConnectivity().getIncidence({D, D - 1}, j);
            for (const Index k : faces)
            {
              const auto& triangulation = m_triangulation[k];
              const auto& vs = mesh.getConnectivity().getPolytope(D - 1, k);
              const auto& simplices = s_simplices[mesh.getGeometry(D - 1, k)].at(triangulation);
              for (const auto& simplex : simplices)
              {
                for (const Index sv : simplex)
                  count[vs[sv]] += 1;
              }
            }

            for (const auto& [k, v] : count)
            {
              // One vertex must be shared by at least 5 triangles
              if (v >= 5)
              {
                compatible = true;
                break;
              }
            }
            break;
          }
        }

        if (!compatible)
          break;
      }

      // If not compatible, then change the triangulation and add the faces in
      // the incident cells to the search queue
      if (!compatible)
      {
        m_triangulation[i] = visitedCount[i] % s_simplices[it->getGeometry()].size();
        for (const auto& j : cells)
        {
          auto cell = mesh.getCell(j);
          const auto& faces = mesh.getConnectivity().getIncidence({D, D - 1}, j);
          for (const Index k : faces)
          {
            if (k == i)
              continue;
            searchQueue.push_back(k);
          }
        }
      }
    }

    // Build the new mesh
    for (auto it = mesh.getCell(); it; ++it)
    {
      const auto& cell = *it;
      const auto& faces = mesh.getConnectivity().getIncidence({D, D - 1}, cell.getIndex());
      for (const Index i : faces)
      {
        const auto& triangulation = m_triangulation[i];
        const auto& simplices = s_simplices[mesh.getGeometry(D - 1, i)].at(triangulation);
        for (const auto& simplex : simplices)
        {
        }
      }
    }

    return m_build.finalize();
  }
}
