#ifndef RODIN_MODELS_DITANCE_EIKONAL_H
#define RODIN_MODELS_DITANCE_EIKONAL_H

#include "Rodin/Geometry/PolytopeIterator.h"
#include "Rodin/Models/Eikonal/FMM.h"

#include "Base.h"
#include "Rodin/Variational/ForwardDecls.h"
#include <functional>

namespace Rodin::Models::Distance
{
  template <class FES, class Data>
  class Eikonal : public Base<Eikonal<FES, Data>>
  {
    public:
      using SolutionType = Variational::GridFunction<FES, Data>;

      Eikonal(SolutionType& u)
        : m_u(u)
      {}

      Eikonal& solve()
      {
        static thread_local const auto s_speed =
          [](const Geometry::Point&) -> Real { return 1.0; };

        auto& u = m_u.get();
        const auto& fes = u.getFiniteElementSpace();
        const auto& interior = this->getInterior();
        const auto& interface = this->getInterface();
        const auto& mesh = fes.getMesh();

        Models::Eikonal::FMM fmm(u, s_speed);

        m_visited.resize(mesh.getVertexCount(), false);

        Geometry::FaceIterator it;
        if (interior.empty())
          it = mesh.getFace();
        else
          it = mesh.getBoundary();
        for (auto it = mesh.getFace(); !it.end(); ++it)
        {
          const auto& face = *it;
          if (interface.find(face.getAttribute()) == interface.end())
            continue;
          for (const auto& vertex : face.getVertices())
          {
            if (m_visited[vertex])
            {
              continue;
            }
            else
            {
              m_seeds.push_back(vertex);
              m_visited[vertex] = true;
            }
          }
        }

        fmm.seed(m_seeds).solve();

        return *this;
      };

      Eikonal& sign()
      {
        auto& u = m_u.get();
        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const size_t d = mesh.getDimension();
        std::fill(m_visited.begin(), m_visited.end(), false);
        for (auto it = mesh.getCell(); !it.end(); ++it)
        {
          const auto& cell = *it;
          const auto cellAttr = cell.getAttribute();
          const bool isInterior = this->getInterior().contains(cellAttr);
          if (isInterior)
          {
            for (const auto& vertex : cell.getVertices())
            {
              decltype(auto) fe = fes.getFiniteElement(d, vertex);
              for (size_t local = 0; local < fe.getCount(); local++)
              {
                const Index global = fes.getGlobalIndex({d, cell.getIndex()}, local);
                if (m_visited[global])
                  continue;
                u[global] *= -1;
                m_visited[global] = true;
              }
            }
          }
        }
        return *this;
      }

    private:
      std::reference_wrapper<SolutionType> m_u;
      std::vector<uint8_t> m_visited;
      std::vector<Index> m_seeds;
  };
}

#endif
