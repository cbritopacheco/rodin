/**
 * @file Eikonal.h
 * @brief Eikonal equation-based distance function computation.
 *
 * This file provides the Eikonal class, which computes distance functions
 * by solving the Eikonal equation using the Fast Marching Method (FMM).
 * The Eikonal equation is given by:
 * @f[
 *   |\nabla d| = 1
 * @f]
 * where @f$ d @f$ is the distance function.
 */
#ifndef RODIN_MODELS_DITANCE_EIKONAL_H
#define RODIN_MODELS_DITANCE_EIKONAL_H

#include "Rodin/Geometry/PolytopeIterator.h"
#include "Rodin/Models/Eikonal/FMM.h"

#include "Base.h"
#include "Rodin/Variational/ForwardDecls.h"
#include <functional>

namespace Rodin::Models::Distance
{
  /**
   * @brief Distance function computation using the Eikonal equation.
   *
   * This class solves the Eikonal equation @f$ |\nabla d| = 1 @f$ to compute
   * distance functions on meshes using the Fast Marching Method. It can
   * optionally compute signed distance functions by specifying interior regions.
   *
   * @tparam FES Finite element space type
   * @tparam Data Data storage type for the solution
   *
   * ## Mathematical Background
   * The Eikonal equation describes the propagation of wavefronts with unit
   * speed, which is equivalent to computing distances. The Fast Marching Method
   * solves this equation efficiently on unstructured meshes.
   *
   * ## Usage Example
   * ```cpp
   * P1 fes(mesh);
   * GridFunction u(fes);
   * Eikonal eikonal(u);
   * eikonal.setInterface(interfaceAttr)
   *        .setInterior(interiorAttr)
   *        .solve()
   *        .sign();
   * ```
   */
  template <class FES, class Data>
  class Eikonal : public Base<Eikonal<FES, Data>>
  {
    public:
      /// Solution type: grid function containing the distance values
      using SolutionType = Variational::GridFunction<FES, Data>;

      /**
       * @brief Constructs an Eikonal distance solver.
       *
       * @param[in,out] u Grid function to store the computed distance values
       */
      Eikonal(SolutionType& u)
        : m_u(u)
      {}

      /**
       * @brief Solves the Eikonal equation to compute the distance function.
       *
       * This method seeds the Fast Marching Method at the interface region
       * and propagates the distance values throughout the domain.
       *
       * @return Reference to this object for method chaining
       */
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

      /**
       * @brief Computes the sign of the distance function.
       *
       * This method assigns negative values to the distance function in the
       * interior region, converting an unsigned distance to a signed distance
       * function. Must be called after solve().
       *
       * @return Reference to this object for method chaining
       *
       * @note The interior region must be set via setInterior() before calling
       * this method for correct signed distance computation.
       */
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
      std::reference_wrapper<SolutionType> m_u;   ///< Reference to the solution grid function
      std::vector<uint8_t> m_visited;              ///< Visited vertices tracker
      std::vector<Index> m_seeds;                  ///< Seed vertices for FMM
  };
}

#endif
