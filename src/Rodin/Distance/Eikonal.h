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
#ifndef RODIN_MODELS_DISTANCE_EIKONAL_H
#define RODIN_MODELS_DISTANCE_EIKONAL_H

#include <functional>

#include "Rodin/Eikonal/FMM.h"

#include "Base.h"
#include "Rodin/Variational/ForwardDecls.h"

namespace Rodin::Distance
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
        const auto& fes  = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const auto& interior  = this->getInterior();
        const auto& interface = this->getInterface();

        Rodin::Eikonal::FMM fmm(u, s_speed);

        m_visited.assign(mesh.getVertexCount(), 0);
        m_seeds.clear();

        for (auto it = mesh.getFace(); !it.end(); ++it)
        {
          const auto& face = *it;
          const auto a = face.getAttribute(); // Optional<Attribute>
          if (!a || !interface.contains(*a))
            continue;

          for (Index vtx : face.getVertices())
          {
            const size_t v = static_cast<size_t>(vtx);
            if (m_visited[v])
              continue;
            m_seeds.push_back(vtx);
            m_visited[v] = 1;
          }
        }

        if (m_seeds.empty())
        {
          Alert::MemberFunctionException(*this, __func__)
            << "No seed vertices found for the interface. "
            << Alert::Raise;
        }

        fmm.seed(m_seeds).solve();
        return *this;
      }

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
        const auto& fes  = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const auto& interior = this->getInterior();

        assert(fes.getSize() == mesh.getVertexCount());
        assert(fes.getVectorDimension() == 1);

        m_visited.assign(mesh.getVertexCount(), 0);

        for (auto it = mesh.getCell(); !it.end(); ++it)
        {
          const auto& cell = *it;
          const auto a = cell.getAttribute();
          if (!a || !interior.contains(*a))
            continue;

          for (const Index vtx : cell.getVertices())
          {
            const size_t v = static_cast<size_t>(vtx);
            if (m_visited[v])
              continue;

            u[vtx] *= -1;      // dof == vertex
            m_visited[v] = 1;
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
