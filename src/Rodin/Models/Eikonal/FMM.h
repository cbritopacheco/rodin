/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MODELS_EIKONAL_FMM_H
#define RODIN_MODELS_EIKONAL_FMM_H

#include "Rodin/Context/Local.h"
#include "Rodin/Variational/P1/P1.h"

namespace Rodin::Models::Distance
{
  template <class Solution, class SpeedFunction, class Context>
  class FMM;

  /**
   * @brief Poisson approximation to the distance function.
   */
  template <class Solution, class SpeedFunction>
  class FMM<Solution, SpeedFunction, Context::Local>
  {
    public:
      using ScalarType = Real;
      using Context = Context::Local;
      using Mesh = Geometry::Mesh<Context>;
      using FES = Variational::P1<ScalarType, Mesh>;

      using SolutionType = Solution;
      using SpeedFunctionType = SpeedFunction;

      enum class Label
      {
        Far,
        Considered,
        Accepted
      };

      template <class Callable>
      FMM(SolutionType& u, Callable&& speed)
        : m_u(u), m_speed(std::forward<Callable>(speed))
      {}

      template <class Pred>
      FMM& setInterface(Pred&& p)
      {
      }

      void solve()
      {
        auto& u = m_u.get();
        const auto& fes = u.getFiniteElementSpace();
        const auto& mesh = fes.getMesh();
        const auto& meshDim = mesh.getDimension();
        m_labels.resize(fes.getSize(), Label::Far);
        u = std::numeric_limits<Real>::infinity();
        return u;
      }

    private:
      std::reference_wrapper<SolutionType> m_u;
      SpeedFunctionType m_speed;
      std::vector<Label> m_labels;
  };
}

#endif




