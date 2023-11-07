/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_MODELS_ELASTICITY_LINEARELASTICITY_H
#define RODIN_MODELS_ELASTICITY_LINEARELASTICITY_H

#include "Rodin/Solver/CG.h"
#include "Rodin/Variational/LinearElasticity.h"

namespace Rodin::Models::Elasticity
{
  /**
   * @brief Linear elasticity system
   */
  template <class FES, class Lambda, class Mu>
  class LinearElasticity
  {
    using FESRange = typename FormLanguage::Traits<FES>::RangeType;
    static_assert(std::is_same_v<FESRange, Math::Vector>);

    public:
      LinearElasticity(const Lambda& l, const Mu& m, const FES& fes)
        : m_lei(l, m, fes)
      {}

      template <class MuDerived, class LambdaDerived>
      LinearElasticity(
          const Variational::FunctionBase<LambdaDerived>& lambda,
          const Variational::FunctionBase<MuDerived>& mu,
          const FES& fes)
        : LinearElasticity(lambda, mu, fes)
      {}

    private:
      Variational::LinearElasticityIntegrator<FES, Lambda, Mu> m_lei;
  };
}

#endif





