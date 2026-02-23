/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file LinearElasticity.h
 * @brief Linear elasticity model for solid mechanics.
 *
 * This file provides the LinearElasticity class, which implements the linear
 * elasticity equations for small deformation analysis in solid mechanics.
 */
#ifndef RODIN_MODELS_ELASTICITY_LINEARELASTICITY_H
#define RODIN_MODELS_ELASTICITY_LINEARELASTICITY_H

#include "Rodin/Solver/CG.h"
#include "Rodin/Variational/LinearElasticity.h"

namespace Rodin::Elasticity
{
  /**
   * @brief Linear elasticity system for small deformations.
   *
   * This class implements the linear elasticity equations in the context of
   * small deformation theory. The governing equations are:
   * @f[
   *   -\nabla \cdot \boldsymbol{\sigma}(\mathbf{u}) = \mathbf{f}
   * @f]
   * where @f$ \mathbf{u} @f$ is the displacement field, @f$ \mathbf{f} @f$ is
   * the body force, and the stress tensor is given by:
   * @f[
   *   \boldsymbol{\sigma} = \lambda (\nabla \cdot \mathbf{u}) \mathbf{I} +
   *   2\mu \boldsymbol{\varepsilon}(\mathbf{u})
   * @f]
   * with the strain tensor:
   * @f[
   *   \boldsymbol{\varepsilon}(\mathbf{u}) = \frac{1}{2}(\nabla \mathbf{u} +
   *   (\nabla \mathbf{u})^T)
   * @f]
   *
   * ## Parameters
   * - @f$ \lambda @f$: Lamé's first parameter (related to bulk modulus)
   * - @f$ \mu @f$: Lamé's second parameter (shear modulus)
   *
   * ## Material Relations
   * The Lamé parameters can be expressed in terms of Young's modulus @f$ E @f$
   * and Poisson's ratio @f$ \nu @f$:
   * @f[
   *   \lambda = \frac{E\nu}{(1+\nu)(1-2\nu)}, \quad
   *   \mu = \frac{E}{2(1+\nu)}
   * @f]
   *
   * @tparam FES Finite element space type (must have vector range)
   * @tparam Lambda Type of Lamé's first parameter
   * @tparam Mu Type of Lamé's second parameter (shear modulus)
   *
   * ## Usage Example
   * ```cpp
   * // Vector-valued FE space for displacement
   * VectorP1 fes(mesh);
   * RealFunction lambda(1.0);
   * RealFunction mu(0.5);
   * LinearElasticity<VectorP1, RealFunction, RealFunction> elasticity(lambda, mu, fes);
   * ```
   */
  template <class FES, class Lambda, class Mu>
  class LinearElasticity
  {
    using FESRange = typename FormLanguage::Traits<FES>::RangeType;
    static_assert(std::is_same_v<FESRange, Math::Vector<Real>>);

    public:
      /**
       * @brief Constructs a linear elasticity model.
       *
       * @param[in] l Lamé's first parameter @f$ \lambda @f$
       * @param[in] m Shear modulus @f$ \mu @f$
       * @param[in] fes Finite element space for the displacement field
       */
      LinearElasticity(const Lambda& l, const Mu& m, const FES& fes)
        : m_lei(l, m, fes)
      {}

      /**
       * @brief Constructs a linear elasticity model with function parameters.
       *
       * Overload accepting FunctionBase objects for spatially varying
       * material parameters.
       *
       * @tparam MuDerived Derived type for shear modulus function
       * @tparam LambdaDerived Derived type for Lamé's first parameter function
       * @param[in] lambda Lamé's first parameter function
       * @param[in] mu Shear modulus function
       * @param[in] fes Finite element space for the displacement field
       */
      template <class MuDerived, class LambdaDerived>
      LinearElasticity(
          const Variational::FunctionBase<LambdaDerived>& lambda,
          const Variational::FunctionBase<MuDerived>& mu,
          const FES& fes)
        : LinearElasticity(lambda, mu, fes)
      {}

    private:
      Variational::LinearElasticityIntegrator<FES, Lambda, Mu> m_lei;  ///< Linear elasticity integrator
  };
}

#endif





