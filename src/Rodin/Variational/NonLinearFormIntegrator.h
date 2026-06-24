/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2025.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file NonLinearFormIntegrator.h
 * @brief CRTP mixin providing the nonlinear form integrator interface.
 *
 * Any class that models a nonlinear weak form contribution should inherit from
 * this mixin and implement:
 *   - `Residual(v)`    — returns a LinearFormIntegrator
 *   - `Tangent(u, v)`  — returns a BilinearFormIntegrator
 *
 * The mixin then provides `operator()(u, v)` for free, which returns a
 * ProblemBody containing both contributions — suitable for direct use in a
 * Problem assignment expression.
 */
#ifndef RODIN_VARIATIONAL_NON_LINEAR_FORM_INTEGRATOR_H
#define RODIN_VARIATIONAL_NON_LINEAR_FORM_INTEGRATOR_H

namespace Rodin::Variational
{
  /**
   * @brief CRTP mixin for nonlinear form integrators.
   *
   * Derived classes must implement:
   * @code
   *   template <class Test>
   *   auto Residual(const Test& v) const;   // LinearFormIntegrator
   *
   *   template <class Trial, class Test>
   *   auto Tangent(const Trial& u, const Test& v) const;  // BilinearFormIntegrator
   * @endcode
   *
   * This mixin then provides:
   * @code
   *   template <class Trial, class Test>
   *   auto operator()(const Trial& u, const Test& v) const;  // ProblemBody
   * @endcode
   *
   * so that `obj(u, v)` can be used directly in a Problem assignment,
   * registering both the residual and tangent contributions at once.
   */
  template <class Derived>
  class NonLinearFormIntegrator
  {
    public:
      /**
       * @brief Returns both residual and tangent contributions as a ProblemBody.
       *
       * Equivalent to `Tangent(u, v) + Residual(v)`.
       */
      template <class Trial, class Test>
      auto operator()(const Trial& u, const Test& v) const
      {
        const Derived& self = static_cast<const Derived&>(*this);
        return self.Tangent(u, v) + self.Residual(v);
      }
  };
}

#endif
