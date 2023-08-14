/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_COMPOSITION_H
#define RODIN_VARIATIONAL_COMPOSITION_H

#include <functional>

#include "ForwardDecls.h"
#include "ScalarFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup CompositionSpecializations Composition Template Specializations
   * @brief Template specializations of the Composition class.
   * @see Composition
   */

  /**
   * @brief Composition of a scalar valued function on the real line and a
   * scalar function on the mesh.
   *
   * Represents the composition of two functions @f$ f : \mathbb{R}
   * \rightarrow \mathbb{R} @f$ and @f$ g : \Omega \rightarrow \mathbb{R} @f$:
   * @f[
   *   (f \circ g)(x) = f(g(x)) \ .
   * @f]
   */
  template <class NestedDerived>
  class Composition<std::function<Scalar(Scalar)>, FunctionBase<NestedDerived>> final
    : public ScalarFunctionBase<Composition<std::function<Scalar(Scalar)>, FunctionBase<NestedDerived>>>
  {
    public:
      using LHS = std::function<Scalar(Scalar)>;
      using RHS = FunctionBase<NestedDerived>;
      using Parent = ScalarFunctionBase<Composition<LHS, RHS>>;

      Composition(LHS f, const RHS& g)
        : m_f(f), m_g(g)
      {}

      Composition(const Composition& other)
        : Parent(other),
          m_f(other.m_f), m_g(other.m_g)
      {}

      Composition(Composition&& other)
        : m_f(std::move(other.m_f)), m_g(std::move(other.m_g))
      {}

      inline
      Scalar getValue(const Geometry::Point& p) const
      {
        return m_f(m_g->getValue(p));
      }

    private:
      LHS m_f;
      RHS m_g;
  };

  template <class NestedDerived>
  Composition(const std::function<Scalar(Scalar)>&, const FunctionBase<NestedDerived>&)
    -> Composition<std::function<Scalar(Scalar)>, FunctionBase<NestedDerived>>;

  /**
   * @brief Composes two functions
   * @param f
   * @param g
   *
   * Represents the composition of two functions @f$ f @f$ and @f$ g @f$:
   * @f[
   *   (f \circ g)(x) = f(g(x))
   * @f]
   *
   */
  template <class Lhs, class Rhs>
  auto compose(Lhs&& f, Rhs&& g)
  {
    return Composition(std::forward<Lhs>(f), std::forward<Rhs>(g));
  }
}

#endif
