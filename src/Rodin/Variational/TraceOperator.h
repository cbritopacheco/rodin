/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_COMPOSITION_H
#define RODIN_VARIATIONAL_COMPOSITION_H

/**
 * @file
 * @brief Trace operator for restricting functions to boundaries and interfaces.
 *
 * This file defines the TraceOperator which restricts functions from a domain
 * to lower-dimensional manifolds (boundaries, interfaces). For a function
 * @f$ u: \Omega \rightarrow \mathbb{R} @f$, the trace is:
 * @f[
 *   \gamma(u) = u|_{\Gamma}
 * @f]
 * where @f$ \Gamma \subset \partial\Omega @f$ is a boundary or interface.
 *
 * ## Trace Theorems
 * For functions in Sobolev spaces:
 * - If @f$ u \in H^1(\Omega) @f$, then @f$ \gamma(u) \in H^{1/2}(\partial\Omega) @f$
 * - The trace operator is bounded and surjective
 *
 * ## Applications
 * - Boundary conditions
 * - Interface continuity in DG methods
 * - Coupling bulk and surface equations
 */

#include <functional>

#include "ForwardDecls.h"
#include "RealFunction.h"

namespace Rodin::Variational
{
  /**
   * @defgroup TraceOperatorSpecializations TraceOperator Template Specializations
   * @brief Template specializations of the TraceOperator class.
   * @see TraceOperator
   */

  /**
   * @ingroup TraceOperatorSpecializations
   * @brief Trace operator
   */
  template <>
  class TraceOperator<FunctionBase> : public FunctionBase
  {
    public:
      TraceOperator(const FunctionBase& fn, Geometry::Attribute attr)
        : m_fn(fn.copy()),
          m_attr(attr)
      {}

      TraceOperator(const TraceOperator& other)
        :  FunctionBase(other),
          m_fn(other.m_fn->copy()),
          m_attr(other.m_attr)
      {}

      TraceOperator(TraceOperator&& other)
        :  FunctionBase(std::move(other)),
          m_fn(std::move(other.m_fn)),
          m_attr(other.m_attr)
      {}

    private:
      std::unique_ptr<FunctionBase> m_fn;
      Geometry::Attribute m_attr;
  };
  // TraceOperator(const FunctionBase&, Geometry::Attribute) -> TraceOperator<FunctionBase>;
}

#endif

