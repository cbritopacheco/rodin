/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
/**
 * @file TraceOperator.h
 * @brief Function trace operator for boundary restriction.
 *
 * This file defines the TraceOperator class, which computes the trace (boundary
 * restriction) of functions to specified mesh boundaries. This is different from
 * the matrix trace in Trace.h - this operator restricts function domains.
 *
 * ## Mathematical Foundation
 * For a function @f$ u : \Omega \to \mathbb{R} @f$ and boundary @f$ \Gamma \subset \partial\Omega @f$,
 * the trace operator computes:
 * @f[
 *   \gamma_\Gamma(u) = u|_\Gamma
 * @f]
 * which is the restriction of @f$ u @f$ to the boundary @f$ \Gamma @f$.
 *
 * ## Trace Theorem
 * For @f$ u \in H^1(\Omega) @f$, the trace @f$ u|_{\partial\Omega} @f$ is well-defined
 * and belongs to @f$ H^{1/2}(\partial\Omega) @f$.
 *
 * ## Applications
 * - Boundary condition evaluation
 * - Surface integrals
 * - Coupling interface conditions
 * - Weak boundary term evaluation
 *
 * ## Usage Example
 * ```cpp
 * // Evaluate function on specific boundary
 * auto u_boundary = TraceOperator(u, boundary_attr);
 * 
 * // Use in boundary integral
 * auto bc = BoundaryIntegral(u_boundary, v).on(boundary_attr);
 * ```
 */
#ifndef RODIN_VARIATIONAL_TRACEOPERATOR_H
#define RODIN_VARIATIONAL_TRACEOPERATOR_H

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
   * @brief Trace (boundary restriction) operator for functions.
   *
   * TraceOperator restricts a function to a specified boundary region,
   * enabling evaluation and integration on boundary surfaces.
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

