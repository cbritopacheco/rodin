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

