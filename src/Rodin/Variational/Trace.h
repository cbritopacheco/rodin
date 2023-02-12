/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_TRACE_H
#define RODIN_VARIATIONAL_TRACE_H

#include "ForwardDecls.h"

#include "ScalarFunction.h"
#include "BasisOperator.h"

namespace Rodin::Variational
{
  /**
   * @defgroup TraceSpecializations
   * @brief Template specializations of the Trace class.
   * @see Trace
   */

  /**
   * @ingroup TraceSpecializations
   * @brief Trace of a FunctionBase instance.
   */
  template <>
  class Trace<FunctionBase> : public ScalarFunctionBase
  {
    public:
      /**
       * @brief Constructs the Trace of the given matrix
       * @param[in] m Square matrix
       */
      Trace(const FunctionBase& m);

      Trace(const Trace& other);

      Trace(Trace&& other);

      FunctionValue getValue(const Geometry::Point& p) const override
      {
        return m_matrix->getValue(p).matrix().trace();
      }

      Trace& traceOf(Geometry::Attribute attrs) override
      {
        ScalarFunctionBase::traceOf(attrs);
        m_matrix->traceOf(attrs);
        return *this;
      }

      Trace* copy() const noexcept override
      {
        return new Trace(*this);
      }
    private:
      std::unique_ptr<FunctionBase> m_matrix;
  };
  Trace(const FunctionBase&) -> Trace<FunctionBase>;
}

#endif
