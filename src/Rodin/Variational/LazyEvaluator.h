/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_VARIATIONAL_LAZYEVALUATOR_H
#define RODIN_VARIATIONAL_LAZYEVALUATOR_H

#include <optional>

#include "Function.h"
#include "RangeShape.h"
#include "ForwardDecls.h"

namespace Rodin::Variational
{
  template <class FunctionDerived>
  class LazyEvaluator<FunctionBase<FunctionDerived>> final
    : public FunctionBase<LazyEvaluator<FunctionBase<FunctionDerived>>>
  {
    public:
      using Parent = FunctionBase<LazyEvaluator<FunctionBase<FunctionDerived>>>;

      constexpr
      LazyEvaluator(const FunctionBase<FunctionDerived>& ref)
        : m_ref(ref)
      {}

      constexpr
      LazyEvaluator(const LazyEvaluator& other)
        : Parent(other),
          m_ref(other.m_ref)
      {}

      constexpr
      LazyEvaluator(LazyEvaluator&& other)
        : Parent(std::move(other)),
          m_ref(std::move(other.m_ref))
      {}

      inline
      constexpr
      RangeShape getRangeShape() const
      {
        return m_ref.get().getRangeShape();
      }

      inline
      constexpr
      auto getValue(const Geometry::Point& p) const
      {
        return m_ref.get().getValue(p);
      }

      inline LazyEvaluator* copy() const noexcept override
      {
        return new LazyEvaluator(*this);
      }

    private:
      std::reference_wrapper<const FunctionBase<FunctionDerived>> m_ref;
  };

  template <class FunctionDerived>
  LazyEvaluator(const FunctionBase<FunctionDerived>&)
    -> LazyEvaluator<FunctionBase<FunctionDerived>>;
}

#endif
