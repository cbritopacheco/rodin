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
  template <class StrictType>
  class LazyEvaluator : public FunctionBase<LazyEvaluator<StrictType>>
  {
    public:
      using Parent = FunctionBase<LazyEvaluator<StrictType>>;

      /**
       * @brief R-Values are not allowed.
       */
      LazyEvaluator(StrictType&&) = delete;

      /**
       * @brief Prevent implicit copies.
       */
      LazyEvaluator(const StrictType& ref) = delete;

      /**
       * @brief Constructs the LazyEvaluator object from a constant reference
       * the data-full object.
       */
      explicit
      constexpr
      LazyEvaluator(std::reference_wrapper<const StrictType> ref)
        : m_ref(ref)
      {}

      /**
       * @brief Copy constructor.
       */
      constexpr
      LazyEvaluator(const LazyEvaluator& other)
        : Parent(other),
          m_ref(other.m_ref)
      {}

      /**
       * @brief Move constructor.
       */
      constexpr
      LazyEvaluator(LazyEvaluator&& other)
        : Parent(std::move(other)),
          m_ref(std::move(other.m_ref))
      {}

      LazyEvaluator& operator=(const LazyEvaluator&) = delete;

      LazyEvaluator& operator=(LazyEvaluator&&) = delete;

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

      inline LazyEvaluator* copy() const noexcept final override
      {
        return new LazyEvaluator(*this);
      }

    private:
      std::reference_wrapper<const StrictType> m_ref;
  };

  /**
   * @brief CTAD for LazyEvaluator.
   */
  template <class StrictType>
  LazyEvaluator(std::reference_wrapper<const StrictType>) -> LazyEvaluator<StrictType>;
}

#endif
