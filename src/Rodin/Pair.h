/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_PAIR_H
#define RODIN_PAIR_H

#include "Tuple.h"

namespace Rodin
{
  template <class L, class R>
  class Pair : public Tuple<L, R>
  {
    public:
      using Parent = Tuple<L, R>;

      using Parent::Parent;

      Pair(const Pair& other)
        : Parent(other)
      {}

      Pair(Pair&& other)
        : Parent(std::move(other))
      {}

      inline
      constexpr
      Pair& operator=(const Pair& other)
      {
        Parent::operator=(other);
        return *this;
      }

      inline
      constexpr
      Pair& operator=(Pair&& other)
      {
        Parent::operator=(std::move(other));
        return *this;
      }

      auto& first()
      {
        return this->template get<0>();
      }

      auto& second()
      {
        return this->template get<1>();
      }

      const auto& first() const
      {
        return this->template get<0>();
      }

      const auto& second() const
      {
        return this->template get<1>();
      }
  };

  template <class L, class R>
  Pair(L, R) -> Pair<L, R>;
}

#endif
