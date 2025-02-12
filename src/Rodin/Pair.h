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
  /**
   * @brief A Pair class representing a tuple with two elements.
   *
   * The Pair class extends the Tuple class for exactly two elements, providing
   * convenient accessors for the first and second elements.
   *
   * @tparam L Type of the first element.
   * @tparam R Type of the second element.
   */
  template <class L, class R>
  class Pair : public Tuple<L, R>
  {
    public:
      /// @brief Parent class type.
      using Parent = Tuple<L, R>;

      /// @brief Inherit constructors from the parent Tuple.
      using Parent::Parent;

      /**
       * @brief Copy constructor.
       *
       * @param other The Pair object to copy from.
       */
      Pair(const Pair& other)
        : Parent(other)
      {}

      /**
       * @brief Move constructor.
       *
       * @param other The Pair object to move from.
       */
      Pair(Pair&& other)
        : Parent(std::move(other))
      {}

      /**
       * @brief Copy assignment operator.
       *
       * @param other The Pair object to copy from.
       * @return Reference to this Pair.
       */
      constexpr
      Pair& operator=(const Pair& other)
      {
        Parent::operator=(other);
        return *this;
      }

      /**
       * @brief Move assignment operator.
       *
       * @param other The Pair object to move from.
       * @return Reference to this Pair.
       */
      constexpr
      Pair& operator=(Pair&& other)
      {
        Parent::operator=(std::move(other));
        return *this;
      }

      /**
       * @brief Retrieves the first element of the Pair.
       *
       * @return Reference to the first element.
       */
      auto& first()
      {
        return this->template get<0>();
      }

      /**
       * @brief Retrieves the second element of the Pair.
       *
       * @return Reference to the second element.
       */
      auto& second()
      {
        return this->template get<1>();
      }

      /**
       * @brief Retrieves the first element of the Pair (const version).
       *
       * @return Const reference to the first element.
       */
      const auto& first() const
      {
        return this->template get<0>();
      }

      /**
       * @brief Retrieves the second element of the Pair (const version).
       *
       * @return Const reference to the second element.
       */
      const auto& second() const
      {
        return this->template get<1>();
      }
  };

  /**
   * @brief Deduction guide for the Pair class.
   *
   * This guide allows the compiler to deduce the template arguments for L and R
   * when constructing a Pair object.
   *
   * @tparam L Type of the first element.
   * @tparam R Type of the second element.
   */
  template <class L, class R>
  Pair(L, R) -> Pair<L, R>;
}

#endif
