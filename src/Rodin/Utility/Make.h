/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_MAKE_H
#define RODIN_UTILITY_MAKE_H

/**
 * @file
 * @brief Defines the Make factory template for perfect forwarding construction.
 */

#include <functional>

namespace Rodin::Utility
{
  /**
   * @brief Factory template for perfect forwarding construction.
   * @ingroup UtilityModule
   * @tparam T The type to construct.
   *
   * The Make template provides a function object that constructs instances
   * of type T using perfect forwarding to preserve value categories of arguments.
   * This is particularly useful in template metaprogramming where construction
   * needs to be performed within template contexts.
   *
   * The Make template serves as a type-safe factory that ensures proper
   * forwarding semantics and can be used in situations where direct construction
   * is not possible or convenient.
   *
   * Example usage:
   * @code{.cpp}
   * auto factory = Make<std::vector<int>>{};
   * auto vec = factory(10, 42);  // Creates std::vector<int>(10, 42)
   * 
   * // Used in template contexts:
   * template<typename T, typename... Args>
   * auto createInstance(Args&&... args) {
   *     return Make<T>{}(std::forward<Args>(args)...);
   * }
   * @endcode
   */
  template <class T>
  struct Make
  {
    /**
     * @brief Constructs an instance of T using perfect forwarding.
     * @tparam Params Parameter pack of constructor argument types.
     * @param params Constructor arguments to forward to T's constructor.
     * @return A new instance of T constructed with the forwarded parameters.
     *
     * This operator uses perfect forwarding to preserve the value categories
     * of the arguments, ensuring optimal performance and correct semantics
     * for both lvalue and rvalue arguments.
     */
    template <class ... Params>
    T operator()(Params&&... params)
    {
      return T{std::forward<Params>(params)...};
    }
  };
}

#endif
