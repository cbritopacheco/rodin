/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_OVERLOADED_H
#define RODIN_UTILITY_OVERLOADED_H

/**
 * @file
 * @brief Defines the Overloaded helper for combining multiple lambdas.
 */

namespace Rodin::Utility
{
  /**
   * @brief Helper type for creating overloaded function objects from lambdas.
   * @ingroup UtilityModule
   * @tparam Ts Variadic template pack of callable types.
   *
   * The Overloaded struct is a utility for combining multiple lambda functions
   * or callable objects into a single overloaded function object. This is
   * particularly useful with std::visit and variant types where different
   * alternatives require different handling logic.
   *
   * The struct inherits from all provided callable types and brings their
   * operator() functions into scope using the 'using' declaration, creating
   * an overload set that can be used for pattern matching on variants.
   *
   * Example usage:
   * @code
   * std::variant<int, std::string, double> var = 42;
   * 
   * auto visitor = Overloaded {
   *     [](int i) { return "integer: " + std::to_string(i); },
   *     [](const std::string& s) { return "string: " + s; },
   *     [](double d) { return "double: " + std::to_string(d); }
   * };
   * 
   * std::string result = std::visit(visitor, var);
   * @endcode
   */
  template <class... Ts>
  struct Overloaded : Ts...
  {
    /// @brief Brings all operator() functions from base types into scope.
    using Ts::operator()...;
  };

  /**
   * @brief Deduction guide for Overloaded.
   * @ingroup UtilityModule
   * @tparam Ts The types of the callable objects.
   *
   * This deduction guide allows the compiler to automatically deduce the
   * template parameters when constructing an Overloaded object from a set
   * of callable objects.
   */
  template<class... Ts> Overloaded(Ts...) -> Overloaded<Ts...>;
}

#endif
