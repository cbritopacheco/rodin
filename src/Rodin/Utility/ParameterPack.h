/*
 *          Copyright Carlos BRITO PACHECO 2021 - 2022.
 * Distributed under the Boost Software License, Version 1.0.
 *       (See accompanying file LICENSE or copy at
 *          https://www.boost.org/LICENSE_1_0.txt)
 */
#ifndef RODIN_UTILITY_PARAMETERPACK_H
#define RODIN_UTILITY_PARAMETERPACK_H

/**
 * @file
 * @brief Defines the ParameterPack class for compile-time parameter pack manipulation.
 */

#include <type_traits>

namespace Rodin::Utility
{
  namespace Internal
  {
    /**
     * @brief Internal implementation for accessing parameter pack elements by index.
     * @ingroup UtilityModule
     * @tparam N The index to access.
     * @tparam Types The parameter pack of types.
     *
     * This internal template recursively traverses the parameter pack to
     * extract the type at the specified index.
     */
    template <std::size_t N, class ... Types>
    struct AtImpl;

    /**
     * @brief Base case specialization for index 0.
     * @ingroup UtilityModule
     * @tparam First The first type in the parameter pack.
     * @tparam Rest The remaining types in the parameter pack.
     *
     * When the index reaches 0, this specialization returns the First type.
     */
    template <class First, class ... Rest>
    struct AtImpl<0, First, Rest...>
    {
      /// @brief The type at index 0.
      using Type = First;
    };

    /**
     * @brief Recursive case for index N > 0.
     * @ingroup UtilityModule
     * @tparam N The index to access (N > 0).
     * @tparam First The first type in the parameter pack.
     * @tparam Rest The remaining types in the parameter pack.
     *
     * This specialization recursively decrements the index and removes
     * the first type from the pack until reaching index 0.
     */
    template <std::size_t N, class First, class ... Rest>
    struct AtImpl<N, First, Rest...>
    {
      /// @brief The type at index N, obtained recursively.
      using Type = typename AtImpl<N - 1, Rest...>::Type;
    };

    /**
     * @brief Internal implementation for checking if all types satisfy a predicate.
     * @ingroup UtilityModule
     * @tparam Predicate A template predicate that evaluates types.
     * @tparam Types The parameter pack of types to check.
     *
     * This internal template recursively applies the predicate to all types
     * in the parameter pack and performs logical AND of the results.
     */
    template <template <class> class Predicate, class ...>
    struct AllImpl;

    /**
     * @brief Base case for a single type.
     * @ingroup UtilityModule
     * @tparam Predicate The predicate to apply.
     * @tparam T1 The single type to check.
     */
    template <template <class> class Predicate, class T1>
    struct AllImpl<Predicate, T1>
    {
      /// @brief Result of applying the predicate to T1.
      static constexpr bool Value = Predicate<T1>::Value;
    };

    /**
     * @brief Recursive case for multiple types.
     * @ingroup UtilityModule
     * @tparam Predicate The predicate to apply.
     * @tparam T1 The first type to check.
     * @tparam T2 The second type to check.
     * @tparam Ts The remaining types to check.
     */
    template <template <class> class Predicate, class T1, class T2, class ... Ts>
    struct AllImpl<Predicate, T1, T2, Ts...>
    {
      /// @brief Logical AND of predicate results for all types.
      static constexpr bool Value = Predicate<T1>::Value && AllImpl<Predicate, T2, Ts...>::Value;
    };
  }

  /**
   * @brief Template metaprogramming utilities for parameter packs.
   * @ingroup UtilityModule
   * @tparam Params The parameter pack of types.
   *
   * The ParameterPack class provides compile-time utilities for manipulating
   * and querying template parameter packs. It offers type-safe access to
   * individual elements by index and predicate-based validation of all types.
   *
   * This utility is essential for template metaprogramming where compile-time
   * type manipulation and validation are required.
   *
   * Example usage:
   * @code{.cpp}
   * using Pack = ParameterPack<int, double, std::string>;
   * 
   * // Get type at index 1 (double)
   * using SecondType = Pack::At<1>;
   * 
   * // Check if all types are arithmetic
   * constexpr bool allArithmetic = Pack::All<std::is_arithmetic>::Value;
   * 
   * // Get pack size
   * constexpr size_t packSize = Pack::Size;  // 3
   * @endcode
   */
  template <class ... Params>
  class ParameterPack
  {
    public:
      /**
       * @brief Type alias for accessing the type at a specific index.
       * @tparam Index The zero-based index of the type to access.
       *
       * This template alias provides compile-time access to the type at the
       * specified index within the parameter pack. Index bounds are not
       * checked at compile time, so out-of-bounds access will result in
       * compilation errors.
       */
      template <std::size_t Index>
      using At = typename Internal::AtImpl<Index, Params...>::Type;

      /**
       * @brief Template alias for checking if all types satisfy a predicate.
       * @tparam Predicate A template predicate that should have a static constexpr bool Value member.
       *
       * This template alias applies the given predicate to all types in the
       * parameter pack and provides a Value member indicating whether all
       * types satisfy the predicate condition.
       */
      template <template <class> class Predicate>
      using All = typename Internal::AllImpl<Predicate, Params...>;

      /// @brief The number of types in the parameter pack.
      static constexpr size_t Size = sizeof...(Params);
  };
}

#endif
